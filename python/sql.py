import sqlite3
from pathlib import Path

from .ref import CustomTypes as CT
from .orc import parse_log_file


def amchi_in_database(amchi: str, connection: sqlite3.Connection):
    """Returns True if amchi string is present in stationaries table of database."""
    cursor = connection.cursor()
    cursor.execute("SELECT id FROM stationary WHERE amchi = ?", (amchi,))
    matches = cursor.fetchall()
    return matches is None


def connect(data_directory: str | Path) -> sqlite3.Connection:
    """Connects to a sql database."""
    Path(data_directory).mkdir(exist_ok=True, parents=True)
    return sqlite3.connect(data_directory / "data.db")


def enumerated_graph_into_database(
    enumerated_graph: CT.NetworkXGraph, data_directory: str | Path, connection: sqlite3.Connection
):
    """Fills sqlite3 database with enumerated reaction graph."""
    cursor = connection.cursor()

    to_submit = {}
    amchi_to_ids = {}
    for amchi, data in enumerated_graph.nodes(data=True):
        if data.get("role") not in ("reactant", "product"):
            continue

        smiles = data.get("smiles")
        xyz = data.get("xyz")

        cursor.execute("SELECT id FROM stationary WHERE amchi = ?", (amchi,))
        result = cursor.fetchone()

        if not result:
            cursor.execute(
                """
                INSERT INTO stationary (amchi, smiles, xyz)
                VALUES (?, ?, ?)
                """,
                (amchi, smiles, xyz),
            )
            amchi_to_ids[amchi] = cursor.lastrowid

            directory = data_directory / amchi
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "guess.xyz").write_text(xyz + "\n")

            to_submit[amchi] = "stationary"

        else:
            amchi_to_ids[amchi] = result[0]

    for amchi, data in enumerated_graph.nodes(data=True):
        if data.get("role") not in ("transition"):
            continue

        cursor.execute("SELECT id FROM transition WHERE amchi = ?", (amchi,))
        result = cursor.fetchone()

        if not result:
            xyz = data.get("xyz")
            scan = data.get("scan")
            scan_string = "\n".join(scan)
            reactants = data.get("reactants")
            products = data.get("products")

            r1 = amchi_to_ids[reactants[0]]
            r2 = amchi_to_ids[reactants[1]] if len(reactants) > 1 else None
            p1 = amchi_to_ids[products[0]]
            p2 = amchi_to_ids[products[1]] if len(products) > 1 else None

            cursor.execute(
                """
                INSERT INTO transition (amchi, reactant_1, reactant_2, product_1, product_2, xyz, scan)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (amchi, r1, r2, p1, p2, xyz, scan_string),
            )

            directory = data_directory / amchi
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "guess.xyz").write_text(xyz + "\n")

            to_submit[amchi] = "transition"

    connection.commit()

    return to_submit


def initialize_database(connection: sqlite3.Connection):
    """Initializes a database for containing data relevant to this project."""
    cursor = connection.cursor()

    cursor.execute("""
                CREATE TABLE IF NOT EXISTS stationary (
                        id INTEGER PRIMARY KEY,
                        amchi TEXT NOT NULL,
                        smiles TEXT NOT NULL,
                        xyz TEXT NOT NULL
                )
                """)

    cursor.execute("""
                CREATE TABLE IF NOT EXISTS transition (
                        id INTEGER PRIMARY KEY,
                        amchi TEXT NOT NULL,
                        reactant_1 INTEGER REFERENCES stationary(id),
                        reactant_2 INTEGER REFERENCES stationary(id),
                        product_1 INTEGER REFERENCES stationary(id),
                        product_2 INTEGER REFERENCES stationary(id),
                        xyz TEXT NOT NULL,
                        scan TEXT NOT NULL
                )
                """)

    cursor.execute("""
                CREATE TABLE IF NOT EXISTS energies (
                    id INTEGER PRIMARY KEY,
                    stationary_id INT REFERENCES stationary(id),
                    transition_id INT REFERENCES transition(id),
                    single_point REAL,
                    zero_point REAL,
                    total_energy REAL,
                    imaginary_frequency REAL,

                    CHECK (
                        (stationary_id is NOT NULL AND transition_id IS NULL)
                        OR
                        (transition_id is NOT NULL AND stationary_id IS NULL)
                    )
                )
                """)

def log_energies(connection: sqlite3.Connection, data_dir: str | Path):
    """Logs energies missing from database."""

    def _parse_energies(amchi: str):
        spc_line = parse_log_file(data_dir / amchi / "calc.log", "FINAL SINGLE POINT ENERGY")
        spc_energy = float(spc_line.split(" ")[-1]) * 627.5095 if spc_line else None

        zpv_line = parse_log_file(data_dir / amchi / "freq.log", "Zero point energy")
        zpv_energy = float(zpv_line.split(" ")[-2]) if zpv_line else None

        tot_energy = None if spc_energy is None or zpv_energy is None else spc_energy + zpv_energy

        return spc_energy, zpv_energy, tot_energy

    def _parse_imaginary(amchi: str):
        imaginary_line = parse_log_file(data_dir / amchi / "freq.log", "***imaginary mode***")
        
        return imaginary_line

    cursor = connection.cursor()

    cursor.execute("SELECT id, amchi FROM stationary")
    stationary_rows = cursor.fetchall()

    for stationary_id, amchi in stationary_rows:
        cursor.execute("""
        SELECT single_point, zero_point, total_energy
        FROM energies
        WHERE stationary_id = ?
        """, (stationary_id,))

        row = cursor.fetchone()
        if row is not None:
            spc_energy, zpv_energy, tot_energy = row

            if spc_energy is None or zpv_energy is None or tot_energy is None:
                spc_energy, zpv_energy, tot_energy = _parse_energies(amchi)

                cursor.execute("""
                UPDATE energies
                SET single_point = ?, zero_point = ?, total_energy = ?
                WHERE stationary_id = ?
                """, (spc_energy, zpv_energy, tot_energy, stationary_id))

        else:
            spc_energy, zpv_energy, tot_energy = _parse_energies(amchi)
                
            cursor.execute("""
            INSERT INTO energies (stationary_id, single_point, zero_point, total_energy)
            VALUES (?, ?, ?, ?)
            """, (stationary_id, spc_energy, zpv_energy, tot_energy))

    cursor.execute("SELECT id, amchi FROM transition")
    transition_rows = cursor.fetchall()

    for transition_id, amchi in transition_rows:
        cursor.execute("""
        SELECT single_point, zero_point, total_energy, imaginary_frequency
        FROM energies
        WHERE transition_id = ?
        """, (transition_id,))

        row = cursor.fetchone()
        if row is not None:
            spc_energy, zpv_energy, tot_energy, imaginary_frequency = row

            if spc_energy is None or zpv_energy is None or tot_energy is None:
                spc_energy, zpv_energy, tot_energy = _parse_energies(amchi)

                cursor.execute("""
                UPDATE energies
                SET single_point = ?, zero_point = ?, total_energy = ?
                WHERE transition_id = ?
                """, (spc_energy, zpv_energy, tot_energy, transition_id))

            if imaginary_frequency is None:
                imag_line = _parse_imaginary(amchi).strip()
                imag_freq = imag_line.split(" ")[3]
                cursor.execute("""
                UPDATE energies
                SET imaginary_frequency = ?
                WHERE transition_id = ?
                """, (imag_freq, transition_id))

        else:
            spc_energy, zpv_energy, tot_energy = _parse_energies(amchi)
            imag_line = _parse_imaginary(amchi).strip()
            imag_freq = imag_line.split(" ")[3]
                
            cursor.execute("""
            INSERT INTO energies (transition_id, single_point, zero_point, total_energy, imaginary_frequency)
            VALUES (?, ?, ?, ?, ?)
            """, (transition_id, spc_energy, zpv_energy, tot_energy, imag_freq))

        connection.commit()

            
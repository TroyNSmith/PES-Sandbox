import sqlite3
from pathlib import Path

from .ref import CustomTypes as CT


def amchi_in_database(amchi: str, connection: sqlite3.Connection):
    """Returns True if amchi string is present in stationaries table of database."""
    cursor = connection.cursor()
    cursor.execute("SELECT id FROM stationary WHERE amchi = ?", (amchi,))
    matches = cursor.fetchall()
    return matches is None


def connect(file_path: str | Path) -> sqlite3.Connection:
    """Connects to a sql database."""
    return sqlite3.connect(file_path)


def enumerated_graph_into_database(
    enumerated_graph: CT.NetworkXGraph, database_path: str | Path
):
    """Fills sqlite3 database with enumerated reaction graph."""
    data_directory = Path(database_path).parent
    connection = sqlite3.connect(database_path)
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
    connection.close()

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
                    zero_point REAL
               
                    CHECK (
                        (stationary_id is NOT NULL AND transition_id IS NULL)
                        OR
                        (transition_id is NOT NULL AND stationary_id IS NULL)
                    )
                )
                """)

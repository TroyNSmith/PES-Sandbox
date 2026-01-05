import sqlite3
from pathlib import Path

from .ref import CustomTypes as CT


def connect(file_path: str | Path) -> sqlite3.Connection:
    """Connects to a sql database."""
    return sqlite3.connect(file_path)


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


def enumerated_graph_into_database(
    enumerated_graph: CT.NetworkXGraph, connection: sqlite3.Connection
):
    """Fills sqlite3 database with enumerated reaction graph."""
    cursor = connection.cursor()

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

    connection.commit()

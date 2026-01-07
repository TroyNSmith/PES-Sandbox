import networkx
from rdkit.Chem import Mol


class CustomTypes:
    """Collects custom type objects for Python scripts."""

    RDMols = Mol | tuple[Mol, ...]

    SMILES = str
    SMILES_Set = tuple[SMILES, ...]

    InChI = str
    InChI_Set = tuple[InChI, ...]

    NetworkXGraph = networkx.Graph
    AutomolGraph = tuple[dict, dict]


class ValenceIdentities:
    """Collects identifiers for valence classes."""

    A_ = "Cv4,Ov2"
    C_ = "Cv4"
    O_ = "Ov2"
    As = "CX4,OX2"
    Cs = "CX4"
    Os = "OX2"
    Ar = "Cv3,Ov1"
    Cr = "Cv3"
    Or = "Ov1"
    Au = "CX3,OX1"
    Cu = "CX3"
    Ou = "OX1"

    Radical_Tokens = [Ar, Cr, Or]

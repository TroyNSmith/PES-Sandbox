import re
from collections.abc import Sequence

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdchem import Mol

from . import ntx
from .ref import CustomTypes as CT
from .ref import ValenceIdentities as VI


def isomorphic(molecules_1: CT.RDMols, molecules_2: CT.RDMols) -> bool:
    """Determines whether two sets of molecules are isomorphic to each other."""
    matcher = ntx.graph_matcher(molecules_1, molecules_2)
    return matcher.is_isomorphic()


def mol_from_smiles(smiles: str | tuple[str], with_coords: bool = False) -> CT.RDMols:
    """Generates RDKit molecule(s) from SMILES string(s)."""
    smiles = (smiles,) if isinstance(smiles, str) else smiles

    mols = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)
        if with_coords:
            raise NotImplementedError
        mols.append(mol)

    return tuple(mols)


def mol_to_smiles(
    molecules: Mol | CT.RDMols, ignore_map_numbers: bool = True
) -> CT.SMILES_Set:
    """Converts RDKit molecule(s) to SMILES string(s)."""
    molecules = (molecules,) if isinstance(molecules, Mol) else molecules

    smiles = []
    for molecule in molecules:
        if ignore_map_numbers:
            molecule = Chem.RemoveHs(molecule)
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)

            Chem.rdmolops.AssignRadicals(molecule)

        smiles.append(Chem.MolToSmiles(molecule, ignoreAtomMapNumbers=True))

    return tuple(smiles)


def mol_to_inchi(molecules: Mol | CT.RDMols) -> CT.InChI_Set:
    """Converts RDKit molecule(s) to standard InChI string(s)."""
    molecules = (molecules,) if isinstance(molecules, Mol) else molecules

    inchis = []
    for molecule in molecules:
        inchis.append(Chem.MolToInchi(molecule))

    return tuple(inchis)


def radicals_from_smarts(reaction_smarts: str) -> set[int]:
    """Infers radical map numbers from reaction smarts."""
    radicals = set()

    atom_re = re.compile(r"\[([^\]:]+):(\d+)\]")
    for valence_identity, map_number in atom_re.findall(reaction_smarts):
        if valence_identity in VI.Radical_Tokens:
            radicals.add(int(map_number))

    return radicals


def run_reaction(
    reactants: CT.RDMols, reaction_smarts: str, isomorphs: bool = False
) -> list[CT.RDMols]:
    """Returns Mol products given reactant Mol(s) and reaction SMARTS."""
    rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts)

    _, rhs_smarts = reaction_smarts.split(">>")
    rhs_radicals = radicals_from_smarts(rhs_smarts)

    if isinstance(reactants, Mol):
        reactants = (reactants,)

    for reactant in reactants:
        for atom in reactant.GetAtoms():
            atom.SetIntProp("molAtomMapNumber", atom.GetIdx())

    product_sets = list(rxn.RunReactants(reactants))
    product_sets = product_sets if isomorphs else unique_molecules(product_sets)

    for products in product_sets:
        for product in products:
            for atom in product.GetAtoms():
                if atom.HasProp("old_mapno"):
                    if int(atom.GetProp("old_mapno")) in rhs_radicals:
                        atom.SetNumRadicalElectrons(1)

    return product_sets


def unique_molecules(molecules: Sequence[CT.RDMols]) -> CT.RDMols:
    "Identifies unique molecules from a sequence."
    unique_set = []

    for molecule in molecules:
        if not any(isomorphic(molecule, m) for m in unique_set):
            unique_set.append(molecule)

    return unique_set


class Reaction_Templates:
    PROTON_TRANSFER = (
        f"([{VI.Ar}:1].[{VI.As}:2]([H:3]))>>([{VI.As}:1]([H:3]).[{VI.Ar}:2])"
    )
    RING_OPENING = (
        f"[{VI.A_}:1]-@[{VI.A_}:2]-[{VI.Ar}:3]>>([{VI.Ar}:1].[{VI.A_}:2]=[{VI.A_}:3])"
    )

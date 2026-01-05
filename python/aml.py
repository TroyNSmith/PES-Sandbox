import automol
import networkx as nx
from rdkit.Chem.rdchem import Mol

from .rdk import mol_to_smiles
from .ref import CustomTypes as CT


def molecular_graph(smiles: str, canonical: str = True) -> CT.AutomolGraph:
    """Generates a molecular AutoMol graph from SMILES string."""
    graph = automol.smiles.graph(smiles)

    if canonical:
        return automol.graph.canonical(graph)

    return graph


def process_rdkit_reaction(reactants: Mol | CT.RDMols, product_sets: list[CT.RDMols]):
    """Processes RDKit reactions through the AutoMol package."""
    enumerated_graph = nx.DiGraph()

    def _add_stationary(molecules, role) -> list:
        amchis_list, smiles_list = [], []

        for molecule in molecules:
            (smiles,) = mol_to_smiles(molecule)
            graph = molecular_graph(smiles)

            amchi = automol.graph.amchi(graph)
            geom = automol.graph.geometry(graph)
            xyz = automol.geom.xyz_string(geom)

            smiles_list.append(smiles)
            amchis_list.append(amchi)

            enumerated_graph.add_node(amchi, smiles=smiles, xyz=xyz, role=role)

        return amchis_list, smiles_list

    (reactant_amchis, reactant_smiles) = _add_stationary(reactants, "reactant")

    for products in product_sets:
        (product_amchis, product_smiles) = _add_stationary(products, "product")

        (reaction,) = reaction_graph(tuple(reactant_smiles), tuple(product_smiles))

        transition_graph = automol.reac.ts_graph(reaction)
        transition_graph = automol.graph.canonical(transition_graph)

        transition_amchi = automol.graph.amchi(transition_graph)
        transition_geom = automol.graph.geometry(transition_graph)
        transition_xyz = automol.geom.xyz_string(transition_geom)

        (_, bonds) = transition_graph
        broken, formed = [], []
        for bond, (order, _) in bonds.items():
            if order == 0.1:
                formed.append(bond)
            elif order == 0.9:
                broken.append(bond)

        dmat_angstrom = automol.geom.distance_matrix(transition_geom) * 0.52917721092

        scans = []
        for broken_bond in broken:
            a, b = broken_bond
            if len(formed) > 0:
                for formed_bond in formed:
                    shared = broken_bond & formed_bond
                    if len(shared) == 1:
                        c = next(iter(formed_bond - shared))
                        dist = dmat_angstrom[b, c]
                        scans.append(f"B {b} {c} = {dist:.3f}, 0.7, 100")
            else:
                scans.append(f"B {a} {b} = {dmat_angstrom[a, b]:.3f}, 2.0, 100")

        enumerated_graph.add_node(
            transition_amchi,
            xyz=transition_xyz,
            scan=scans,
            role="transition",
            reactants=reactant_amchis,
            products=product_amchis,
        )

        for r_amchi in reactant_amchis:
            enumerated_graph.add_edge(r_amchi, transition_amchi, kind="reactant_to_ts")

        for p_amchi in product_amchis:
            enumerated_graph.add_edge(transition_amchi, p_amchi, kind="ts_to_product")

    return enumerated_graph


def reaction_graph(
    reactants: CT.SMILES | CT.SMILES_Set,
    products: CT.SMILES | CT.SMILES_Set,
    include_stereo: bool = False,
    canonical: bool = True,
) -> CT.AutomolGraph:
    """Generates a reaction Automol graph from SMILES strings."""
    reactants = (reactants,) if isinstance(reactants, str) else reactants
    products = (products,) if isinstance(products, str) else products

    graph = automol.reac.from_smiles(reactants, products, include_stereo)

    return graph

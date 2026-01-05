import automol
import networkx as nx
from rdkit.Chem.rdchem import Mol

from .rdk import mol_to_smiles
from .ref import CustomTypes as CT


def molecular_graph(smiles: str) -> CT.AutomolGraph:
    """Generates a molecular AutoMol graph from SMILES string."""
    graph = automol.smiles.graph(smiles)

    return graph


def process_rdkit_reaction(reactants: Mol | CT.RDMols, product_sets: list[CT.RDMols]):
    """Processes RDKit reactions through the AutoMol package."""
    enumerated_graph = nx.DiGraph()

    reactant_amchis = []
    reactant_smiles_list = []

    for reactant in reactants:
        (reactant_smiles,) = mol_to_smiles(reactant)
        reactant_smiles_list.append(reactant_smiles)

        reactant_graph = molecular_graph(reactant_smiles)
        reactant_graph = automol.graph.canonical(reactant_graph)

        reactant_amchi = automol.graph.amchi(reactant_graph)
        reactant_amchis.append(reactant_amchi)

        reactant_geom = automol.graph.geometry(reactant_graph)
        reactant_xyz = automol.geom.xyz_string(reactant_geom)

        enumerated_graph.add_node(
            reactant_amchi, smiles=reactant_smiles, xyz=reactant_xyz, role="reactant"
        )

    for products in product_sets:
        product_amchis = []
        product_smiles_list = []

        for product in products:
            (product_smiles,) = mol_to_smiles(product)
            product_smiles_list.append(product_smiles)

            product_graph = molecular_graph(product_smiles)
            product_graph = automol.graph.canonical(product_graph)

            product_amchi = automol.graph.amchi(product_graph)
            product_amchis.append(product_amchi)

            product_geom = automol.graph.geometry(product_graph)
            product_xyz = automol.geom.xyz_string(product_geom)

            enumerated_graph.add_node(
                product_amchi, smiles=product_smiles, xyz=product_xyz, role="product"
            )

        (reaction,) = reaction_graph(
            tuple(reactant_smiles_list), tuple(product_smiles_list)
        )

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
) -> CT.AutomolGraph:
    """Generates a reaction Automol graph from SMILES strings."""
    reactants = (reactants,) if isinstance(reactants, str) else reactants
    products = (products,) if isinstance(products, str) else products

    graph = automol.reac.from_smiles(reactants, products, include_stereo)

    return graph

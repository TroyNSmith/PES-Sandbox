import automol
import networkx as nx
from rdkit.Chem.rdchem import Mol

from .rdk import mol_to_smiles
from .ref import CustomTypes as CT


def canonical_enantiomer(graph: CT.AutomolGraph):
    """Returns canonical amchi and graph of stationary species."""
    amchi = automol.graph.amchi(graph)

    if not automol.amchi.is_canonical_enantiomer(amchi):
        amchi = automol.amchi.canonical_enantiomer(amchi)
        graph = stationary_graph(amchi=amchi)

    return amchi, graph


def process_rdkit_reaction(reactants: Mol | CT.RDMols, product_sets: list[CT.RDMols]):
    """Processes RDKit reactions through the AutoMol package."""
    enumerated_graph = nx.DiGraph()

    # == Bulky definitions ==
    def _add_stationary(molecules, role) -> list:
        amchis_list, smiles_list = [], []

        for molecule in molecules:
            (smiles,) = mol_to_smiles(molecule)
            amchi, graph = canonical_enantiomer(stationary_graph(smiles))

            geom = automol.graph.geometry(graph)
            xyz = automol.geom.xyz_string(geom)
            enumerated_graph.add_node(amchi, smiles=smiles, xyz=xyz, role=role)

            smiles_list.append(automol.amchi.smiles(amchi))
            amchis_list.append(amchi)

        return amchis_list, smiles_list

    def _add_transition(reaction: CT.AutomolGraph):
        graph = transition_graph(reaction)

        amchi = automol.graph.amchi(graph)
        geom = automol.graph.geometry(graph)
        xyz = automol.geom.xyz_string(geom)

        def _build_scan():
            formed = automol.graph.ts.forming_bond_keys(graph)
            broken = automol.graph.ts.breaking_bond_keys(graph)

            dmat_angstrom = automol.geom.distance_matrix(geom) * 0.529177

            for broken_bond in broken:
                a, b = broken_bond
                if len(formed) > 0:
                    for formed_bond in formed:
                        shared = broken_bond & formed_bond
                        if len(shared) == 1:
                            c = next(iter(formed_bond - shared))
                            dist = dmat_angstrom[b, c]
                            active_atoms = f"{b} {c}"
                            scan = f"scan B {active_atoms} = {dist:.3f}, 0.7, 100"
                            
                else:
                    active_atoms = f"{a} {b}"
                    scan = f"scan B {active_atoms} = {dmat_angstrom[a, b]:.3f}, 2.0, 100"

            return scan, active_atoms

        scan, active_atoms = _build_scan()

        enumerated_graph.add_node(
            amchi,
            xyz=xyz,
            scan=scan,
            active_atoms=active_atoms,
            role="transition",
            reactants=reactant_amchis,
            products=product_amchis,
        )

        return amchi

    # == Workflow ==
    (reactant_amchis, reactant_smiles) = _add_stationary(reactants, "reactant")

    for products in product_sets:
        (product_amchis, product_smiles) = _add_stationary(products, "product")

        for reaction in reaction_graphs(tuple(reactant_smiles), tuple(product_smiles)):
            transition_amchi = _add_transition(reaction)

            for r_amchi in reactant_amchis:
                enumerated_graph.add_edge(r_amchi, transition_amchi)

            for p_amchi in product_amchis:
                enumerated_graph.add_edge(transition_amchi, p_amchi)

    return enumerated_graph


def reaction_graphs(
    reactants: CT.SMILES | CT.SMILES_Set,
    products: CT.SMILES | CT.SMILES_Set,
) -> tuple[CT.AutomolGraph, ...]:
    """Generates a reaction Automol graph from SMILES strings."""
    reactants = (reactants,) if isinstance(reactants, str) else reactants
    products = (products,) if isinstance(products, str) else products

    return automol.reac.from_smiles(reactants, products)


def stationary_graph(
    smiles: str = None, amchi: str = None, canonical: str = True
) -> CT.AutomolGraph:
    """Generates a stationary AutoMol graph from SMILES or AMChI string."""
    assert smiles is not None or amchi is not None, (
        "A SMILES or AMChI string must be provided to create a graph."
    )

    graph = (
        automol.smiles.graph(smiles)
        if smiles is not None
        else automol.amchi.graph(amchi)
    )
    if canonical:
        return automol.graph.canonical(graph)

    return graph


def transition_graph(
    reaction: CT.AutomolGraph, canonical: str = True
) -> CT.AutomolGraph:
    """Generates a transition AutoMol graph from AutoMol reaction."""
    graph = automol.reac.ts_graph(reaction)

    if canonical:
        return automol.graph.canonical(graph)

    return graph

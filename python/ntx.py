import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

from .ref import CustomTypes as CT


def graph_matcher(molecules_1: CT.RDMols, molecules_2: CT.RDMols):
    """Determines isomorphism of one molecule onto another."""
    graph_1 = networkx_molecular_graph(molecules_1)
    graph_2 = networkx_molecular_graph(molecules_2)

    def _node_match(node_1, node_2) -> bool:
        return node_1["symbol"] == node_2["symbol"]

    return GraphMatcher(graph_1, graph_2, node_match=_node_match)


def networkx_molecular_graph(mols: CT.RDMols) -> CT.NetworkXGraph:
    """Generates a NetworkX graph from RDKit molecule."""
    graph = nx.Graph()

    for idx, mol in enumerate(mols):
        for atom in mol.GetAtoms():
            key = (idx, atom.GetIdx())
            graph.add_node(key, symbol=atom.GetSymbol())

        for bond in mol.GetBonds():
            key_1 = (idx, bond.GetBeginAtomIdx())
            key_2 = (idx, bond.GetEndAtomIdx())
            graph.add_edge(key_1, key_2)

    return graph

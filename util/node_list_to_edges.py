import networkx
from util.pairwise import pairwise


def node_list_to_edges(node_list: list, graph: networkx.Graph=None):
    """
    Takes a list of nodes and bundles them pairwise to edges.
    If a networkx.Graph is passed as an argument, an Exception is raised if an
    edge should be generated which is not contained in the graph.

    Example:
        node_list_to_edges([1, 2, 3, 4, 5])
        # returns
        [(1, 2), (2, 3), (3, 4), (4, 5)]
    """
    edge_list = []
    for node1, node2 in pairwise(node_list):
        edge = (node1, node2)
        if edge not in graph.edges():
            edge = (node2, node1)
            if edge not in graph.edges():
                raise Exception(
                    "Edge " + str(edge) + " does not exist in graph")
        edge_list.append(edge)
    return edge_list

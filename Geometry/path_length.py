import networkx as nx


def is_path_length_leq_k(graph, node1, node2, k):
    """
    Checks whether the shortest path length of two nodes node1 and node2
    in a networkX.Graph graph is less than or equals to k.

    :param graph: networkX.Graph
    :param node1: node ID
    :param node2: node ID
    :param k: int
    :return: True if path length is less than or equals k, False otherwise
    """
    try:
        path_length = nx.shortest_path_length(graph, node1, node2)
        if path_length > k or path_length is -1:
            return False
    except nx.NetworkXNoPath:
        return False
    return True

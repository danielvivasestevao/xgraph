import networkx as nx
from Cluster.tools import apply_clustering


def cluster_maximal_independent_set(graph: nx.Graph):
    """
    :param graph: A networkX.Graph
    """

    graph = graph.copy()

    # maximal independent set; consists of all cluster heads
    mis_set = set(nx.maximal_independent_set(graph))
    nodes_set = set(graph.nodes())
    not_clustered = nodes_set - mis_set

    adj = dict()

    for head in mis_set:
        neighbors = set(graph.neighbors(head))
        adj[head] = neighbors
        not_clustered -= neighbors

    if not_clustered:
        print("Something went hella wrong during MIS clustering")

    return apply_clustering(graph, adj)

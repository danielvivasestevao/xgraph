import networkx as nx
import random
from Cluster.tools import apply_clustering


def cluster_naive(graph, k):
    """
    :param graph: A networkx.Graph
    :param k: range that should be clustered
    """
    # print("Clustering...")

    graph = graph.copy()

    if k is 0:
        return graph

    not_clustered = set(graph.nodes())
    if len(not_clustered) < 3:
        return graph
    # heads = set()
    adj = dict()  # brauche ich fuer spaeter, um alle Kanten neu zu verteilen

    while not_clustered:
        # print("--Not clustered nodes: ", not_clustered)
        head = random.sample(not_clustered, 1)  # get a non-clustered node
        head = head[0]
        # print("--Added cluster head:  ", head)
        # heads.add(head)  # turn node into a cluster head
        not_clustered -= set([head])

        # get all k-hop distant nodes of cluster head
        nbrs_k = nx.single_source_shortest_path_length(graph, head, cutoff=k)
        nbrs_k_set = set()
        [nbrs_k_set.add(node) for node in nbrs_k if (0 < nbrs_k[node] <= k) &
         (node in not_clustered)]
        # print("- Neighbors:        ", nbrs_k_set)
        adj[head] = nbrs_k_set
        not_clustered -= nbrs_k_set

    return apply_clustering(graph, adj)

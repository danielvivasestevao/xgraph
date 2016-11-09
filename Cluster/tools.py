import networkx


def edge_exists(graph, cluster1, cluster2):
    for node1 in cluster1:
        for node2 in cluster2:
            try:
                graph[node1][node2]  # check if edge exists
            except KeyError:
                continue  # edge does not exist
            # edge exists
            return True
    return False


def apply_clustering(graph: networkx.Graph, cluster_dict):
    """
    Takes a graph and a dictionary of [cluster head - cluster member list]
    and changes the graph to represent the clustered graph, i.e. only display
    the cluster heads and their connectivity.
    :param graph: a networkx.Graph
    :param cluster_dict: cluster heads as keys, cluster members as values
    :return:
    """
    # create adjacency matrix for clusters
    cluster_adj = dict()

    if len(cluster_dict) > 1:
        for cluster_head1 in cluster_dict:
            for cluster_head2 in cluster_dict:
                if cluster_head1 == cluster_head2:
                    continue
                cluster1 = cluster_dict[cluster_head1] | {cluster_head1}
                cluster2 = cluster_dict[cluster_head2] | {cluster_head2}
                if edge_exists(graph, cluster1, cluster2):
                    try:
                        cluster_adj[cluster_head1].append(cluster_head2)
                    except KeyError:
                        cluster_adj[cluster_head1] = [cluster_head2]

        # remove all nodes that are not cluster heads
        for node in graph.nodes():
            if node not in cluster_adj:
                graph.remove_node(node)

        # remove all old edges
        graph.remove_edges_from(graph.edges())

        for n1 in cluster_adj:
            for n2 in cluster_adj[n1]:
                graph.add_edge(n1, n2)

    else:  # only one cluster head
        for node in graph.nodes():
            if node not in cluster_dict:
                graph.remove_node(node)

    return graph

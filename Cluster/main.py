from Cluster.maximal_independent_set import cluster_maximal_independent_set
from Cluster.naive import cluster_naive

"""
Adding a new clustering algorithm:
-add its name to the cluster_types list at the end of this file
-add an entry for it in the get_cluster_alg function which returns its
 clustering function for its name
"""


class ClusterTypeNotFoundError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


def get_cluster_alg(cluster_type: str):
    if cluster_type == "NAIVE":
        return cluster_naive
    elif cluster_type == "MIS":
        return cluster_maximal_independent_set
    raise ClusterTypeNotFoundError


cluster_types = [
    "NAIVE",
    "MIS"
]
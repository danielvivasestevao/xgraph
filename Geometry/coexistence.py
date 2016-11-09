import datetime
import os
from pickle import dump, load

import networkx as nx

from Geometry import stats
from Geometry.path_length import is_path_length_leq_k
from Geometry.triangles import get_nodes_in_triangles_dict
from util import Files


class Coexistence:
    """
    Calculates for all nodes of a graph nx.Graph() if it lies inside a triangle
    and if so, whether the node fulfills the weak or strong coexistence
    property for a triangle. If a node lies inside more than one triangle,
    the strongest fulfillment is stored.

    It is possible to look up for each individual triangle-node-tuple whether
    it fulfills the strong or weak redundancy property for any given k.

    Example:
    triangle = (node1, node2, node3)
    coexistence_obj[(triangle, node)] = {'strong': {1: 0, 2: 1},
                                         'weak': {1: 1, 2: 1}}
    """
    def __init__(self, k, graph=None, coe_event_list=None):
        self._range = k  # range of maximum allowed path length to be tested
        self._dict = {}
        self._stats_dict_a = {}  # absolute
        self._stats_dict_p = {}  # percentage

        event_counter = -1  # logging

        if graph is not None:
            if nx.number_connected_components(graph) > 1:
                raise Exception("Graph contains subgraphs")

            print(datetime.datetime.now().date().isoformat(),
                  datetime.datetime.now().time().isoformat(),
                  ": Calc coexistence...",
                  flush=True)  # logging
            print(datetime.datetime.now().date().isoformat(),
                  datetime.datetime.now().time().isoformat(),
                  ": \tCalc nodes in triangles...",
                  flush=True)  # logging

            if coe_event_list:
                self._dict = get_nodes_in_triangles_dict(coe_event_list=coe_event_list)
            else:
                self._dict = get_nodes_in_triangles_dict(graph=graph)

            event_number = len(self._dict)  # logging

            for triangle_node_tuple in self._dict:
                event_counter += 1  # logging
                if event_counter % 10000 == 0:
                    print(datetime.datetime.now().date().isoformat(),
                          datetime.datetime.now().time().isoformat(),
                          ": \tchecked", event_counter, "of",
                          event_number, "node-in-triangle-events",
                          flush=True)  # logging
                triangle = triangle_node_tuple[0]
                node = triangle_node_tuple[1]
                for i in range(1, k+1):
                    # Sobald für ein i strong gilt, gilt für alle j>1
                    # auch strong; dasselbe für weak
                    if (i > 1 and self._dict[(triangle, node)]
                            ['strong'][i-1] == 1):
                        self._dict[(triangle, node)]['strong'][i] = 1
                        self._dict[(triangle, node)]['weak'][i] = 1
                        continue
                    else:
                        if check_strong_coexistence(graph, node,
                                                    triangle, i):
                            self._dict[(triangle, node)]['strong'][i]\
                                = 1
                            self._dict[(triangle, node)]['weak'][i]\
                                = 1
                            continue
                        else:
                            self._dict[(triangle, node)]['strong'][i]\
                                = 0

                    if (i > 1 and self._dict[(triangle, node)]
                            ['weak'][i-1] == 1):
                        self._dict[(triangle, node)]['weak'][i] = 1
                    else:
                        if (check_weak_coexistence(graph, node,
                                                   triangle, i)):
                            self._dict[(triangle, node)]['weak'][i] = 1
                        else:
                            self._dict[(triangle, node)]['weak'][i] = 0

            print(datetime.datetime.now().date().isoformat(),
                  datetime.datetime.now().time().isoformat(),
                  ": Calc coexistence DONE",
                  flush=True)  # logging

    def get_stats(self, percent=True):
        stats_dict = self._stats_dict_a
        if percent:
            stats_dict = self._stats_dict_p

        if stats_dict:
            return stats_dict
        else:
            if percent:
                self._stats_dict_p = stats.get_stats(self, percent)
                return self._stats_dict_p
            else:
                self._stats_dict_a = stats.get_stats(self, percent)
                return self._stats_dict_a

    def get_dict(self):
        return self._dict

    def save(self, file_path: "relative path to safe file"):
        """
        Save Coexistence object to a binary file. The file can be loaded via
        Coexistence.load().

        :param file_path: Relative path to the Coexistence file.
        """
        # cut off file name from path and make each directory in it (if needed)
        Files.ensure_dir(os.path.dirname(file_path))
        with open(file_path, "wb") as f:
            try:
                dump(self, f, 4)  # pickle.dump
            except MemoryError:
                dump(self, f)

    @staticmethod
    def load(file_path: "relative path to file"):
        """
        Load a Coexistence object from a binary file. The file must have been
        saved via Coexistence.save().

        :param file_path: Relative path to the Coexistence file.
        :return: The Coexistence object from the input file.
        """
        with open(file_path, "rb") as f:
            return load(f)


def check_strong_coexistence(graph, node, triangle, k):
    return check_coexistence(graph, node, triangle, k, weak=False)


def check_weak_coexistence(graph, node, triangle, k):
    return check_coexistence(graph, node, triangle, k, weak=True)


def check_coexistence(graph, node, triangle, k, weak=True):
    """
    Checks for a node and a triangle (three line segments defined by the
    positions of two nodes each) whether they fulfill the weak/strong
    coexistence property, i.e. the node has a path with k or less hops to
    [at least one node of the triangle (weak)]/[all nodes of the triangle
    (strong)].

    :param graph: the networkx.Graph containing the node and the nodes of the
    triangle
    :param node: a node
    :param triangle: three nodes with edges forming a triangle
    :param k: maximum allowed path length
    :param weak: True if you want to check weak coexistence, False if you
    want to check strong coexistence
    :return: True if property is fulfilled, False otherwise
    """
    b1 = is_path_length_leq_k(graph, node, triangle[0], k)
    b2 = is_path_length_leq_k(graph, node, triangle[1], k)
    b3 = is_path_length_leq_k(graph, node, triangle[2], k)

    if weak:
        return b1 or b2 or b3
    else:
        return b1 and b2 and b3

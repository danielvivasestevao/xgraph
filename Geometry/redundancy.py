import datetime
import os
from pickle import dump, load

import networkx as nx

from Geometry import stats
from Geometry.intersections import get_intersections_dict
from Geometry.path_length import is_path_length_leq_k
from util import Files


class Redundancy:
    """
    Calculates for all intersections of a graph nx.Graph() whether the two
    edges or rather the four nodes of the two intersecting edges fulfill the
    weak or strong redundancy property.

    It is possible to look up for each individual intersection whether it
    fulfills the strong or weak redundancy property for any given k.

    Example:
    intersection = (('A', 'B'), ('C', 'D'))
    redundancy_obj[intersection] = {{strong: {1: 0}}, {weak: {1: 1}}}
    """
    def __init__(self, k, graph=None, intersection_list=None):
        self._edge_cut = 0  # number of intersecting edges in the graph
        # self.intersections = []  # list of tuples of intersecting edges
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
                  ": Calc redundancy...",
                  flush=True)  # logging
            print(datetime.datetime.now().date().isoformat(),
                  datetime.datetime.now().time().isoformat(),
                  ": \tCalc intersections...",
                  flush=True)  # logging

            if intersection_list:
                self._dict = get_intersections_dict(intersections=intersection_list)
            else:
                self._dict = get_intersections_dict(graph=graph)

            event_number = len(self._dict)  # logging

            for edge_tuple in self._dict:
                event_counter += 1  # logging
                if event_counter % 100000 == 0:
                    print(datetime.datetime.now().date().isoformat(),
                          datetime.datetime.now().time().isoformat(),
                          ": \tcalculated", event_counter, "of",
                          event_number, "intersections",
                          flush=True)  # logging
                for i in range(1, k+1):
                    # Sobald für ein i strong gilt, gilt für alle j>1
                    # auch strong; dasselbe für weak
                    if i > 1 and self._dict[edge_tuple]['strong'][i-1] == 1:
                        self._dict[edge_tuple]['strong'][i] = 1
                        self._dict[edge_tuple]['weak'][i] = 1
                        continue
                    else:
                        if check_strong_redundancy(graph, edge_tuple[0],
                                                   edge_tuple[1], i):
                            self._dict[edge_tuple]['strong'][i] = 1
                            self._dict[edge_tuple]['weak'][i] = 1
                            continue
                        else:
                            self._dict[edge_tuple]['strong'][i] = 0

                    if i > 1 and self._dict[edge_tuple]['weak'][i-1] == 1:
                        self._dict[edge_tuple]['weak'][i] = 1
                    else:
                        if check_weak_redundancy(graph, edge_tuple[0],
                                                 edge_tuple[1], i):
                            self._dict[edge_tuple]['weak'][i] = 1
                        else:
                            self._dict[edge_tuple]['weak'][i] = 0

            self._edge_cut = len(self._dict)
            print(datetime.datetime.now().date().isoformat(),
                  datetime.datetime.now().time().isoformat(),
                  ": Calc redundancy DONE",
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
        Save Redundancy object to a binary file. The file can be loaded via
        Redundancy.load().

        :param file_path: Relative path to the Redundancy file.
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
        Load a Redundancy object from a binary file. The file must have been
        saved via Redundancy.save().

        :param file_path: Relative path to the Redundancy file.
        :return: The Redundancy object from the input file.
        """
        with open(file_path, "rb") as f:
            return load(f)


def check_strong_redundancy(graph, line1, line2, k):
    return check_redundancy(graph, line1, line2, k, weak=False)


def check_weak_redundancy(graph, line1, line2, k):
    return check_redundancy(graph, line1, line2, k, weak=True)


def check_redundancy(graph, line1, line2, k, weak=True):
    """
    Checks for two line segments defined by the positions of four nodes whether
    they fulfill the weak/strong redundancy property, i.e. both line segments
    are edges in the graph and they intersect and at least one node has a path
    with k or less hops to [at least one node (weak)]/[both nodes (strong)] of
    the respective other line segment.

    IMPORTANT NOTE: This function does *not* check whether the two line
    segments exist in the graph and whether they intersect or not. This has to
    be done beforehand.

    :param graph: the networkx.Graph containing the nodes of line1 and line2
    :param line1: a tuple of nodes which defines a line segment
    :param line2: a tuple of nodes which defines a line segment
    :param k: maximum allowed path length
    :param weak: True if you want to check weak redundancy, False if you
    want to check strong redundancy
    :return: True if property is fulfilled, False otherwise
    """

    b00 = is_path_length_leq_k(graph, line1[0], line2[0], k)
    b01 = is_path_length_leq_k(graph, line1[0], line2[1], k)
    b10 = is_path_length_leq_k(graph, line1[1], line2[0], k)
    b11 = is_path_length_leq_k(graph, line1[1], line2[1], k)

    if weak:
        return b00 or b01 or b10 or b11
    else:
        return b00 and b01 or b10 and b11 or b00 and b10 or b01 and b11

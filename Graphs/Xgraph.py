# GUI
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,\
    NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import tkinter as Tk

import networkx as nx

from util.dict_funcs import first_one
from util.node_list_to_edges import node_list_to_edges

# Xgraph, XgraphList
import Cluster.main
from Geometry import *
from Geometry.Lognormal.lognorm_edge import lognorm_edge
from Geometry import intersections
from Geometry.triangles import get_triangles
from matplotlib.pyplot import show, clf, savefig
from multiprocessing import Pool
import numpy
import os
from pickle import dump, load
import platform
from PPP.poisson_pp import spatial_poisson_gen
import shutil
import subprocess
from util import dict_funcs, ProcessTimer, Files
from uuid import uuid4, UUID
import warnings
import Xgraph_testing.rust.rust2python as r2p


# TODO include example code in method docs


class Xgraph:
    """
    A wrapper class for networkx.Graph() which expands it by an ID and a myriad
    of useful functions for the reactive spanner project.

    See the PyDoc of the __init__ function for more information.
    """

    def __init__(self, gps: bool=False, lognormal: bool=False,
                 path_loss_exponent: float=None,
                 variance: float=None, transmit_power: float=None,
                 path_loss_one_meter: float=None, poisson_point: bool=False,
                 lam: int=None, dim_x: int=None, dim_y: int=None):
        """
        Constructor. Called when a new Xgraph object is generated.

        :param gps: If True, the node positions are treated as gps coordinates;
         consequently, distances between nodes are calculated based on the
         haversine formula. If false, node positions are treated as given in
         meter on a cartesian plane.
        :param lognormal: If True, edges from each newly added node to all
         other nodes are added automatically if the RSSI value of the
         connection between two nodes is higher than -95 dBm when calculated
         with the lognormal shadowing interference model.
         If True, the following parameters must not be None:
         path_loss_exponent, variance, transmit_power, path_loss_one_meter.
        :param path_loss_exponent: Path loss exponent of the lognormal
         shadowing interference model.
        :param variance: Variance of the lognormal shadowing interference
         model. Please not that [standard deviation = sqrt(variance)].
        :param transmit_power: Transmit power of the simulated nodes in dBm.
         Needed in the lognormal shadowing interference model.
        :param path_loss_one_meter: The reference path loss over one meter in
         dB. Used in the lognormal shadowing interference model.
        :param poisson_point: If True, nodes will be generated according to a
         Poisson Point process.
         If True, the following parameters must not be None:
         lam, dim_x, dim_y.
        :param lam: Lambda value of the Poisson Point process, i.e. the
         expected number of nodes over an area of size [dim_x * dim_y].
        :param dim_x: x-axis of the simulation area for the Poisson Point
         process. [in meter if gps is False]
        :param dim_y: y-axis of the simulation area for the Poisson Point
         process. [in meter if gps is False]
        """
        self._graph = nx.Graph()  # networkx Graph with undirected edges
        self._id = uuid4()
        self._red = None  # redundancy object
        self._coe = None  # coexistence object
        self._cluster = dict()
        self._error_flag = False

        # used to save networkx.get_node_attributes(graph, 'pos') output
        # and prevent unnecessary calls
        self._node_attributes = None

        # Properties
        self._properties = dict()
        self._properties["gps"] = gps

        self._properties["lognormal"] = lognormal
        self._properties["ple"] = path_loss_exponent
        self._properties["variance"] = variance
        self._properties["transmit_power"] = transmit_power
        self._properties["path_loss_one_meter"] = path_loss_one_meter
        if lognormal:
            if (path_loss_exponent is None or variance is None or
                    transmit_power is None or path_loss_one_meter is None):
                raise Exception("Graph is set as lognormal, but no "
                                "path loss exponent, variance, transmit power "
                                "or path loss over one meter given")

        self._properties["ppp"] = poisson_point
        self._properties["lam"] = lam
        if dim_x is None or dim_x > 0:
            self._properties["dim_x"] = dim_x
        else:
            raise Exception("Dimension has to be greater than 0")
        if dim_y is None:  # If only dim_x is given, area is square
            self._properties["dim_y"] = dim_x
        else:
            self._properties["dim_y"] = dim_y
        if poisson_point:
            if lam is None or dim_x is None:
                raise Exception("Point positions are created according to "
                                "a Poisson Point Process, but no dimensions "
                                "or lambda are given")
            self._generate_pp_graph()
            if lognormal:
                counter = 0
                try:
                    while nx.edge_connectivity(self._graph) == 0:
                        self._graph.clear()
                        self._generate_pp_graph()
                        counter += 1
                        if counter > 100:
                            raise Exception("Chance of creating unconnected "
                                            "graph is too high.")
                except MemoryError:
                    print("Number of edges too high:",
                          len(self._graph.edges()))
                    self._error_flag = True

    # TODO __eq__ - Check redundancy/coexistence?
    def __eq__(self, other):
        """
        Two Xgraph objects are equal if they have the same properties and their
        underlying networkx.Graph objects have the same nodes and edges,
        including node attributes and edge weights.
        Xgraph objects with different IDs can be equal.

        :param other: object class
        :return: bool
        """

        if isinstance(other, Xgraph):
            if not self._eq_properties(other):
                return False
            if self._graph.node != other._graph.node:
                return False
            if self._graph.edge != other._graph.edge:
                return False
        return True

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def _add_node_lognormal(self, node_id: int, pos: (float, float)):
        """
        Add a node to the graph, and calculate whether edges from this node to
        any other node in the graph exist according to the lognormal shadowing
        interference model.

        :param node_id: ID of the added node as an int.
        :param pos: Position of the added node as a tuple of floats.
        """
        pos1 = pos
        self._graph.add_node(node_id, pos=pos1)
        for node in self._graph.nodes_iter():
            if node == node_id:  # do not add edge to node itself
                continue
            pos2 = self._graph.node[node]["pos"]
            rssi = lognorm_edge(self._properties["ple"],
                                self._properties["variance"],
                                self._properties["transmit_power"],
                                self._properties["path_loss_one_meter"],
                                pos1, pos2, self._properties["gps"])
            if rssi:
                self._graph.add_edge(node, node_id, rssi=rssi)

    def _eq_properties(self, other):
        """
        :param other: Xgraph object
        :return: True if the properties of this Xgraph and other are identical;
         False otherwise.
        """
        for prop in self._properties:
                if self.get_property(prop) != other.get_property(prop):
                    return False
        return True

    def add_node(self, node_id: int, pos: (float, float)):
        """
        Adds a node to the graph.
        If the graph is a lognormal graph, edges from the added node to all
        nodes of the graph are calculated.

        :param node_id: ID of the added node as an int.
        :param pos: Position of the added node as a tuple of floats.
        """
        try:
            self._node_attributes = None
        except AttributeError:  # for older Xgraph versions
            pass
        if self._properties["lognormal"]:
            self._add_node_lognormal(node_id, pos)
        else:
            self._graph.add_node(node_id,
                                 pos=(float(pos[0]), float(pos[1])))
        # print("degree(", node_id, "):", self._graph.degree(node_id))

    def avrg_node_degree(self):
        """
        :return: Average node degree of the graph, i.e. the average number of
        edges adjacent to a node.
        """
        if not self._graph.nodes():
            return 0
        degree_dict = self._graph.degree()
        degree_list = list(degree_dict.values())
        a = numpy.array(degree_list)
        return a.mean()

    def calc_coexistence(self, k: int=None):
        """
        Calculate coexistence of the Xgraph up to a path length of k hops.
        If no k is given, calculate coexistence up to the Xgraph's diameter.
        This method uses a Rust binary program to perform all point-in-triangle
        tests. For a pure Python version, use calc_coexistence_py.

        :param k: range of maximum allowed path length to be tested
        """
        if platform.system() == "Linux":
            binary_name = "redcoe"
        else:
            binary_name = "redcoe.exe"

        point2pos = self.get_node_positions()
        point2pos_file = "point2pos_" + str(self._id) + ".txt"
        point_pos_tup_str = []
        for id in point2pos:
            point_pos_tup_str.append(str(id))
            point_pos_tup_str.append(str(point2pos[id][0]))
            point_pos_tup_str.append(str(point2pos[id][1]))
        point_pos_tup_str = ' '.join(point_pos_tup_str)
        with open(point2pos_file, "w") as f:
            f.write(point_pos_tup_str)

        triangle_list_file = "triangle_list_" + str(self._id) + ".txt"
        triangle_list_str = []
        for triangle in list(get_triangles(self._graph, True, point2pos)):
            triangle_list_str.append(str(triangle.node1))
            triangle_list_str.append(str(triangle.node2))
            triangle_list_str.append(str(triangle.node3))
        triangle_list_str = ' '.join(triangle_list_str)
        with open(triangle_list_file, "w") as f:
            f.write(triangle_list_str)

        if not k:
            k = self.diameter()

        p = ProcessTimer.ProcessTimer("Calculate coexistence events with Rust")
        p.start()
        coe_events = \
            r2p.coe_event_list_str2list(
                subprocess
                    .check_output(
                        [os.path.join(os.path.dirname(__file__), "..", "Rust_files", binary_name),
                         "false", point2pos_file, triangle_list_file])
                    .decode(encoding='utf-8'),
                self.get_node_positions()
            )
        p.end()

        self._coe = coexistence.Coexistence(
            k, self._graph, coe_event_list=coe_events)

        os.remove(point2pos_file)
        os.remove(triangle_list_file)

    def calc_coexistence_py(self, k: int=None):
        """
        Calculate coexistence of the Xgraph up to a path length of k hops.
        If no k is given, calculate coexistence up to the Xgraph's diameter.

        :param k: range of maximum allowed path length to be tested
        """
        if k:
            self._coe = coexistence.Coexistence(k, self._graph)
        else:
            self._coe = coexistence.Coexistence(self.diameter(),
                                                self._graph)

    def calc_redundancy(self, k: int=None):
        """
        Calculate redundancy of the Xgraph up to a path length of k hops.
        If no k is given, calculate redundancy up to the Xgraph's diameter.
        This method uses a Rust binary program to calculate all intersections
        in the graph. For a pure Python version, use calc_redundancy_py.

        :param k: range of maximum allowed path length to be tested
        """
        if platform.system() == "Windows":
            binary_name = "redcoe.exe"
        else:
            binary_name = "redcoe"

        point2pos = self.get_node_positions()
        point2pos_file = "point2pos_" + str(self._id) + ".txt"
        point_pos_tup_str = []
        for id in point2pos:
            point_pos_tup_str.append(str(id))
            point_pos_tup_str.append(str(point2pos[id][0]))
            point_pos_tup_str.append(str(point2pos[id][1]))
        point_pos_tup_str = ' '.join(point_pos_tup_str)
        with open(point2pos_file, "w") as f:
            f.write(point_pos_tup_str)

        edge_list_file = "edge_list_" + str(self._id) + ".txt"
        edge_list_str = []
        for edge_tup in self._graph.edges():
            edge_list_str.append(str(edge_tup[0]))
            edge_list_str.append(str(edge_tup[1]))
        edge_list_str = ' '.join(edge_list_str)
        with open(edge_list_file, "w") as f:
            f.write(edge_list_str)

        if not k:
            k = self.diameter()

        p = ProcessTimer.ProcessTimer("Calculate redundancy events with Rust")
        p.start()
        intersections = \
            r2p.edge_list_str2list(
                subprocess.check_output([
                    os.path.join(
                        os.path.dirname(__file__), "..", "Rust_files", binary_name),
                    "true", point2pos_file, edge_list_file])
                .decode(encoding='utf-8')
            )
        p.end()

        self._red = redundancy.Redundancy(
            k, self._graph, intersection_list=intersections)

        os.remove(point2pos_file)
        os.remove(edge_list_file)

    def calc_redundancy_py(self, k: int=None):
        """
        Calculate redundancy of the Xgraph up to a path length of k hops.
        If no k is given, calculate redundancy up to the Xgraph's diameter.

        :param k: range of maximum allowed path length to be tested
        """
        if k:
            self._red = redundancy.Redundancy(k, self._graph)
        else:
            self._red = redundancy.Redundancy(self.diameter(),
                                              self._graph)

    def cluster(self, cluster_type: str, k: int=None):
        """
        Clusters the graph and saves the clustered graph in the Xgraph's
        cluster dictionary. The clustered graph is an Xgraph object with
        default parameter values.

        Example:
            xg.cluster("naive", 2)  # xg is an Xgraph object
            xg.cluster("mis")  # 1-hop maximal independent set clustering

        :param cluster_type: The clustering algorithm.
        :param k: The range that should be clustered (if the clustering
         algorithm supports it).
        """
        cluster_type = cluster_type.strip(" ").upper()
        cluster_alg = Cluster.main.get_cluster_alg(cluster_type)

        try:
            clustered_graph = cluster_alg(self._graph, k)
            clustered_xg = Xgraph()
            clustered_xg._graph = clustered_graph
            try:
                self._cluster[cluster_type][k] = clustered_xg
            except KeyError:
                self._cluster[cluster_type] = dict()
                self._cluster[cluster_type][k] = clustered_xg
        except TypeError:
            clustered_graph = cluster_alg(self._graph)
            clustered_xg = Xgraph()
            clustered_xg._graph = clustered_graph
            self._cluster[cluster_type] = clustered_xg

    def diameter(self):
        """
        :return: The Xgraph's diameter, i.e. the length of the 'longest
        shortest path' in the graph.
        """
        if not self._graph.nodes():
            raise Exception("Diameter: Xgraph " + str(self.get_id()) +
                            " has no nodes.")
        return nx.diameter(self._graph)

    def _generate_pp_graph(self):
        """
        Adds nodes with IDs and positions according to a Poisson Point
        Process to the graph. The parameters for the graph generation are the
        ones used in the __init__ method.
        """
        node_pos_list = spatial_poisson_gen(self._properties["lam"],
                                            self._properties["dim_x"],
                                            self._properties["dim_y"])
        i = 0
        for node_pos in node_pos_list:
            self.add_node(i, (node_pos[0], node_pos[1]))
            i += 1

    def get_cluster(self, cluster_type: str, k: int=None):
        """
        If it has been computed already, return the clustered variant of this
        graph.

        Example:
            xg.get_cluster("naive", 2)  # returns an Xgraph object

        :param cluster_type: used clustering algorithm
        :param k: clustering range
        :return: Xgraph object which is a clustered version of this
         Xgraph object's underlying networkx.Graph
        """
        cluster_type = cluster_type.strip(" ").upper()
        if k:
            try:
                return self._cluster[cluster_type][k]
            except TypeError:
                return self._cluster[cluster_type]
        else:
            return self._cluster[cluster_type]

    def get_coexistence(self):
        """
        :return: The coexistence object of the Xgraph. Read the coexistence
        dictionary by calling the get_dict() method on the coexistence object.
        """
        if self._coe:
            return self._coe
        else:
            warnings.warn("No Coexistence object generated yet", UserWarning)

    def get_coexistence_stats(self, percent=True):
        """
        :param percent: False for total values; otherwise percentage values.
        :return: A dictionary of the coexistence data.
        """
        if self._coe:
            return self._coe.get_stats(percent=percent)
        else:
            warnings.warn("No Coexistence object generated yet", UserWarning)

    def get_error_flag(self):
        return self._error_flag

    def get_graph(self):
        """
        :return: The underlying networkx.Graph object of the Xgraph.
        """
        warnings.warn("Manipulating the graph will NOT result in an update of "
                      "the overlying Xgraph's properties.", UserWarning)
        return self._graph

    def get_id(self):
        """
        :return: The Xgraph's ID (generated via UUID4).
        """
        return self._id

    def get_node_positions(self):
        """
        :return: A dictionary with every node ID as a key and the position of
          the node as its value. The node position is given as a float tuple.
        """
        try:
            if not self._node_attributes:
                self._node_attributes = nx.get_node_attributes(self._graph, "pos")
            return self._node_attributes
        except AttributeError:  # for older versions of Xgraph
            return nx.get_node_attributes(self._graph, "pos")

    def get_properties(self):
        warnings.warn("Manipulating properties will NOT result in an update"
                      "of the underlying graph.\n"
                      "Do only edit these values if you know what "
                      "you are doing.")
        return self._properties

    def get_property(self, property_name: str):
        """
        Get a property of the Xgraph. Properties are: gps, lognormal, ple,
        variance, transmit_power, path_loss_one_meter, ppp, lam, dim_x, dim_y.
        For an attribute, use the corresponding get-method.

        :param property_name: name of the property as a string
        """
        try:
            return self._properties[property_name]
        except KeyError:
            raise KeyError("Property not found. Properties are: gps, "
                           "lognormal, ple, variance, transmit_power, "
                           "path_loss_one_meter, ppp (Poisson Point Process), "
                           "lam, dim_x, dim_y. "
                           "For an attribute, use the corresponding "
                           "get-method.")

    def get_redundancy(self):
        """
        :return: The redundancy object of the Xgraph. Read the redundancy
        dictionary by calling the get_dict() method on the redundancy object.
        """
        if self._red:
            return self._red
        else:
            warnings.warn("No Redundancy object generated yet", UserWarning)

    def get_redundancy_stats(self, percent=True):
        """
        :param percent: False for total values; otherwise percentage values.
        :return: A dictionary of the redundancy data.
        """
        if self._red:
            return self._red.get_stats(percent=percent)
        else:
            warnings.warn("No Redundancy object generated yet", UserWarning)

    def node_density(self):
        """
        :return: Node density, i.e. number of nodes per square unit of area.
        """
        if (self._properties["dim_x"] is None or
                self._properties["dim_y"] is None or
                self._properties["dim_x"] == 0 or
                self._properties["dim_y"] == 0):
            raise Exception("Node density: Xgraph " + str(self.get_id()) +
                            "has dim_x not set or area is 0.")
        return self.number_of_nodes() /\
            (self._properties["dim_x"] * self._properties["dim_y"])

    @staticmethod
    def load(file_path: "relative path to file"):
        """
        Load an Xgraph object from a binary file. The file must have been saved
        via Xgraph.save().

        :param file_path: Relative path to the Xgraph file.
        :return: The Xgraph object from the input file.
        """
        with open(file_path, "rb") as f:
            return load(f)

    def number_of_coe_events(self):
        """
        :return: The number of coexistence events, i.e. the number of times a
        node of the graph lies within a triangle of three other connected nodes
        of the graph.
        """
        if self._coe:
            return len(self.get_coexistence()._dict)
        else:
            warnings.warn("No Coexistence object generated yet", UserWarning)

    def number_of_intersections(self):
        """
        :return: The number of edge intersections of the graph.
        """
        return len(intersections.simple_bentley_ottmann(self._graph))

    def number_of_edges(self):
        """
        :return: The number of edges of the graph.
        """
        return nx.number_of_edges(self._graph)

    def number_of_nodes(self):
        """
        :return: The number of nodes of the graph.
        """
        return len(self._graph)

    def number_of_red_events(self):
        """
        :return: The number of redundancy events, i.e. the number of times two
        edges of the graph intersect.
        """
        if self._red:
            return len(self.get_redundancy()._dict)
        else:
            warnings.warn("No Redundancy object generated yet", UserWarning)

    def open_in_gui(self):
        """
        Opens the graph in a GUI to inspect single coexistence and redundancy
        values manually. Coexistence and redundancy of the graph have to be
        calculated beforehand.
        """
        XgraphGUI(self)

    def plot_graph(self, with_labels: bool=True,
                   file_path: "relative path to safe file"=None):
        """
        Plot the underlying networkX.Graph of the Xgraph.
        To save the plot as an image file, use Xgraph.save_graph_plot().

        :param with_labels: If True, node IDs are displayed in the plot.
        :param file_path: Path to the image file if it is to be saved. The
         method save_graph_plot should be used for this instead.
        """
        nx.draw(self._graph,
                alpha=1,
                node_color="black",
                node_size=10,
                pos=self.get_node_positions(),
                with_labels=with_labels)
        if file_path:
            savefig(file_path + ".png", format="PNG")
            clf()
        else:
            show()
        # clf()  # closes figure immediately when called in console mode

    def plot_graph_frank(self, with_labels: bool = True,
                         file_path: "relative path to safe file" = None):
        """
        Plot the underlying networkX.Graph of the Xgraph.
        To save the plot as an image file, use Xgraph.save_graph_plot().

        :param with_labels: bool
        :param file_path: if a path is given, a png file of the plotted
         graph is saved to that file instead of plotting it in a window
        """
        nx.draw(self.get_graph(),
                alpha=1,
                node_color="red",
                node_size=5,
                pos=self.get_node_positions(),
                label_size=5,
                with_labels=with_labels)
        if file_path:
            savefig(file_path + ".png", format="PNG")
            clf()
        else:
            show()

        # clf()  # closes figure immediately when called in console mode

    def pretty_print(self):
        """
        Prints the Xgraph's ID, number of nodes, edges, redundancy and
        coexistence events on the console.
        """
        pretty_str_list = []
        pretty_str_list.extend(list(map
                                    (str,
                                     ["Xgraph", self._id, ":",
                                      "\t", self.number_of_nodes(), "Nodes |",
                                      "\t", self.number_of_edges(), "Edges"])))
        if self._red:
            pretty_str_list.extend(["|", "\t",
                                    str(self.number_of_red_events()),
                                    "Redundancy Events "])
        if self._coe:
            pretty_str_list.extend(["|", "\t",
                                    str(self.number_of_coe_events()),
                                    "Coexistence Events"])
        print(" ".join(pretty_str_list))

    def save(self, file_path: "relative path to safe file"):
        """
        Save an Xgraph object to a binary file. The file can be loaded via
        Xgraph.load().

        :param file_path: Relative path to the Xgraph file.
        """
        # cut off file name from path and make each directory in it (if needed)
        Files.ensure_dir(os.path.dirname(file_path))
        with open(file_path, "wb") as f:
            dump(self, f, 4)  # pickle.dump

    def save_graph_plot(self, file_path: "relative path to safe file",
                        with_labels=True):
        """
        Save a plot of the underlying networkX.Graph of the Xgraph to an image
        file.

        :param file_path: Path to the image file.
        :param with_labels: If True, node IDs are displayed in the plot.
        """
        self.plot_graph(with_labels=with_labels, file_path=file_path)


# TODO refactor using __getitem__ so that _small_mem checking is done in __getitem__ -> causes problems with __iter__, which does not automatically adapt the __getitem__ func
# TODO refactor the xg and xgl saving clusterf!ck with self.__xgraph_dir
class XgraphList(list):
    """
    List of Xgraph objects.
    """

    def __init__(self, graph_list: "list of Xgraph objects" = None,
                 small_mem: bool = False):
        super().__init__()
        self._id = uuid4()
        self._values = dict()

        self._small_mem = small_mem  # flag for systems with limited RAM
        # attributes for _small_mem
        self.__xgraph_dir = None
        self.__xgraph_props = None
        if self._small_mem:
            self._set_small_mem()

        if graph_list:
            for xg in graph_list:
                if (self._small_mem and not isinstance(xg, Xgraph) and not
                   isinstance(xg, str) and not isinstance(xg, UUID)):
                    raise Exception("Parameter is not an Xgraph object or an "
                                    "Xgraph object's ID (UUID or UUID string)."
                                    )
                if not isinstance(xg, Xgraph):
                    raise Exception("Parameter contains elements which are "
                                    "not Xgraph objects.")
                else:
                    self.append(xg)

    def __add__(self, other):
        self.extend(other)

    def __iadd__(self, other):
        self.extend(other)

    def __avrg_numpy_func(self, func_name: str):
        if func_name in self._values:
            return self._values[func_name]

        avrg = self.__generic_numpy_func(func_name).mean()
        self._values[func_name] = avrg
        return avrg

    def __calc_geometry(self, k: int=None, type_: bool=True):
        """
        Calculate coexistence or redundancy of all Xgraph objects in this list
        up to a path length of k hops.
        If no k is given, calculate coexistence or redundancy up to the
        diameter of the Xgraph object with the biggest diameter out of all
        Xgraph objects in this list.

        This method is to be accessed using either calc_coexistence() or
        calc_redundancy().

        :param k: range of maximum allowed path length to be tested
        :param type_: if True, calculate redundancy; else coexistence
        """
        msg = "calc_coexistence XgraphList "
        if type_:
            msg = "calc_redundancy XgraphList "

        if k is None:
            k = self.max_diameter()

        t = ProcessTimer.ProcessTimer(msg + str(self._id))
        t.start()

        calc_geometry_xgraphlist_multiprocessing(self, k, type_)

        t.end()

    def _delete_from_disk(self):
        if not self._small_mem:
            return
        warnings.warn("The whole directory will be removed, even non-Xgraph "
                      "files!")
        os.remove(str(os.path.abspath(os.curdir)) + "/" +
                  str(self.get_id()) + ".xgl")
        shutil.rmtree(self.__xgraph_dir)

    def __generic_numpy_func(self, func_name: str):
        if not self:
            raise Exception("The XgraphList is empty.")

        func_name_dict = {
            "avrg_diameter": "diameter",
            "avrg_node_number": "number_of_nodes",
            "avrg_node_degree": "avrg_node_degree",
            "avrg_node_density": "node_density",
            "avrg_number_of_coe_events": "number_of_coe_events",
            "avrg_number_of_red_events": "number_of_red_events",
            "max_diameter": "diameter"
        }
        if self._small_mem:
            a = self.__get_xg_function_list(func_name_dict[func_name])
        else:
            a = numpy.array([getattr(xg, func_name_dict[func_name])()
                             for xg in self])
        return a

    # TODO description
    def __get_geometry_stats(self, percent=True, type_: bool=True):
        """

        :param percent:
        :param type_: redundancy if True, coexistence otherwise
        :return:
        """
        geometry_stats_list = list()
        for xg in self:
            if self._small_mem:
                xg = Xgraph.load(self.__xgraph_dir + "/" + str(xg) + ".xg")
            if type_:
                geometry_stats = xg.get_redundancy_stats(percent=percent)
            else:
                geometry_stats = xg.get_coexistence_stats(percent=percent)
            geometry_stats_list.append(geometry_stats)
        return dict_funcs.stats_mean(geometry_stats_list)

    def __get_xg_function_list(self, func_name: str):
        """
        Used only in case of _small_mem = True.
        Creates a list of values returned from using a certain function on
        all Xgraph objects loaded from the Xgraph IDs in this list.

        :param func_name: name of an Xgraph function with a return value
        :return: numpy.array
        """
        return numpy.array([getattr(
            Xgraph.load(self.__xgraph_dir + "/" + str(xg_id) + ".xg"),
            func_name)()
                            for xg_id in self])

    def __max_numpy_func(self, func_name: str):
        if func_name in self._values:
            return self._values[func_name]

        max_ = self.__generic_numpy_func(func_name).max()
        self._values[func_name] = max_
        return max_

    def __reset_values(self):
        if self._values:
            self._values = dict()

    def _set_small_mem(self):
        """
        Sets the a flag for systems with limited memory. Instead of loading
        all Xgraph objects into memory, only their IDs are saved in the list,
        while the objects themselves are saved as files to the disk. Whenever
        an Xgraph object is needed, its respective file is loaded into memory
        until the task is finished.
        """
        self.__xgraph_dir = str(os.path.abspath(os.curdir)) +\
                            "/Xgraph_files_" + str(self.get_id())
        Files.ensure_dir(self.__xgraph_dir)

        id_list = list()
        for xg in self:
            id_list.append(xg.get_id())
            file_path = self.__xgraph_dir + str(xg.get_id()) + ".xg"
            # TODO Should the file be saved anyways?
            if not Files.file_exists(file_path):
                xg.save(self.__xgraph_dir + str(xg.get_id()) + ".xg")
            self.remove_graph(xg)
        self._small_mem = True
        self.extend(id_list)
        self.save(str(self.get_id()) + ".xgl")

    def append(self, xgraph):
        """
        Appends an Xgraph object to the XgraphList.

        :param xgraph: Xgraph object to be appended to the list.
         If _small_mem is set, xgraph can also be an Xgraph object's ID (as
         UUID or UUID string).
        """
        if self._small_mem:
            if isinstance(xgraph, Xgraph):
                super().append(xgraph.get_id())
                xgraph.save(self.__xgraph_dir + "/" +
                            str(xgraph.get_id()) + ".xg")
            elif isinstance(xgraph, str) or isinstance(xgraph, UUID):
                super().append(str(xgraph))
            else:
                raise Exception("Parameter is not an Xgraph object or an "
                                "Xgraph object's ID (UUID or UUID string).")
            # TODO Problem: saves XgraphList within xgraph_dir, creates second xgraph_dir folder
            # self.save(str(self.get_id()) + ".xgl")
        else:
            if not isinstance(xgraph, Xgraph):
                raise Exception("Parameter is not an Xgraph object.")
            super().append(xgraph)

        self.__reset_values()

    def extend(self, iterable):
        """
        Appends a list of Xgraph objects to the list.

        :param iterable: List of Xgraph objects.
        """
        [self.append(obj) for obj in iterable]

    def are_properties_homogeneous(self):
        """
        :return: True if all Xgraph objects in this list have the same
         property values; False otherwise.
        """
        if len(self) == 0:
            return True

        if self._small_mem:
            xg1 = Xgraph.load(self.__xgraph_dir + "/" + str(self[0]) + ".xg")
            for xg_id in self:
                xg = Xgraph.load(self.__xgraph_dir + "/" + str(xg_id) + ".xg")
                if not xg1._eq_properties(xg):
                    return False
            return True

        xg1 = self[0]
        for xg in self:
            if not xg1._eq_properties(xg):
                return False
        return True

    def avrg_number_of_coe_events(self):
        """
        :return: The mean number of coexistence events, i.e. intersections, of
         all Xgraph objects in this list.
        """
        return self.__avrg_numpy_func("avrg_number_of_coe_events")

    def avrg_number_of_red_events(self):
        """
        :return: The mean number of redundancy events, i.e. intersections, of
         all Xgraph objects in this list.
        """
        return self.__avrg_numpy_func("avrg_number_of_red_events")

    def avrg_diameter(self):
        """
        :return: The average diameter of all Xgraph objects in this list.
        """
        return self.__avrg_numpy_func("avrg_diameter")

    def avrg_node_degree(self):
        """
        :return: The mean node degree of all Xgraph objects in this list.
         More precisely, it calculates the average node degree of every Xgraph
         in the list, sums them up and calculates their mean.
        """
        return self.__avrg_numpy_func("avrg_node_degree")

    def avrg_node_density(self):
        """
        :return: The mean node density of all Xgraph objects in this list.
         More precisely, it calculates the node density of every Xgraph
         in the list, sums them up and calculates their mean.
        """
        return self.__avrg_numpy_func("avrg_node_density")

    def avrg_node_number(self):
        """
        :return: The mean number of nodes of all Xgraph objects in this list.
        """
        return self.__avrg_numpy_func("avrg_node_number")

    def calc_coexistence(self, k: int=None):
        """
        Calculate coexistence of all Xgraph objects in this list up to a path
        length of k hops.
        If no k is given, calculate coexistence up to the diameter of the
        Xgraph object with the biggest diameter out of all Xgraph objects in
        this list.

        :param k: range of maximum allowed path length to be tested
        """
        self.__calc_geometry(k, False)
        # calc_geometry_xgraphlist_multiprocessing(self, None, False)

    def calc_redundancy(self, k: int=None):
        """
        Calculate redundancy of all Xgraph objects in this list up to a path
        length of k hops.
        If no k is given, calculate redundancy up to the diameter of the
        Xgraph object with the biggest diameter out of all Xgraph objects in
        this list.

        :param k: range of maximum allowed path length to be tested
        """
        self.__calc_geometry(k)

    # TODO description
    def check_for_errors(self):
        error_dict = dict()
        for xg in self:
            if self._small_mem:
                xg = Xgraph.load(self.__xgraph_dir + "/" + str(xg) + ".xg")
            error_dict[str(xg.get_id())] = xg.get_error_flag()
        return error_dict

    def get_coexistence_stats(self, percent=True):
        """
        :param percent: False for total values; otherwise percentage values.
        :return: A dictionary of the mean coexistence data of all Xgraph
         objects in this list.
        """
        return self.__get_geometry_stats(percent=percent, type_=False)

    def get_id(self):
        """
        :return: The ID of this XgraphList object.
        """
        return self._id

    def get_redundancy_stats(self, percent=True):
        """
        :param percent: False for total values; otherwise percentage values.
        :return: A dictionary of the mean redundancy data of all Xgraph
         objects in this list.
        """
        return self.__get_geometry_stats(percent=percent, type_=True)

    def get_small_mem(self):
        return self._small_mem

    def get_xgraph_dir(self):
        if self._small_mem:
            return self.__xgraph_dir
        else:
            print("_small_mem not set")
            return None

    @staticmethod
    def load(file_path: "relative path to file"):
        """
        Load an XgraphList object from a binary file. The file must have been
        saved via XgraphList.save().

        :param file_path: Relative path to the XgraphList file.
        :return: The XgraphList object from the input file.
        """
        with open(file_path, "rb") as f:
            xgl = load(f)

            if xgl.get_small_mem():
                for xg_id in xgl:
                    if not Files.file_exists(xgl.get_xgraph_dir() +
                                                     "/" + str(xg_id) + ".xg"):
                        warnings.warn("Xgraph file " + str(xg_id) +
                                      " was removed from this list because "
                                      "no corresponding file existed.")
                        xgl.remove(xg_id)

            if not xgl.are_properties_homogeneous():
                warnings.warn(
                    "The Xgraph objects in this list do not have homogeneous "
                    "properties!")
            return xgl

    @staticmethod
    def load_dir(small_mem: bool, folder_path: "relative path to file"):
        """
        Create an XgraphList object and load all Xgraph files from a directory
        into it.
        Load an XgraphList object from a binary file. The file must have been
        saved via XgraphList.save().

        :param small_mem: flag for limited memory
        :param folder_path: Relative path to the folder containing only
         XgraphList files. Must contain an "/" at the end.
        :return: An XgraphList object containing all Xgraph objects in the
         given folder.
        """
        xgl = XgraphList(small_mem=small_mem)
        for file in os.listdir(folder_path):
            xg = Xgraph.load(folder_path + file)
            xgl.append(xg)
        return xgl

    def max_diameter(self):
        """
        :return: The diameter of the Xgraph object in this list with the
         biggest diameter out of all Xgraph objects in the list.
        """
        return self.__max_numpy_func("max_diameter")

    def pretty_print(self):
        """
        Prints the ID, number of nodes, edges, redundancy and
        coexistence events of each Xgraph in this list on the console.
        """
        print("XgraphList", self._id)
        print("[")

        for xg in self:
            if self._small_mem:
                print(xg)
            else:
                xg.pretty_print()
        print("]")

    def remove_graph(self, graph):
        """
        Remove a graph from this list.

        :param graph: The Xgraph object (in case of _small_mem: its UUID ID,
         possibly as a string) which is to be removed from the list.
        """
        del_id = str(graph)
        if not self._small_mem:
            del_id = graph.get_id()

        del_graph = []
        for g in self:
            if self._small_mem and str(g) == del_id:
                del_graph.append(g)
            else:
                if g.get_id() == del_id:
                    del_graph.append(g)
        if del_graph:
            [self.remove(g) for g in del_graph]
            self.__reset_values()

    def save(self, file_path: "relative path to safe file"):
        # TODO small_mem: save Xgraphs in a sub directory
        """
        Save an XgraphList object to a binary file. The file can be loaded via
        XgraphList.load().

        :param file_path: Relative path to the XgraphList file.
        """
        # cut off file name from path and make each directory in it (if needed)
        Files.ensure_dir(os.path.dirname(file_path))
        with open(file_path, "wb") as f:
            dump(self, f, 4)  # pickle.dump
        if self._small_mem:
            for xg_id in self:
                xg = Xgraph.load(self.__xgraph_dir + "/" + str(xg_id) + ".xg")
                xg.save(os.path.dirname(file_path) + "/"
                        "Xgraph_files_" + str(self.get_id()) + "/" +
                        str(xg_id) + ".xg")



# with help from
# http://stackoverflow.com/questions/14305337/animated-networkx-graph-in-tkcanvas-background-color
class XgraphGUI:
    """
    A GUI for Xgraph objects with already calculated coexistence and
    redundancy values.
    All coexistence and redundancy events (node in triangle, intersecting
    edges) are displayed in lists and the corresponding nodes and edges in the
    graph representation will be highlighted. Red color indicates the actual
    nodes and edges of the corresponding event, while cyan color indicates
    edges which connect two nodes as a shortest path.
    """

    def __init__(self, xg: Xgraph):
        self.graph = xg.get_graph()
        self.red_dict = xg.get_redundancy()._dict
        self.red_events = list(self.red_dict.keys())
        self.coe_dict = xg.get_coexistence()._dict
        self.coe_events = list(self.coe_dict.keys())

        # setting up the GUI frame
        self.root = Tk.Tk()
        self.root.wm_title("Window title")
        self.root.wm_protocol('WM_DELETE_WINDOW', self.root.quit)

        # set background color
        self.figure = plt.figure(figsize=(5, 4))
        # self.figure.set_facecolor('w')

        # get subplot
        self.subplot = self.figure.add_subplot(111)
        plt.axis('off')

        nx.draw_networkx(self.graph,
                         pos=nx.get_node_attributes(self.graph, "pos"),
                         ax=self.subplot, alpha=1, node_color="black",
                         node_size=10, with_labels=False)
        # xlim = a.get_xlim()
        # ylim = a.get_ylim()

        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.root)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self.root)
        toolbar.update()
        toolbar.pack(side=Tk.TOP, fill=Tk.X)
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        # Listbox with redundancy events + scrollbar (order is important)
        self.lb1 = Tk.Listbox(self.root, width=70)
        for event in self.red_events:
            self.lb1.insert(Tk.END, (event, first_one(self.red_dict[event])))
        self.lb1.bind("<<ListboxSelect>>", self.highlight_redundancy)

        scroll1 = Tk.Scrollbar(self.root, command=self.lb1.yview)
        self.lb1.config(yscrollcommand=scroll1.set)

        self.lb1.pack(side="left")
        scroll1.pack(side="left", fill=Tk.Y)

        # Listbox with coexistence events + scrollbar (order is important)
        self.lb2 = Tk.Listbox(self.root, width=70)
        for event in self.coe_events:
            self.lb2.insert(Tk.END, (event, first_one(self.coe_dict[event])))
        self.lb2.bind("<<ListboxSelect>>", self.highlight_coexistence)

        scroll2 = Tk.Scrollbar(self.root, command=self.lb2.yview)
        self.lb2.config(yscrollcommand=scroll2.set)

        scroll2.pack(side="right", fill=Tk.Y)
        self.lb2.pack(side="right")

        Tk.mainloop()

    def highlight_coexistence(self, _):
        event = self.coe_events[self.lb2.curselection()[0]]

        triangle = event[0]
        point1 = triangle[0]
        point2 = triangle[1]
        point3 = triangle[2]
        node = event[1]
        node_list1 = [point1, point2, point3, node]

        path1_nodes = nx.shortest_path(self.graph, point1, node)
        path2_nodes = nx.shortest_path(self.graph, point2, node)
        path3_nodes = nx.shortest_path(self.graph, point3, node)
        paths = [path1_nodes, path2_nodes, path3_nodes]
        node_list2 = (set(path1_nodes) | set(path2_nodes) | set(path3_nodes) -
                      set(node_list1))

        edge1 = (point1, point2)
        edge3 = (point1, point3)
        edge2 = (point2, point3)
        edge_list1 = [edge1, edge2, edge3]

        edge_list2 = []
        for path in paths:
            edge_list2.extend(node_list_to_edges(path, self.graph))

        self.highlight(node_list1, edge_list1, node_list2, edge_list2)

    def highlight_redundancy(self, _):
        event = self.red_events[self.lb1.curselection()[0]]

        point1 = event[0][0]
        point2 = event[0][1]
        point3 = event[1][0]
        point4 = event[1][1]
        node_list1 = [point1, point2, point3, point4]

        edge1 = (point1, point2)
        edge2 = (point3, point4)
        edge_list1 = [edge1, edge2]

        path_1_nodes = nx.shortest_path(self.graph, point1, point3)
        path_2_nodes = nx.shortest_path(self.graph, point1, point4)
        path_3_nodes = nx.shortest_path(self.graph, point2, point3)
        path_4_nodes = nx.shortest_path(self.graph, point2, point4)
        paths = [path_1_nodes, path_2_nodes, path_3_nodes,
                 path_4_nodes]
        node_list2 = (set(path_1_nodes) | set(path_2_nodes) |
                      set(path_3_nodes) | set(path_4_nodes)) - set(node_list1)

        edge_list2 = []
        for path in paths:
            edge_list2.extend(node_list_to_edges(path, self.graph))

        self.highlight(node_list1, edge_list1, node_list2, edge_list2)

    def highlight(self, node_list1, edge_list1,
                  node_list2=None, edge_list2=None):
        self.subplot.cla()

        uncolored_edges = list(
            set(self.graph.edges()) - set(edge_list1) - set(edge_list2))
        self.draw_edges(uncolored_edges)
        if edge_list2:
            self.draw_edges(edge_list2, "cyan", 2)
        self.draw_edges(edge_list1, "red", 2)

        uncolored_nodes = list(
            set(self.graph.node) - set(node_list1) - set(node_list2))
        self.draw_nodes(uncolored_nodes)
        if node_list2:
            self.draw_nodes(node_list2, "cyan", 15)
        self.draw_nodes(node_list1, "red", 15)

        # nx.draw_networkx_labels(self.graph, pos=pos)

        # a.set_xlim(xlim)
        # a.set_ylim(ylim)
        plt.axis("off")
        self.canvas.draw()

    def draw_edges(self, edge_list, color=None, width=None):
        if not color:
            color = "black"

        if not width:
            width = 1

        pos = nx.get_node_attributes(self.graph, "pos")

        nx.draw_networkx_edges(self.graph, pos=pos, edgelist=edge_list,
                               edge_color=color, width=width, ax=self.subplot,
                               with_labels=False)

    def draw_nodes(self, node_list, color=None, size=None):
        if not color:
            color = "black"

        if not size:
            size = 10

        pos = nx.get_node_attributes(self.graph, "pos")

        nx.draw_networkx_nodes(self.graph, pos=pos, nodelist=node_list,
                               node_color=color, ax=self.subplot,
                               alpha=1, node_size=size, with_labels=False)


def worker_(xg: Xgraph, k: int=None, type_: bool=True):
    """
    Worker function for multiprocessing stuff. Calculates
    coexistence/redundancy of an Xgraph's underlying graph and writes it to a
    binary file with the id of the Xgraph (".coe"/".red") as its filename.
    The file is written to the current directory.

    :param k: range of maximum allowed path length to be tested
    :param xg: Xgraph object
    :param type_: if True, calculate redundancy; else coexistence
    """
    if k is None:
        k = xg.diameter()

    warnings.filterwarnings("ignore")  # ignore the get_graph() warning

    if type_:
        g = redundancy.Redundancy(k, xg.get_graph())
        file_extension = ".red"
    else:
        g = coexistence.Coexistence(k, xg.get_graph())
        file_extension = ".coe"

    warnings.resetwarnings()  # enable warnings again

    with open(str(xg.get_id()) + file_extension, "wb") as f:
        dump(g, f)


def calc_geometry_xgraphlist_multiprocessing(
        xgl: XgraphList, k: int=None, type_: bool=True):
    """
    Calculates coexistence or redundancy of each Xgraph object in the given
    XgraphList object in multiprocessing mode.

    At any given time, each core of the processor calculates coexistence or
    redundancy of one Xgraph object of the XgraphList object (if there are
    any); after calculation, the Xgraph object is saved as a binary file
    (including its respective coexistence or redundancy object; see worker()
    for more information). When all calculations have terminated, each Xgraph
    object is loaded from its binary file, assigned to its respective Xgraph in
    the XgraphList object, and the binary file is deleted.

    :param xgl: XgraphList object
    :param k: range of maximum allowed path length to be tested
    :param type_: if True, calculate redundancy; else coexistence
    :return:
    """
    file_extension = ".coe"
    if type_:
        file_extension = ".red"

    pool = Pool()
    for xg in xgl:
        if xgl.get_small_mem():
            xg = Xgraph.load(xgl.get_xgraph_dir() + "/" + str(xg) + ".xg")
        r = pool.apply_async(worker_, args=(xg, k, type_))
        r.get()
    pool.close()
    pool.join()

    for xg in xgl:
        if xgl.get_small_mem():
            xg = Xgraph.load(xgl.get_xgraph_dir() + "/" + str(xg) + ".xg")
        file_name = str(xg.get_id()) + file_extension
        with open(file_name, "rb") as f:
            if type_:
                xg._red = load(f)
            else:
                xg._coe = load(f)
            if xgl.get_small_mem():
                xg.save(xgl.get_xgraph_dir() + "/" + str(xg.get_id()) + ".xg")
        os.remove(file_name)

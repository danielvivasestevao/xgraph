import networkx as nx
from shapely.geometry import LineString


def intersect(line1, line2, graph, node_pos_dict=None):
    """
    Checks whether two line segments defined by the positions of four nodes
    intersect.

    :param line1: tuple of two node IDs
    :param line2: tuple of two node IDs
    :param graph: networkx.Graph
    :param node_pos_dict: a dictionary containing all nodes of a graph as keys
    and their respective positions as values; used for reasons of
    efficiency when calculating a big amount of intersections (>10000)
    :return: boolean (intersection (True) or no intersection (False))
    """
    try:
        node0_0_pos = node_pos_dict[line1[0]]
    except TypeError:
        node_pos_dict = nx.get_node_attributes(graph, 'pos')
        node0_0_pos = node_pos_dict[line1[0]]
    node0_1_pos = node_pos_dict[line1[1]]
    node1_0_pos = node_pos_dict[line2[0]]
    node1_1_pos = node_pos_dict[line2[1]]
    linestring = LineString
    """ TODO generate a set of all linestrings instead of converting always two
            lines on-the-fly"""
    l1 = linestring([node0_0_pos, node0_1_pos])
    l2 = linestring([node1_0_pos, node1_1_pos])
    result = l1.intersects(l2)
    if (result and line1[0] == line2[0] or line1[0] == line2[1] or
            line1[1] == line2[0] or line1[1] == line2[1]):
        return False
    # intersection = l1.intersection(l2)
    # intersection.x / intersection.y
    return result


def get_node_data(edge, node_pos_dict):
    nd0 = NodeData(0, edge, node_pos_dict)
    nd1 = NodeData(1, edge, node_pos_dict)
    if nd0.node is nd1.node:
        raise Exception("One edge returns equal nodes for both ends.")
    return nd0, nd1


def get_intersections_dict(graph=None, intersections: [((int, int), (int, int))]=None):
    _dict = {}
    if not intersections:
        if not graph:
            raise Exception("get_intersections_dict needs either a "
                            "graph or a list of intersections")
        intersections = simple_bentley_ottmann(graph)
    for e in intersections:
        _dict[e] = {'strong': {}, 'weak': {}}
    return _dict


def simple_bentley_ottmann(graph):
    """
    A simplified version of the Bentley-Ottman sweep line algorithm.

    This version iterates through the whole graph from left to right (low x
    value to high x value), keeps track of which edges have been started (their
    leftmost point has been seen) and checks if a new edge intersects with any
    other edge that has not been closed (their rightmost point has not been
    reached yet) yet.

    :param graph: networkx.Graph
    :return: list of edge tuples
    """
    # get all nodes of the graph
    point_id_2_pos = nx.get_node_attributes(graph, 'pos')

    nd_list = []  # list of NodeData

    for e in graph.edges():
        nd = get_node_data(e, point_id_2_pos)
        nd_list.append(nd[0])
        nd_list.append(nd[1])

    # sort list by x values
    nd_list.sort(key=lambda x: point_id_2_pos[x.node][0])

    active_edges = set()
    intersections = set()

    for nd in nd_list:
        if nd.side is 0:
            e1 = nd.edge
            for e2 in active_edges:
                if intersect(e1, e2, graph, point_id_2_pos):
                    if (e2, e1) not in intersections:
                        intersections.add((e1, e2))
            active_edges.add(e1)
        elif nd.side is 1:
            active_edges.remove(nd.edge)

    return list(intersections)


class NodeData:
    """
    Data structure which keeps track of an edge's "start" (left) and "end"
    (right) node.
    """
    def __init__(self, side, edge, node_pos_dict):
        """
        :param side: 0 for left, 1 for right
        :param edge: an edge defined by a tuple of nodes, defined by their IDs
        :param node_pos_dict: a dictionary which maps a nodes' ID to its
          position (tuple of floats)
        """
        self.side = side
        self.edge = edge

        node_0_x_pos = node_pos_dict[edge[0]][0]
        node_1_x_pos = node_pos_dict[edge[1]][0]

        if node_0_x_pos <= node_1_x_pos:
            if side is 0:
                self.node = edge[0]
            elif side is 1:
                self.node = edge[1]
        elif node_0_x_pos > node_1_x_pos:
            if side is 0:
                self.node = edge[1]
            elif side is 1:
                self.node = edge[0]

import networkx as nx


class Triangle:
    """
    Data structure which keeps track of a triangle's leftmost and rightmost
    node. NOTE: This class does NOT check whether the triangle actually
    exists (in a graph).
    """
    def __init__(self, node1, node2, node3, node_pos_dict):
        """
        :param node1: first node of the triangle (given by node ID)
        :param node2: second node of the triangle (given by node ID)
        :param node3: third node of the triangle (given by node ID)
        :param node_pos_dict: dictionary mapping a node ID to the node's
          position in the graph
        """
        l = [node1, node2, node3]
        # sort by node's x value and break ties with the node's ID
        l.sort(key=lambda x: (node_pos_dict[x][0], x))
        self.node1 = l[0]  # leftmost node
        self.node2 = l[1]
        self.node3 = l[2]  # rightmost node
        self.node_pos_dict = node_pos_dict  # creates a reference, not a copy

    def __eq__(self, other):
        if (isinstance(other, self.__class__) and
                    self.node1 == other.node1 and
                    self.node2 == other.node2 and
                    self.node3 == other.node3 and
                    self.node_pos_dict[self.node1] ==
                    other.node_pos_dict[other.node1] and
                    self.node_pos_dict[self.node2] ==
                    other.node_pos_dict[other.node2] and
                    self.node_pos_dict[self.node3] ==
                    other.node_pos_dict[other.node3]):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # equivalent items should hash the same! Make it __eq__ compatible
        key1 = self.get_triangle()
        key2 = (self.node_pos_dict[self.node1],
                self.node_pos_dict[self.node2],
                self.node_pos_dict[self.node3])
        return hash(key1) + hash(key2)

    def __lt__(self, other):
        return self.node_pos_dict[self.node1][0] <\
               other.node_pos_dict[other.node1][0]

    def __le__(self, other):
        return self.node_pos_dict[self.node1][0] <=\
               other.node_pos_dict[other.node1][0]

    def __ge__(self, other):
        return self.node_pos_dict[self.node1][0] >=\
               other.node_pos_dict[other.node1][0]

    def __gt__(self, other):
        return self.node_pos_dict[self.node1][0] >\
               other.node_pos_dict[other.node1][0]

    def __str__(self):
        return "T(" + str(self.node1) + ", " + str(self.node2) + ", " + \
            str(self.node3) + ")"

    def contains_node(self, node):
        return node_in_triangle(node, self.get_triangle(),
                                node_pos_dict=self.node_pos_dict)

    def does_exist(self, graph: nx.Graph):
        return triangle_exists(graph, self.get_triangle())

    def get_leftmost_node(self):
        return self.node1

    def get_middle_node(self):
        return self.node2

    def get_rightmost_node(self):
        return self.node3

    def get_triangle(self):
        return self.node1, self.node2, self.node3


def get_triangles(graph, as_triangle_object=False, node_pos_dict=None):
    """
    http://www.cl.cam.ac.uk/~cm542/teaching/2010/stna-pdfs/stna-lecture8.pdf
    page 44

    Returns a generator with all triangles in the graph.
    A triangle is a node triple (n1, n2, n3) where
    -all three nodes are interconnected via edges in the graph
    -all three nodes are not colinear
    :param graph: nx.Graph object
    :param as_triangle_object: True if the triangles are to be yielded in form
      of Triangle objects; otherwise triples of node IDs are returned
    :param node_pos_dict: a dictionary mapping each node's ID to its position
      in form of a tuple (x_coord, y_coord)
    :return: all triangles of the given graph
    """
    if not node_pos_dict:
        node_pos_dict = nx.get_node_attributes(graph, "pos")
    for n1 in graph.nodes_iter():
        neighbors_n1 = set(graph[n1])
        for n2 in filter(lambda x: x > n1, graph.nodes_iter()):
            neighbors_n2 = set(graph[n2])
            common_neighbors = neighbors_n1 & neighbors_n2
            for n3 in filter(lambda x: x > n2, common_neighbors):
                if (triangle_exists(graph, (n1, n2, n3)) and
                      not is_triangle_flat((n1, n2, n3), node_pos_dict)):
                    if as_triangle_object and node_pos_dict:
                        yield Triangle(n1, n2, n3, node_pos_dict)
                    else:
                        yield n1, n2, n3


def get_triangle_list(graph, as_triangle_object=False, node_pos_dict=None):
    return list(get_triangles(graph, as_triangle_object, node_pos_dict))


def triangle_exists(graph, triangle):
    """
    Checks whether three nodes are connected via edges to a triangle in a
    graph.

    :param graph: a networkx graph
    :param triangle: a triple of nodes, e.g. (1, 2, 3)
    :return: a boolean indicating whether the triangle exists in form of nodes
    and edges in the graph
    """
    e1 = graph.has_edge(triangle[0], triangle[1])
    e2 = graph.has_edge(triangle[1], triangle[2])
    e3 = graph.has_edge(triangle[2], triangle[0])
    return e1 and e2 and e3


def is_triangle_flat(triangle, node_pos_dict):
    """
    Checks if a triangle is flat, i.e. whether all three nodes of the triangle
    are collinear, i.e. lying on a line.
    :param triangle: triangle in form of a node triple
    :param node_pos_dict: a dictionary mapping node IDs to their positions in
      form of a tuple (x_coord, y_coord)
    :return:
    """
    n1_pos = node_pos_dict[triangle[0]]
    n2_pos = node_pos_dict[triangle[1]]
    n3_pos = node_pos_dict[triangle[2]]
    return are_colinear(n1_pos, n2_pos, n3_pos)


def are_colinear (pos1, pos2, pos3):
    for i in range(2):
        if pos1[i] == pos2[i] == pos3[i]:
            return True
    return False


def point_in_triangle(vertex_pos_1, vertex_pos_2, vertex_pos_3, point_pos):
    """
    https://totologic.blogspot.de/2014/01/accurate-point-in-triangle-test.html
    Check if a point lies within a triangle using barycentric coordinates.
    :param vertex_pos_1: first point of the triangle in form of a tuple (x,y)
    :param vertex_pos_2: second point of the triangle in form of a tuple (x,y)
    :param vertex_pos_3: third point of the triangle in form of a tuple (x,y)
    :param point_pos: point we are trying to locate in form of a tuple (x,y)
    :return: True if point_pos lies inside of the triangle; False otherwise
    """
    # if the node lies exactly on an edge of the triangle, it is disregarded
    for i in range(2):
        if (are_colinear(vertex_pos_1, vertex_pos_2, point_pos) or
             are_colinear(vertex_pos_2, vertex_pos_3, point_pos) or
             are_colinear(vertex_pos_1, vertex_pos_3, point_pos)):
            return False

    denominator = (vertex_pos_2[1] - vertex_pos_3[1]) *\
                  (vertex_pos_1[0] - vertex_pos_3[0]) +\
                  (vertex_pos_3[0] - vertex_pos_2[0]) *\
                  (vertex_pos_1[1] - vertex_pos_3[1])
    """
    The denominator is zero if all three nodes of the triangle lie on one line.
    In that case the area of the triangle is zero, and thus cannot have a node
    within itself.
    """
    if denominator == 0:
        return False
    a = ((vertex_pos_2[1] - vertex_pos_3[1]) *
         (point_pos[0] - vertex_pos_3[0]) +
         (vertex_pos_3[0] - vertex_pos_2[0]) *
         (point_pos[1] - vertex_pos_3[1])) / denominator
    b = ((vertex_pos_3[1] - vertex_pos_1[1]) *
         (point_pos[0] - vertex_pos_3[0]) +
         (vertex_pos_1[0] - vertex_pos_3[0]) *
         (point_pos[1] - vertex_pos_3[1])) / denominator
    c = 1 - a - b

    return 0 <= a <= 1 and 0 <= b <= 1 and 0 <= c <= 1


def node_in_triangle(node, triangle, graph=None, node_pos_dict=None):
    """
    :param node: a node, defined by its ID
    :param triangle: tuple (triple) of three nodes, defined by their IDs
    :param graph: networkx.Graph
    :param node_pos_dict: a dictionary mapping a node's ID to its position in
      form of a float tuple (x_coord, y_coord)
    :return: true if the position of node lies within the triangle, false
      otherwise
    """
    if not graph and not node_pos_dict:
        raise Exception("No graph or node-position dictionary given")
    # Nodes that are vertices of the triangle don't lie inside the triangle.
    if node in triangle:
        return False
    try:
        vertex1 = node_pos_dict[triangle[0]]
    except TypeError:
        node_pos_dict = nx.get_node_attributes(graph, 'pos')
        vertex1 = node_pos_dict[triangle[0]]
    vertex2 = node_pos_dict[triangle[1]]
    vertex3 = node_pos_dict[triangle[2]]
    node_pos = node_pos_dict[node]

    return point_in_triangle(vertex1, vertex2, vertex3, node_pos)


def get_nodes_in_triangles_dict(graph=None, coe_event_list=None):
    """
    Returns a dictionary with a key (triangle, node) for every event of a node
    of the graph lying within the triangle in the graph.
    The values are all empty coexistence dictionaries:
    {'strong': {}, 'weak': {}}
    :param graph: networkx.Graph
    :return: empty coexistence dictionaries
    """
    if not graph and not coe_event_list:
        raise Exception("get_nodes_in_triangles_dict needs either a "
                "graph or a list of coexistence events")

    _dict = {}

    if coe_event_list:
        for e in coe_event_list:
            _dict[e] = {'strong': {}, 'weak': {}}
        return _dict

    node_pos_dict = nx.get_node_attributes(graph, 'pos')

    # get a list of tuples (node_id, (x_coord, y_coord)), sorted by x_coord
    # and break ties with y_coord
    node_position_tuple_list = list(node_pos_dict.items())
    node_position_tuple_list.sort(key=lambda x: (x[1][0], x[1][1]))
    # get a list of only node IDs sorted by their x_coord
    node_list_sorted = [x[0] for x in node_position_tuple_list]

    # get a list of all triangles and sort it by the x_coord of the leftmost
    # node of each triangle
    triangles = get_triangle_list(
        graph, as_triangle_object=True, node_pos_dict=node_pos_dict)
    triangles.sort()

    for node in node_list_sorted:
        for triangle in triangles:
            if (node_pos_dict[triangle.get_rightmost_node()][0] <
                 node_pos_dict[node][0]):
                continue
            if (node_pos_dict[triangle.get_leftmost_node()][0] >
                 node_pos_dict[node][0]):
                break
            if node in triangle.get_triangle():
                continue
            if triangle.contains_node(node):
                _dict[(triangle.get_triangle(), node)] =\
                    {'strong': {}, 'weak': {}}

    return _dict

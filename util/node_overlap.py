from Graphs.Xgraph import *
from util.math_funcs import isclose
from operator import itemgetter
from util.ProcessTimer import ProcessTimer

"""
This code has been written to check if a bug from phase1 still exists.
When generating graphs in multiprocessing mode with the code from phase1,
some graphs generated at the same time would sometimes have nodes with the same
position data. The nodes would thus 'overlap' if you layed plots of the graphs
upon each other.
"""


def sort_nodes(node_dict):
    """
    :param node_dict: dictionary of nodes in the form of
     {node_id: {'pos': (x_value, y_value)}, ...}
     with node_id: int, x_value and y_value: float
    :return: a list of node IDs sorted by x-values first and y-values second
    """
    # convert the node_dict to a list of triples: (node_id, x_val, y_val)
    node_list = [(node_id,
                  node_dict[node_id]['pos'][0],  # x-value
                  node_dict[node_id]['pos'][1])  # y-value
                 for node_id in node_dict]
    return sorted(node_list, key=itemgetter(1, 2))


def node_overlap(xgl: XgraphList):
    """
    Checks whether Xgraph objects in an XgraphList have nodes with the same
    position information.

    :param xgl: XgraphList object
    :return: dictionary with the number of overlapping nodes of two graphs with
     those graphs's IDs in tuple form as the dictionary key
    """
    p = ProcessTimer("calculating node_overlap")
    p.start()

    closeness_dict = dict()
    c = -1

    for xg1 in xgl:
        for xg2 in xgl:
            c += 1
            p2 = ProcessTimer("     calculating overlap no. " + str(c))
            p2.start()

            print("Comparison no.", c)

            if xg1 == xg2:
                continue

            nodes1 = xg1.get_graph().node
            nodes2 = xg2.get_graph().node

            nodes1_s = sort_nodes(nodes1)
            nodes2_s = sort_nodes(nodes2)

            closeness = 0

            for node_triple1 in nodes1_s:
                for node_triple2 in nodes2_s:
                    # if x-val of n2 is too small to be equal, continue
                    if round(node_triple2[1]) < round(node_triple1[1]):
                        continue
                    # if x-val of n2 is too big, the x-vals of all the next
                    # nodes are too big, too; n1 is done
                    if round(node_triple2[1]) > round(node_triple1[1]):
                        break
                    if (isclose(node_triple1[1], node_triple2[1]) and
                            isclose(node_triple1[2], node_triple2[2])):
                        closeness += 1

                    closeness_dict[(xg1.get_id(), xg2.get_id())] = closeness

            p.end()

    p.end()
    return closeness_dict


if __name__ == "__main__":
    xgl = XgraphList()
    file_path = "../Xgraph_testing/strasbourg_2000on2600-1/"
    for i in range(0, 50):
         xg = Xgraph.load(file_path + "graph" + str(i) + ".xg")
         xgl.append(xg)
    with open(file_path + "closeness_test.txt", "w") as f:
        cl_d = node_overlap(xgl)
        for key in cl_d:
            f.write(str(key) + ": " + str(cl_d[key]) + "\n")

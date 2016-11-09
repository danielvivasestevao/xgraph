import re
from Geometry.triangles import Triangle

def edge_list_file2list(file_path_to_edge_list):
    with open(file_path_to_edge_list, "r") as f:
        edge_list = [eval(line.rstrip("\n")) for line in f]
    return edge_list


# [(Edge { left: 93, right: 142 }, Edge { left: 137, right: 97 }), (Edge { left: 93, right: 142 }, Edge { left: 137, right: 98 })]\n
def edge_list_str2list(edge_list_str):
    # Returns the list backwards, which does not really matter,
    # because edges come in tuples and their order within a tuple does not matter.
    numbers = re.compile('\d+(?:\.\d+)?')
    number_list = numbers.findall(edge_list_str)
    edge_list = []
    while number_list:
        edge_list.append(
            ((int(number_list.pop()), int(number_list.pop())),
             (int(number_list.pop()), int(number_list.pop())))
        )
    return edge_list


def coe_event_list_str2list(coe_event_list_str, node_pos_dict):
    numbers = re.compile('\d+(?:\.\d+)?')
    number_list = numbers.findall(coe_event_list_str)
    event_list = []
    while number_list:
        t3 = int(number_list.pop())
        t2 = int(number_list.pop())
        t1 = int(number_list.pop())
        n = int(number_list.pop())
        event_list.append((Triangle(t1, t2, t3, node_pos_dict).get_triangle(), n))
    return event_list


def triangle_list_str2list(triangle_list_str, node_pos_dict):
    numbers = re.compile('\d+(?:\.\d+)?')
    number_list = numbers.findall(triangle_list_str)
    triangle_list = []
    while number_list:
        t3 = int(number_list.pop())
        t2 = int(number_list.pop())
        t1 = int(number_list.pop())
        triangle_list.append(Triangle(t1, t2, t3, node_pos_dict))
    return triangle_list

# s = "[(Edge { left: 93, right: 142 }, Edge { left: 137, right: 97 }), (Edge { left: 93, right: 142 }, Edge { left: 137, right: 98 }),]\n"
# print(edge_list_str2list(s))
# t = "[(45, Triangle { left: 16, middle: 196, right: 24 }), (45, Triangle { left: 16, middle: 196, right: 24 })"
# print(coe_event_list_str2list(t))
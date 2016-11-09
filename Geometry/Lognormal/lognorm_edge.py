from math import radians, cos, sin, asin, sqrt, log
import numpy
import random


def get_distance(pos1, pos2):
    pos_a_array = numpy.array(pos1)
    pos_b_array = numpy.array(pos2)
    return numpy.linalg.norm(pos_a_array - pos_b_array)


def get_distance_gps(pos1, pos2):
    # convert value from haversine formula from kilometer to meter
    return haversine(pos1[0], pos1[1], pos2[0], pos2[1]) * 1000


# http://stackoverflow.com/questions/4913349/
# haversine-formula-in-python-bearing-and-distance-between-two-gps-points
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    Tested with distance Girod-Montabaur (6.18 kilometers):

    distance = lognormal_graph.haversine(7.909531000000015, 50.45158379999999,
                                         7.825795299999982, 50.4358385)
    print(distance)  # 6.182644268597313
    print(lognormal_graph.kilometer_to_meter(distance))  # 6182.644268597313
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles
    return c * r


# given the distance, it calculates the path loss according to
# deterministic channel model
# (friis equation extended by path loss exponent)
def distance_to_attenuation_db(path_loss_onemeter,
                               path_loss_exponent,
                               distance):
    return path_loss_onemeter + 10 * path_loss_exponent * log(distance, 10)


def lognorm_edge(path_loss_exponent: float, variance: float,
                 transmit_power: float, path_loss_one_meter: float,
                 pos1: (float, float), pos2: (float, float),
                 gps: bool=False):
    stddev = sqrt(variance)

    calc_distance = get_distance
    if gps:
        calc_distance = get_distance_gps

    lognorm = random.normalvariate(0, stddev)
    attenuation_db = lognorm + distance_to_attenuation_db(
        path_loss_one_meter, path_loss_exponent, calc_distance(pos1, pos2))
    rssi = transmit_power - attenuation_db

    if rssi >= -95:
        return rssi
    return False

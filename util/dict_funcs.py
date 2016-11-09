import numpy as np


def first_one(dict_):
    """
    Takes a dictionary of an event in a coexistence or redundancy object and
    returns for 'weak' and 'strong' the lowest k value which fulfills it.

    Example:
        dict_ = { 'weak': {1: 0, 2: 1, 3:1}, 'strong': {1:0, 2:0, 3:1} }
        first_one(dict_)  # returns {'weak': 2, 'strong': 3}
    """

    return_dict = {}
    for key in dict_:  # 'weak' or 'strong'
        for k in range(1, len(dict_[key])):  # from lowest to highest k
            if dict_[key][k] == 1:
                return_dict[key] = k
                break
    return return_dict


def dict_mean(dict_list):
    """
    Takes a list of dictionaries with similar key values and returns one
    dictionary that has for each key the mean of the other dictionaries'
    values of the respective key.
    """
    mean_dict = {}
    for key in dict_list[0]:
        value_list = []
        for dict_ in dict_list:
            value_list.append(dict_[key])
        mean_dict[key] = np.array(value_list).mean()
    return mean_dict


def stats_mean(stats_list):
    """
    Takes a list of dictionaries of redundancy/coexistence stats and returns
    one dictionary that has for each key the mean of the other dictionaries'
    values of the respective key.
    """

    # find the dictionary with the highest k
    biggest_dict = stats_list[0]
    for stats_dict in stats_list:
        if len(stats_dict["weak"]) > len(biggest_dict["weak"]):
            biggest_dict = stats_dict

    mean_dict = {"strong": {}, "weak": {}}
    for strength in biggest_dict:
        for key in biggest_dict[strength]:
            value_list = []
            for stats_dict in stats_list:
                try:
                    value_list.append(stats_dict[strength][key])
                except KeyError:
                    value_list.append(1)
            mean_dict[strength][key] = np.array(value_list).mean()

    return mean_dict

def get_containing_elements(l: list, ele):
    """
    :param l: a list of collections (tuples, lists, sets, dicts)
    :param ele: an element
    :return: a list of all collections in l which contain ele
    """
    for coll in l:
        if ele in coll:
            yield coll

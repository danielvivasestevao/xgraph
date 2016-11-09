# http://stackoverflow.com/a/32160518/4466895

import itertools
from shapely.ops import izip


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return izip(a, b)

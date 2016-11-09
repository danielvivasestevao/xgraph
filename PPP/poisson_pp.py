import numpy as np


def spatial_poisson_gen(lam, dim_x, dim_y=None):
    """Returns an array of tuples (coordinates) according to a Poisson
    point process (spatial Poisson process).

    :param lam: lambda
    :param dim_x: horizontal dimension of the field
    :param dim_y: vertical dimension of the filed
    :return:
    """
    n = np.random.poisson(lam)  # warning does not make sense
    x = np.random.rand(n)
    y = np.random.rand(n)
    if dim_y is None:
        dim_y = dim_x
    for i in range(n):
        x[i] *= dim_x
        y[i] *= dim_y
    return np.vstack(([x.T], [y.T])).T

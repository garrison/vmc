import operator

def product(iterable):
    return reduce(operator.mul, iterable)

def add_hc(x):
    """Returns the value plus its complex conjugate"""
    return 2 * x.real

def ensure_real(x):
    """Returns a real number as a real data type

    >>> ensure_real(23+0j)
    23.0
    """
    assert x.imag == 0
    return x.real

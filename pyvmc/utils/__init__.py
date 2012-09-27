import operator

def product(iterable):
    return reduce(operator.mul, iterable)

def add_hc(x):
    """Returns the value plus its complex conjugate"""
    return 2 * x.real

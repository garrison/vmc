import operator

def product(iterable):
    return reduce(operator.mul, iterable)

import operator

def product(iterable):
    return reduce(operator.mul, iterable)

def kronecker(a, b):
    return 1 if (a == b) else 0

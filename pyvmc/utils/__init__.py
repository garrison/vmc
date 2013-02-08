import operator
from collections import Sequence

import numpy

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

def average(seq):
    assert isinstance(seq, (Sequence, numpy.ndarray))
    length = len(seq)
    assert length > 0
    return sum(seq) / float(length)

def stddevmean(seq):
    assert isinstance(seq, (Sequence, numpy.ndarray))
    if len(seq) > 1:
        x = average(seq)
        variance = sum(abs(x - b) ** 2 for b in seq) / (len(seq) - 1)
        return numpy.sqrt(variance / len(seq))
    else:
        return 0.0

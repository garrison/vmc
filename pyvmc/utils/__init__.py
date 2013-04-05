import operator
from collections import Sequence
from functools import reduce
from math import isinf

import numpy

def ensure_real(x):
    """Returns a real number as a real data type

    >>> ensure_real(23+0j)
    23.0
    """
    assert x.imag == 0
    return x.real

def stddevmean(seq):
    assert isinstance(seq, (Sequence, numpy.ndarray))
    if len(seq) > 1:
        x = numpy.average(seq)
        variance = sum(abs(x - b) ** 2 for b in seq) / (len(seq) - 1)
        return numpy.sqrt(variance / len(seq))
    else:
        return 0.0

def average_and_stddevmean(seq):
    return numpy.average(seq), stddevmean(seq)

def is_square_matrix(mat, N=None):
    assert isinstance(mat, numpy.ndarray)
    if N is None:
        return mat.ndim == 2 and mat.shape[0] == mat.shape[1]
    else:
        return mat.ndim == 2 and mat.shape[0] == mat.shape[1] == N

# these are deprecated, as they already exist in numpy
product = numpy.product
average = numpy.average
is_finite = numpy.isfinite

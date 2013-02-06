"""Tools for iterating over a set of parameters (e.g. relative energy scales)
"""

from collections import OrderedDict
from itertools import product
from math import sqrt

sqrt2 = sqrt(2)
inf = float('inf')

def sqrt2_exp(v):
    """Raise the square root of 2 to a power

    This function exists because if we do this the naive way we get additional
    roundoff error."""
    # assume v is integer or negative infinity
    if v == -inf:
        return 0.0
    return 2 ** (v // 2) * (sqrt2 ** (v % 2))

def iterate_values(n_parameters, n_steps=None, exp_func=None):
    if n_steps is None:
        n_steps = 4
    if exp_func is None:
        exp_func = sqrt2_exp

    return (tuple(exp_func(b) for b in a)
            for a in product([-inf] + range(n_steps - 1), repeat=n_parameters)
            if 0 in a)

def iterate_parameters(parameters, n_steps=None, exp_func=None):
    """Iterate over some set of dimensionful parameters.

    This assumes that the states are equivalent if rescaled by a common factor.
    For instance, we don't consider the state where all values are equal to 2,
    as this is equivalent to the state where all values are equal to 1.

    This also skips the state where all values are zero.

    A dict is returned for each state.  The keys of this dictionary are taken
    from the list of `parameters` given.

    >>> for parameters in iterate_parameters(['t', 'J'], n_steps=5):
    ...  print('t = {t:.3f} ; J = {J:.3f}'.format(**parameters))
    ...
    t = 0.000 ; J = 1.000
    t = 1.000 ; J = 0.000
    t = 1.000 ; J = 1.000
    t = 1.000 ; J = 1.414
    t = 1.000 ; J = 2.000
    t = 1.000 ; J = 2.828
    t = 1.414 ; J = 1.000
    t = 2.000 ; J = 1.000
    t = 2.828 ; J = 1.000
    """
    return (OrderedDict(zip(parameters, values))
            for values in iterate_values(len(parameters), n_steps, exp_func))

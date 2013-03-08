"""Tools for iterating over a set of parameters (e.g. relative energy scales)
"""

from __future__ import division

from collections import namedtuple
from itertools import product

inf = float('inf')

def nroot_exp(n):
    """Return a function taises the nth root of 2 to a power

    This function exists because if we do this the naive way we get additional
    roundoff error."""
    nroot = 2 ** (1 / n)

    def f(v):
        # assume v is integer or negative infinity
        if v == -inf:
            return 0.0
        return 2 ** (v // n) * (nroot ** (v % n))

    return f

def iterate_parameters(parameters, n_steps=4, exp_func=nroot_exp(2)):
    """Iterate over some set of dimensionful parameters.

    This assumes that the states are equivalent if rescaled by a common factor.
    For instance, we don't consider the state where all values are equal to 2,
    as this is equivalent to the state where all values are equal to 1.

    This also skips the state where all values are zero.

    A dict is returned for each state.  The keys of this dictionary are taken
    from the list of `parameters` given.

    >>> for parameters in iterate_parameters(['t', 'J'], n_steps=5):
    ...  print('t = {t:.3f} ; J = {J:.3f}'.format(**parameters._asdict()))
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
    ParameterSet = namedtuple("ParameterSet", parameters, rename=True)

    return (ParameterSet._make(exp_func(b) for b in a)
            for a in product([-inf] + range(n_steps - 1), repeat=len(parameters))
            if 0 in a)

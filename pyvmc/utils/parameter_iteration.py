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
    return (dict(zip(parameters, values))
            for values in iterate_values(len(parameters), n_steps, exp_func))

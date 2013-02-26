import json
from functools import partial
from collections import OrderedDict
from fractions import Fraction

from pyvmc.core.boundary_conditions import BoundaryCondition

def tuplize(x):
    if isinstance(x, dict):
        return {k: tuplize(v) for k, v in six.iteritems(x)}
    elif isinstance(x, list):
        return tuple(tuplize(a) for a in x)
    else:
        return x

class CustomEncoder(json.JSONEncoder):
    """This custom JSON encoder can handle complex and Fraction and BoundaryCondition types.
    """

    def default(self, obj):
        if isinstance(obj, complex):
            return OrderedDict([
                ("__class__", "complex"),
                ("real", obj.real),
                ("imag", obj.imag),
            ])
        if isinstance(obj, Fraction):
            return OrderedDict([
                ("__class__", "Fraction"),
                ("numerator", obj.numerator),
                ("denominator", obj.denominator),
            ])
        if isinstance(obj, BoundaryCondition):
            return obj.p
        return super(CustomEncoder, self).default(obj)

class CustomDecoder(json.JSONDecoder):
    """This custom JSON decoder can handle complex and Fraction and BoundaryCondition types.
    """

    def __init__(self, *args, **kwargs):
        super(CustomDecoder, self).__init__(*args, object_hook=self.object_hook, **kwargs)

    restore = {
        'complex': lambda d: complex(d["real"], d["imag"]),
        'Fraction': lambda d: Fraction(d["numerator"], d["denominator"]),
    }

    @staticmethod
    def object_hook(d):
        if "__class__" not in d:
            return d
        return CustomDecoder.restore[d["__class__"]](d)

dumps = partial(json.dumps, cls=CustomEncoder)
dump = partial(json.dump, cls=CustomEncoder)

loads = partial(json.loads, cls=CustomDecoder)
load = partial(json.loads, cls=CustomDecoder)

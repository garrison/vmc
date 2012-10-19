import json
from functools import partial
from collections import OrderedDict
from fractions import Fraction

class CustomEncoder(json.JSONEncoder):
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
        return super(ComplexEncoder, self).default(obj)

dumps = partial(json.dumps, cls=CustomEncoder)
dump = partial(json.dump, cls=CustomEncoder)

loads = json.loads
load = json.load
import json
from functools import partial

class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, complex):
            return [obj.real, obj.imag]
        return super(ComplexEncoder, self).default(obj)

dumps = partial(json.dumps, cls=ComplexEncoder)
dump = partial(json.dump, cls=ComplexEncoder)

loads = json.loads
load = json.load

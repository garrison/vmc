import os
from errno import EEXIST

def mkdir_p(path):
    'Run the equivalent of "mkdir -p"'
    # http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != EEXIST or not os.path.isdir(path):
            raise

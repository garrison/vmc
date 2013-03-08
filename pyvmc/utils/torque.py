import re

def make_flag_property(pattern, post=lambda x: x, default=None):
    r = re.compile(pattern)
    def property_func(self):
        for line in self.lines:
            m = r.match(line)
            if m is not None:
                return post(m.group(1))
        return default
    return property(property_func)

def get_array_ids(s):
    try:
        # it might be a single number, in which case it represents an array
        return tuple(range(int(s)))
    except ValueError:
        rv = []
        for x in s.split(','):
            try:
                # it might be a array id
                rv.append(int(x))
            except ValueError:
                # it must be a range
                begin, end = x.split('-')
                rv.extend(range(int(begin), int(end) + 1))
        return tuple(rv)

class PBSFlags(object):
    def __init__(self, pbs_filename):
        with file(pbs_filename) as f:
            self.lines = [line[4:].strip() for line in f.readlines() if line.startswith('#PBS')]

    walltime = make_flag_property(r'^-l walltime=(.*)$')
    nodes = make_flag_property(r'^-l nodes=(\d+):ppn=\d+$', int, default=1)
    ppn = make_flag_property(r'^-l nodes=\d+:ppn=(\d+)$', int, default=1)
    array_ids = make_flag_property(r'^-t (.*)$', get_array_ids)

if __name__ == "__main__":
    flags = PBSFlags('/home/garrison/vmc-studies/electron-ring/2leg/rho-2_3/dmetal-rho-two-thirds-scan.py')
    print(flags.walltime)
    print(flags.nodes)
    print(flags.ppn)
    print(flags.array_ids)

    print get_array_ids("5")
    print get_array_ids("6")
    print get_array_ids("10-13,50, 70-80")
    print get_array_ids("20,40,60")

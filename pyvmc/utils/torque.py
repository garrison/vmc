import re
import os
import logging
from itertools import islice

logger = logging.getLogger(__name__)

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
        with open(pbs_filename) as f:
            self.lines = [line[4:].strip() for line in f.readlines() if line.startswith('#PBS')]

    walltime = make_flag_property(r'^-l\s+walltime=(.*)$')
    nodes = make_flag_property(r'^-l\s+nodes=(\d+):ppn=\d+$', int, default=1)
    ppn = make_flag_property(r'^-l\s+nodes=\d+:ppn=(\d+)$', int, default=1)
    array_ids = make_flag_property(r'^-t\s+(.*)$', get_array_ids)

def run_pbs_thread(do_work, iterable, src_filename):
    if "PBS_JOBID" not in os.environ:
        # not running under PBS, so we only have one thread
        do_work(iterable)
        return

    # fixme: some day use python threading instead
    # (http://docs.python.org/2/library/queue.html)

    array_id = int(os.environ.get("PBS_ARRAYID", 0))

    pbs_flags = PBSFlags(src_filename)
    assert pbs_flags.nodes == 1
    ppn = pbs_flags.ppn
    array_ids = pbs_flags.array_ids

    step = ppn * len(array_ids)

    # fork into processes numbered in the range [0, ppn)
    pid = 0
    for core_id in range(ppn - 1):
        pid = os.fork()
        if pid:
            break
    else:
        core_id = ppn - 1

    offset = array_ids.index(array_id) * ppn + core_id
    logger.info("offset %d of %d", offset, step)

    try:
        do_work(islice(iterable, offset, None, step))
    finally:
        if pid:
            os.waitpid(pid, 0)

def get_job_id():
    """Returns a (somewhat) unique identifier based on the current job or process id"""
    if "PBS_JOBID" in os.environ:
        return "job{}".format(re.match(r"^(\d+)", os.environ["PBS_JOBID"]).group(1))
    else:
        return "pid{}".format(os.getpid())

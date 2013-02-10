from functools import partial
import logging
import time
import resource

# RUSAGE_THREAD is new in python 3.2
getrusage = partial(resource.getrusage, getattr(resource, "RUSAGE_THREAD", resource.RUSAGE_SELF))

class log_rusage(object):
    """This is meant to be used as a context manager."""

    def __init__(self, logger, message=""):
        self.logger = logger
        self.message = message

    def __enter__(self):
        self.time = time.time()
        self.rusage = getrusage()

    def __exit__(self, exc_type, exc_value, traceback):
        rusage = getrusage()
        d = {
            'message': self.message,
            'walltime': time.time() - self.time,
            'utime': rusage.ru_utime - self.rusage.ru_utime,
            'stime': rusage.ru_stime - self.rusage.ru_stime,
        }
        logfunc = self.logger.info if (exc_type is None) else self.logger.error
        logfunc("%s", "{message} (utime {utime:.3f}; stime {stime:.3f}; walltime {walltime:.3f})".format(**d))

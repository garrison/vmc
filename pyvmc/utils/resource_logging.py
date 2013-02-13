from functools import partial, wraps
from sys import exc_info
import logging
import time
import resource

try:
    from contextlib import ContextDecorator
except ImportError:
    # ContextDecorator is new in python 3.2
    class ContextDecorator(object):
        def __call__(self, func):
            @wraps(func)
            def new_func(*args, **kwargs):
                self.__enter__()
                try:
                    rv = func(*args, **kwargs)
                except:
                    self.__exit__(*exc_info)
                    raise
                else:
                    self.__exit__(None, None, None)
                    return rv
            return new_func

# RUSAGE_THREAD is new in python 3.2
getrusage = partial(resource.getrusage, getattr(resource, "RUSAGE_THREAD", resource.RUSAGE_SELF))

class log_rusage(ContextDecorator):
    """This is meant to be used as a context manager.

    It also works as a decorator."""

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

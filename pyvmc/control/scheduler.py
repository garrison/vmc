import os
import logging

from twisted.internet import reactor, defer, protocol

logger = logging.getLogger(__name__)

class VmcProcessProtocol(protocol.ProcessProtocol):
    def __init__(self, d):
        self.d = d
        self.out = []
        self.err = []

    def outReceived(self, data):
        self.out.append(data)

    def errReceived(self, data):
        self.err.append(data)

    def processEnded(self, reason):
        output = ''.join(self.out)
        err = ''.join(self.err)
        logger.debug('OUTPUT: %s', output)
        logger.debug('ERROR: %s', err)
        if reason.value.exitCode != 0:
            logger.error("nonzero exit value: %s\n%s", output, err)
            raise Exception("process failed with exit code: %s" % reason.value.exitCode)
        self.d.callback((output, err))

class Scheduler(object):
    _cores_in_use = 0
    _n_cores = 1

    def __init__(self):
        self.queue = defer.DeferredQueue()

    def _run_command(self, d, path, input_str):
        logger.debug("INPUT: %s", input_str)
        logger.info("Beginning task; %d left in queue", len(self.queue.pending))
        process_protocol = VmcProcessProtocol(d)
        transport = reactor.spawnProcess(process_protocol, path, env=os.environ)
        transport.write(input_str)
        transport.closeStdin()

    def submit(self, path, input_str):
        d = defer.Deferred()
        self.queue.put((d, path, input_str))
        return d

    def _command_complete(self, args):
        self._cores_in_use -= 1
        assert self._cores_in_use >= 0
        self._start_processes()
        return args

    def _start_run(self, args):
        d, path, input_str = args
        self._run_command(*args)
        d.addCallback(self._command_complete)

    def _start_processes(self):
        while self._cores_in_use < self._n_cores:
            self._cores_in_use += 1
            self.queue.get().addCallback(self._start_run)

    def use_cores(self, n_cores=None):
        if n_cores is None:
            # try to guess
            import os
            from multiprocessing import cpu_count
            try:
                n_cores = int(os.environ["PBS_NUM_PPN"])
            except KeyError:
                n_cores = cpu_count()
        assert n_cores > 0
        self._n_cores = n_cores

    def run(self, n_cores=None):
        if n_cores is not None:
            self.use_cores(n_cores)
        self._start_processes()
        reactor.run()

default_scheduler = Scheduler()

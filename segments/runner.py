import multiprocessing
import time

from .fragment_set import FragmentSet
from .processor import FragmentSetProcessor


class Runner(object):

    def __init__(self, exp, keys, workers = 8):
        self._workers = workers
        self.experiment = exp
        self.keys = keys

    def _run(self, key):
        proc = FragmentSetProcessor(self.experiment, FragmentSet(key))
        proc.barcodePrune = False
        proc.autosave()
        proc.run()
        proc.saveResults()
        proc.classifier.saveReads(key)

    def _workerDispatch(self):
        from .util import _speedup
        assert(_speedup)
        workers = []
        for key in self.keys:
            worker = multiprocessing.Process(target = self._run, args = (key,))
            workers.append(worker)
            worker.start()
            while len(workers) >= self._workers:
                time.sleep(1)
                workers = [ w for w in workers if w.is_alive() ]
        for worker in workers:
            worker.join()

    def runAll(self):
        self._workerDispatch()

    def runSingle(self, diagram = False, allDiagrams = False):
        for key in self.keys:
            proc = FragmentSetProcessor(self.experiment, FragmentSet(key))
            proc._diagram = bool(diagram)
            proc._allDiagrams = bool(allDiagrams)
            proc.barcodePrune = False
            proc.autosave()
            proc.run()
            proc.saveResults()
            proc.classifier.saveReads(key)

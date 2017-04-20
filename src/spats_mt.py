import thread
import Queue
import multiprocessing

class Spats(object):

    def __init__(self, target_path, data_r1_path, data_r2_path, output_folder):
        self.target_path = target_path
        self.data_r1_path = data_r1_path
        self.data_r2_path = data_r2_path
        self.output_folder = output_folder
        self.num_workers = mulitprocessing.cpu_count()

    def setup(self):
        if self.workers:
            return
        self.tasks = Queue.Queue()
        self.workers = ...
        # TODO: have a generic set of workers prepared which can handle tasks
        # (such as indexing task, or counting task)

    def indexTargets(self):
        if self.targets:
            return
        # TODO: use generic workers as above
        targets = Queue.Queue()
        threads = []
        def worker(name, seq):
            target = Target(name, seq)
            target.index()
            targets.put(target)

        for name, seq in zip(fasta_parse(self.target_path)):
            thd = threading.Thread(target = worker, args = (name, seq))
            threads.append(thd)
            thd.start()
        for thd in threads:
            thd.join()
        assert(len(targets) == len(threads))
        self.targets = list(targets)

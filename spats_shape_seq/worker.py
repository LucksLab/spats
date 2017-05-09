
import multiprocessing
import Queue
import sys

from pair import Pair
from util import _debug, _warn


class SpatsWorker(object):
    '''Manages multiprocessing aspects of Spats.
    '''

    def __init__(self, run, processor, pair_db, result_set_id = None):
        self._run = run
        self._processor = processor
        self._pair_db = pair_db
        self._result_set_id = result_set_id
        self._workers = []

    def _worker(self, worker_id):
        try:
            processor = self._processor
            self._pair_db.worker_id = worker_id
            writeback = bool(self._result_set_id)
            pair = Pair()
            batches = 0
            while True:
                pairs = self._pairs_to_do.get()
                if not pairs:
                    break
                results = []
                for lines in pairs:
                    pair.set_from_data('', lines[1], lines[2], lines[0])
                    processor.process_pair(pair)
                    if writeback:
                        results.append((lines[3],
                                        pair.target.rowid if pair.target else None,
                                        pair.site if pair.has_site else -1,
                                        pair.mask.chars if pair.mask else None,
                                        pair.multiplicity,
                                        pair.failure))
                                        
                if writeback:
                    self._results.put(results)

                batches += 1
                if not self._run.quiet:
                    sys.stdout.write('.')
                    sys.stdout.flush()

            self._pairs_done.put((processor.counters._counts, [(m.total, m.kept, m.count_data()) for m in processor._masks]))
        except:
            print "**** Worker exception, aborting..."
            raise

    def _createWorkers(self, num_workers):
        if 1 == num_workers:
            return
        for i in range(num_workers):
            worker = multiprocessing.Process(target = self._worker, args = (i,))
            self._workers.append(worker)
            worker.start()
        if not self._run.quiet:
            print "Created {} workers".format(num_workers)

    def _joinWorkers(self):

        if not self._workers:
            self._worker(0)
            return

        for w in self._workers:
            self._pairs_to_do.put(None) # just a dummy object to signal we're done
        for w in self._workers:
            w.join()

    def run(self, db_iterator):

        num_workers = max(1, self._run.num_workers or multiprocessing.cpu_count())
        self._pairs_to_do = multiprocessing.Queue(maxsize = 2 * num_workers)
        self._pairs_done = multiprocessing.Queue()
        self._results = multiprocessing.Queue()

        self._createWorkers(num_workers)

        quiet = self._run.quiet
        more_pairs = True
        pair_db = self._pair_db
        writeback = bool(self._result_set_id)
        num_batches = 0
        if writeback:
            result_set_id = self._result_set_id

        def put_batch():
            pair_info = next(db_iterator)
            self._pairs_to_do.put(pair_info)
            if not quiet:
                sys.stdout.write('^')
                sys.stdout.flush()

        def get_results():
            results = self._results.get_nowait()
            pair_db.add_results(self._result_set_id, results)
            if not quiet:
                sys.stdout.write('v')
                sys.stdout.flush()

        while more_pairs:
            try:
                cur_count = 0
                while cur_count < num_workers or num_batches < 2 * num_workers:
                    put_batch()
                    num_batches += 1
                    cur_count += 1
                if writeback:
                    while True:
                        get_results()
                        num_batches -= 1
            except StopIteration:
                more_pairs = False
            except Queue.Empty:
                pass

        if writeback:
            while num_batches > 0:
                results = self._results.get()
                pair_db.add_results(self._result_set_id, results)
                num_batches -= 1
                if not quiet:
                    sys.stdout.write('v')
                    sys.stdout.flush()

        self._joinWorkers()

        if not self._run.quiet:
            print "\nAggregating data..."

        processor = self._processor
        try:
            targets = { t.name : t for t in processor._targets.targets }
            while 1 < num_workers:
                data = self._pairs_done.get_nowait()
                their_counters = data[0]
                our_counters = processor.counters._counts
                for key, value in their_counters.iteritems():
                    if key != "_counts":
                        our_counters[key] = our_counters.get(key, 0) + value
                for i in range(len(processor._masks)):
                    m = processor._masks[i]
                    d = data[1][i]
                    m.total += d[0]
                    m.kept += d[1]
                    m.update_with_count_data(d[2], targets)
        except Queue.Empty:
            pass

        processor.counters.total_pairs = self._pair_db.count()
        processor.counters.unique_pairs = self._pair_db.unique_pairs()


import multiprocessing
import Queue
import sys

from pair import Pair
from parse import SamWriter
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

    def _make_result(self, ident, pair, tagged = False):
        res = [ ident,
                pair.target.rowid if pair.target else None,
                pair.mask.chars if pair.mask else None,
                pair.site if pair.has_site else -1,
                pair.end if pair.has_site else -1,
                pair.mutations[0] if pair.mutations else -1,
                pair.multiplicity,
                pair.failure ]
        if tagged:
            res.append(pair.tags)
        return res

    def _worker(self, worker_id):
        try:
            processor = self._processor
            if self._pair_db:
                self._pair_db.worker_id = worker_id
            writeback = bool(self._result_set_id)
            tagged = processor.uses_tags
            use_quality = self._run._parse_quality
            pair = Pair()
            while True:
                pairs = self._pairs_to_do.get()
                if not pairs:
                    break
                results = []
                for lines in pairs:
                    pair.set_from_data(lines[3], str(lines[1]), str(lines[2]), lines[0])
                    if use_quality:
                        pair.r1.quality = str(lines[4])
                        pair.r2.quality = str(lines[5])
                    processor.process_pair(pair)
                    #if pair.failure:
                    #    print('FAIL: {}'.format(pair.failure))
                    if writeback:
                        results.append(self._make_result(lines[3], pair, tagged))

                if writeback:
                    self._results.put(results)

                if not self._run.quiet:
                    sys.stdout.write('.')#str(worker_id))
                    sys.stdout.flush()

            self._pairs_done.put(processor.counters.count_data())
        except:
            print("**** Worker exception, aborting...")
            raise

    def _createWorkers(self, num_workers):
        for i in range(num_workers):
            worker = multiprocessing.Process(target = self._worker, args = (i,))
            self._workers.append(worker)
            worker.start()
        if not self._run.quiet:
            print("Created {} workers".format(num_workers))

    def _joinWorkers(self):
        for w in self._workers:
            w.join()

    def run(self, pair_iterator):

        num_workers = max(1, self._run.num_workers or multiprocessing.cpu_count())

        if 1 == num_workers:
            self.run_simple(pair_iterator)
            return

        self._pairs_to_do = multiprocessing.Queue(maxsize = 2 * num_workers)
        self._pairs_done = multiprocessing.Queue()
        self._results = multiprocessing.Queue()

        self._createWorkers(num_workers)

        quiet = self._run.quiet
        more_pairs = True
        pair_db = self._pair_db
        writeback = bool(self._result_set_id)
        num_batches = 0
        total = 0
        if writeback:
            result_set_id = self._result_set_id

        def put_batch():
            pair_info = next(pair_iterator)
            self._pairs_to_do.put(pair_info)
            if not quiet:
                sys.stdout.write('^')
                sys.stdout.flush()
            return len(pair_info)

        def write_results():
            all_results = []
            num_batches = 0
            try:
                while True:
                    all_results.extend(self._results.get(True, 0.01))
                    num_batches += 1
                    if not quiet:
                        sys.stdout.write('v')
                        sys.stdout.flush()
            except Queue.Empty:
                pass
            if all_results:
                pair_db.add_results(self._result_set_id, all_results)
            return num_batches

        while more_pairs:
            try:
                cur_count = 0
                while cur_count < num_workers or num_batches < 2 * num_workers:
                    total += put_batch()
                    num_batches += 1
                    cur_count += 1
                if writeback:
                    num_batches -= write_results()
            except StopIteration:
                more_pairs = False
            except Queue.Empty:
                pass

        if writeback:
            while num_batches > 0:
                num_batches -= write_results()

        # put signal objects to indicate we're done
        for i in range(num_workers):
            self._pairs_to_do.put(None)

        processor = self._processor
        targets = { t.name : t for t in processor._targets.targets }
        accumulated = 0

        def accumulate_counts():
            num_accumulated = 0
            try:
                while 1 < num_workers:
                    data = self._pairs_done.get_nowait()
                    processor.counters.update_with_count_data(data)
                    num_accumulated += 1
                    if not quiet:
                        sys.stdout.write('x')
                        sys.stdout.flush()
            except Queue.Empty:
                pass
            return num_accumulated

        while accumulated < num_workers:
            accumulated += accumulate_counts()

        if not self._run.quiet:
            print("\nAggregating data...")

        self._joinWorkers()

        processor.counters.total_pairs = total
        #if self._pair_db:
        #    processor.counters.unique_pairs = self._pair_db.unique_pairs()


    def run_simple(self, pair_iterator):

        quiet = self._run.quiet
        run_limit = self._run._run_limit
        more_pairs = True
        pair_db = self._pair_db
        writeback = bool(self._result_set_id)
        sam = bool(self._run.generate_sam)
        use_quality = self._run._parse_quality
        total = 0
        if writeback:
            result_set_id = self._result_set_id

        processor = self._processor
        if self._pair_db:
            self._pair_db.worker_id = 0
        tagged = processor.uses_tags
        pair = Pair()

        if sam:
            sam_writer = SamWriter(self._run.generate_sam, processor._targets.targets)

        while more_pairs:
            try:
                while True:
                    pair_info = next(pair_iterator)
                    if not quiet:
                        sys.stdout.write('^')
                        sys.stdout.flush()
                    results = []
                    for lines in pair_info:
                        pair.set_from_data(lines[3], str(lines[1]), str(lines[2]), lines[0])
                        if use_quality:
                            pair.r1.quality = str(lines[4])
                            pair.r2.quality = str(lines[5])

                        try:
                            processor.process_pair(pair)
                        except:
                            print("**** Error processing pair: {} / {}".format(pair.r1.original_seq, pair.r2.original_seq))
                            raise
                        if sam:
                            sam_writer.write(pair)

                        total += pair.multiplicity
                        if writeback:
                            results.append(self._make_result(lines[3], pair, tagged))

                    if not quiet:
                        sys.stdout.write('v')
                        sys.stdout.flush()

                    if results:
                        pair_db.add_results(self._result_set_id, results)
                        if not quiet:
                            sys.stdout.write('.')
                            sys.stdout.flush()

                    if run_limit and total > run_limit:
                        raise StopIteration()

            except StopIteration:
                more_pairs = False

        if not self._run.quiet:
            print("\nAggregating data...")

        processor.counters.total_pairs = total
        if self._pair_db:
            processor.counters.unique_pairs = self._pair_db.unique_pairs()

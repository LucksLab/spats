"""For reads analysis, first create a :class:`.ReadsData` instance,
passing it the path to use for its database (.db) file:

.. code-block:: python

    from spats_shape_seq.reads import ReadsData, ReadsAnalyzer
    data = ReadsData('~/path/to/reads.db')

Then, parse the reads data from the original input files: pass the
target file you use for analysis (single target only -- without linker
-- for cotrans experiments), and the R1/R2 input paths:

.. code-block:: python

    data.parse('~/path/to/target.fa', '~/path/to/data/R1.fastq', '~/path/to/data/R2.fastq')

Note that, by default, this only parses a sample of 100k of the R1/R2
pairs. This should be sufficient for analyzing the reads, and is
required for reasonable performance of the interactive UI.

Once you have the data parsed, it needs to be analyzed (tagged); this
is similar to a normal :class:`.spats.Spats` run, but more information
is gathered for non-matching pairs:

.. code-block:: python

    analyzer = ReadsAnalyzer(data)

Configure the run, if desired, based on the :class:`.run.Run`; for
example, to change the linker from the default for a cotrans
experiment:

.. code-block:: python

    analyzer.run.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'

Then, perform the analysis:

.. code-block:: python

    analyzer.process_tags()

Once that completes, your reads data will be ready for analysis.

"""

import os
import time

from spats_shape_seq import Spats
from spats_shape_seq.db import PairDB
from spats_shape_seq.util import reverse_complement


class ReadsData(object):
    """Encapsulates the data for reads analysis.

       :param db_path: the path for the reads data
    """

    def __init__(self, db_path):
        self.db_path = db_path
        self._pair_db = None

    def parse(self, target_path, r1_paths, r2_paths, sample_size = 100000, show_progress_every = 1000000):
        """Used to parse target and r1/r2 reads data for reads processing.

           :param target_path: Path to targets file, must be in FASTA format.

           :param r1_paths: List of paths to R1 reads files, must be in FASTQ format.

           :param r2_paths: List of paths to R2 reads files, must be in FASTQ format.

           :param sample_size: the number of samples to use for analysis. Samples will be (more or less) uniformly
                               randomly selected from the population of pairs.

           :param show_progress_every: default show a '.' for every 1 million pairs parsed. Set to 0 to disable output.
        """
        start = time.time()
        print("Parsing to db...")
        pair_db = self.pair_db
        if show_progress_every:
            pair_db.show_progress_every = show_progress_every
        pair_db.wipe()
        pair_db.add_targets_table(target_path)
        total = 0
        for r1_path, r2_path in zip(r1_paths, r2_paths):
            total += pair_db.parse(r1_path, r2_path, sample_size = sample_size)
        print("Sampled {} records in {:.1f}s".format(total, time.time() - start))

    @property
    def pair_db(self):
        """Access the underlying :class:`.db.PairDB`.
        """
        if not self._pair_db:
            self._pair_db = PairDB(self.db_path)
        return self._pair_db




class ReadsAnalyzer(object):
    """Performs the analysis/tagging required for the reads analysis.

       :param reads_data: the :class:`.ReadsData` with the input data.
       :param cotrans: pass `True` for cotrans-style experiments.
    """

    def __init__(self, reads_data, cotrans = False):
        self._reads_data = reads_data
        self._pair_db = reads_data.pair_db
        s = Spats()
        s.run._p_use_tag_processor = True
        s.run.cotrans = cotrans
        s.run.writeback_results = True
        s.run.result_set_name = "tags"
        s.run.allow_indeterminate = True
        s.run.allowed_target_errors = 2
        s.run.allowed_adapter_errors = 2
        s.run.num_workers = 8
        self._spats = s
        self._extra_tag_targets = []
        self._plugins = {}

    @property
    def run(self):
        """Provides access to the :class:`.run.Run` which is used to configure tag analysis.
        """
        return self._spats.run

    def addTagTarget(self, name, tag):
        self._extra_tag_targets.append((name, tag))
        self.run._p_extra_tags = self._extra_tag_targets

    def addTagPlugin(self, tag, handler):
        self._plugins[tag] = handler

    def process_tags(self):
        """Processes the tags in the input data for analysis.
        """
        s = self._spats
        pair_db = self._pair_db
        s.loadTargets(pair_db)
        s.run.pair_length = pair_db.pair_length()
        p = s._processor
        for target in pair_db.targets():
            p.addTagTarget(target[0], target[1])
            p.addTagTarget(target[0] + "_rc", reverse_complement(str(target[1])))
        p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
        p.addTagTarget("adapter_b", s.run.adapter_b)
        if s.run.cotrans:
            p.addTagTarget("linker_cotrans", s.run.cotrans_linker)
            p.addTagTarget("linker_cotrans_rc", reverse_complement(s.run.cotrans_linker))
        for tag in self._extra_tag_targets:
            p.addTagTarget(tag[0], tag[1])
        for tag, handler in self._plugins.iteritems():
            p.addTagPlugin(tag, handler)

        s.process_pair_db(pair_db, batch_size = 10000)
        self.result_set_id = pair_db.result_set_id_for_name(s.run.result_set_name)
        pair_db.count_tags(self.result_set_id)
        pair_db.store_run(s.run)

    def tag_counts(self):
        """Returns a dictionary of ``{ tag: count }`` of tags analyzed.
        """
        return self._pair_db.tag_counts(self.result_set_id)

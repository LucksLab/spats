

spats_tool
==========

The ``spats_tool`` command can be used to perform most of the standard
SPATS tasks. The first thing to do to use the ``spats_tool`` is to
designate a directory for your experiment. In general, this should
probably be a new directory for each new data set; for organizational
purposes, it's best not to have too many other files in the
folder. The directory needs to be on the same filesystem as your data
files.

As an example, we'll create a folder for a new "testing" experiment:

.. code-block:: bash

    $ mkdir testing_jjb_2017_08_10
    $ cd testing_jjb_2017_08_10

The tool works using a ``spats.config`` file, which is a text file
that provides the information required for the various tool
functions. The primary information is about the locations of the
target and input files. Here's a sample:

.. code-block:: text

    [spats]
    target = /path/to/data/my_target.fa
    r1 = /path/to/data/experiment_data_R1.fastq.gz
    r2 = /path/to/data/experiment_data_R2.fastq.gz

The ``target``, ``r1``, and ``r2`` parameters are
required for most tools. If this is a cotrans experiment,
set ``cotrans = True``; then provide the paths to the
target/R1/R2 files. Note that you can use ``.gz`` files for r1/r2, in
which case they will be decompressed on the fly, and cleaned up after
the run completes.

The ``spats_tool`` command must be run in the experiment directory you
created with the ``spats.config`` file.

spats_tool run
--------------

The basic command is ``run``, which performs the SPATS run to compute
site reactivities:

.. code-block:: bash

    $ spats_tool run
    :using native cotrans processor
    :decompressing /projects/b1044/.../EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq.gz
    :decompress R1 @ 41.51s
    :decompressing /projects/b1044/.../EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq.gz
    :decompress R1 @ 93.62s
    :wrote output to run.spats
    :run complete @ 134.99s

All ``spats_tool`` work in the experiment directory and update the
``spats.log`` file there; for example, in this case, it looks like:

.. code-block:: text

 2017/08/10 13:22 : run, 134.99s
   - ** removing previous run.spats
   - using native cotrans processor
   - decompressing /projects/b1044/.../EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq.gz
   - decompress R1 @ 47.76s
   - decompressing /projects/b1044/.../EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq.gz
   - decompress R1 @ 93.00s
   - wrote output to run.spats
   - run complete @ 134.99s

As the output and log indicate, the results of the run are written to
the ``run.spats`` file, which is a sqlite-DB file that can be used by
``spats_tool`` and the visualization tool. All tools append to this
log file, so you have a record of all analyses performed, including
date/time stamps. To do matrix visualization, copy the ``run.spats``
file to your Mac laptop and open it with the visualization tool.


spats_tool reads
----------------

The ``reads`` command analyzes the experimental data and creates a
``reads.spats`` file, which can be used with the visualization tool to
analyze the quality of the data.

.. code-block:: bash

    $ spats_tool reads
    :** removing previous reads.spats
    :using native reads
    Lookup table: 1076 R1 entries, 121 R2 entries.
    Lookup table: 1076 R1 entries, 121 R2 entries.
    Processing pairs...
    Created 8 workers
    ^^^^^^^^^.v........vvvvvvvvxxxxxxxx
    Aggregating data...
    Successfully processed 3640 properly paired fragments:
      ...
    :tags processed to reads.spats
    :reads complete @ 59.38s

To do reads visualization, copy the ``reads.spats`` file to your Mac
laptop and open it with the visualization tool.


spats_tool validate
-------------------

The ``validate`` command is used after the ``run`` command: it re-runs
the SPATS analysis on the input data, using a different (slower)
algorithm, and then verifies that the results match.

.. code-block:: bash

    $ spats_tool validate
    Processing pairs...
    Created 20 workers
    ^^^...
    Aggregating data...
    Successfully processed 2257112 properly paired fragments:
      ...
    Total time: (195.5s)
    Original results (native algorithm) validated using find_partial algorithm, 17402 registered sites match.
    :Validation pass
    :validate complete @ 195.68s

Any mismatches will result in an error, and should be reported as a
bug!


spats_tool dump
-------------------

The ``dump`` command is used to access the raw data and dump it to
CSV. Requires a dump type -- options:

- ``spats_tool dump reads``: dumps the tags data for the reads
  analysis to ``reads.csv``

- ``spats_tool dump run``: dumps the treated/untreated count, beta,
  theta, and rho values from the run analysis to CSV files named for
  the corresponding targets.


..
   Miscellaneous
   -------------

   Other commands:

   - ``spats_tool viz``: used to launch the visualization tool UI (Mac only)

   - ``spats_tool help``: used to show usage help



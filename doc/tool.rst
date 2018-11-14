

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
    $ spats_tool init

The tool works using a ``spats.config`` file, which is a text file
that provides the information required for the various tool
functions. ``spats_tool init`` will create a default config file, 
be sure to edit it after it's created. The rest of the information is
about the locations of the target and input files. Here's a sample:

.. code-block:: text

    [spats]
    preseq = my_preseq_file.fsa  # optional
    target = my_target.fa
    r1 = /path/to/data/experiment_data_R1.fastq.gz
    r2 = /path/to/data/experiment_data_R2.fastq.gz
    cotrans = False
    
    [metadata]
    ...

The ``preseq`` parameter is optional and only required for the ``pre`` tool.
The ``target``, ``r1``, and ``r2`` parameters are
required for most tools. If this is a cotrans experiment,
set ``cotrans = True``; then provide the paths to the
target/R1/R2 files. Note that you can use ``.gz`` files for r1/r2, in
which case they will be decompressed on the fly, and cleaned up after
the run completes.

If you wish to set any configuration for the SPATS run or reads
analysis (see below), you can set parameters according to the
:class:`.run.Run` documentation; for example,
``minimum_target_match_length = 12``.

The ``spats_tool`` command must be run in the experiment directory you
created with the ``spats.config`` file.


spats_tool pre
--------------

The first command you may wish to run is ``pre``, which takes the
information from a ``*.fsa`` (ABIF) file and creates a plot.

.. code-block:: bash

    $ spats_tool pre
    :pre-sequencing data processed to pre.spats
    :pre complete @ 0.02s

This extracts the data and stores it a format from which it can be
easily inspected and plotted.


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
``spats_tool`` to dump results and create plots. All tools append to this
log file, so you have a record of all analyses performed, including
date/time stamps.

.. To do matrix visualization, copy the ``run.spats``
.. file to your Mac laptop and open it with the visualization tool.


..
  spats_tool reads
  ----------------

  The ``reads`` command analyzes the experimental data and creates a
  ` `reads.spats`` file, which can be used with the visualization tool to
  analyze the quality of the data.

  . . code-block:: bash

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


spats_tool dump
-------------------

The ``dump`` command is used to access the raw data and dump it to
CSV. Requires a dump type -- options:

- ``spats_tool dump reads``: dumps the tags data for the reads
  analysis to ``reads.csv``

- ``spats_tool dump run``: dumps the treated/untreated count, beta,
  theta, and rho values from the run analysis to CSV files named for
  the corresponding targets.


spats_tool nb
-------------

To see plots and perform further analyses on any data, run
``spats_tool nb``. This brings up a Jupyter Notebook, which is
automatically updated after any analysis runs. For example, after
``spats_tool pre``, a plot of the presequencing data will be added to
the notebook.

The first cell in the notebook contains experiment metadata.

By default, the code require to create plots is hidden. To show the
code, click on the ``In [ ]:`` prompt above a given plot or
output. From there, you can edit the code to change plot styling or
parameters. Click the "Run" button at the top of the page to re-run a
cell after making changes.

Note that you can copy, edit, delete, and rearrange cells. To create a
new cell, hit the "+" button. By default, cells are code cells; they
can also be markdown, to change this use the dropdown menu to the
right of the "Run" button. For example, you may wish to create a new
markdown cell to write down notes about the experiment.

You may wish to edit the look and feel of the plots. Plots are created
using ``matplotlib``. (Support for plots in ``R`` is coming soon.) You
may wish to consult the `tutorial
<https://matplotlib.org/users/pyplot_tutorial.html>`_, the `cheat sheet
<https://s3.amazonaws.com/assets.datacamp.com/blog_assets/Python_Matplotlib_Cheat_Sheet.pdf>`_,
or the `documentation <https://matplotlib.org/contents.html>`_.

..
  spats_tool validate
  -------------------

  The ``validate`` command is used after the ``run`` command: it re-runs
  the SPATS analysis on the input data, using a different (slower)
  algorithm, and then verifies that the results match.

  . . code-block:: bash

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


..
   Miscellaneous
   -------------

   Other commands:

   - ``spats_tool viz``: used to launch the visualization tool UI (Mac only)

   - ``spats_tool help``: used to show usage help



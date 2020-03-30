
Quickstart Guide
================

This gives a quick overview of running spats. For more detail, see the
`spats_tool <tool.html>`_ documentation or the `detailed options
<reference.html>`_.

First, create a folder for your experiment, and run ``spats_tool init``:

.. code-block:: bash

    $ mkdir testing_2017_08_10
    $ cd testing_2017_08_10
    $ spats_tool init

This creates a default ``spats.config`` file. The minimum required
options for this file are as follows:

.. code-block:: text

    [spats]
    target = my_target.fa
    r1 = /path/to/data/experiment_data_R1.fastq.gz
    r2 = /path/to/data/experiment_data_R2.fastq.gz

See the :class:`.run.Run` documentation for detailed information about
all available options.

You may now wish to do a reads analysis:

.. code-block:: bash

  $ spats_tool reads

Use ``spats_tool dump reads`` to dump the data for the reads analysis
to ``reads.csv``.

The basic command is ``spats_tool run``, which performs the SPATS run
to compute site reactivities:

.. code-block:: bash

    $ spats_tool run

Use ``spats_tool dump run`` to dump the data for the run to CSV files.

For more details, see the `spats_tool <tool.html>`_ documentation or
the `detailed options <reference.html>`_.

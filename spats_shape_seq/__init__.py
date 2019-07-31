"""Typically, you'll create a :class:`.spats.Spats` instance:

.. code-block:: python

    from spats_shape_seq import Spats
    spats = Spats()

Then, configure the run, if desired, based on the :class:`.run.Run`:

.. code-block:: python

    spats.run.minimum_target_length = 10
    ...

Add the .fa file with the target(s):

.. code-block:: python

    spats.addTargets(path_to_targets_fasta)

That sets up the Spats processor. To process input data, call
:meth:`.spats.Spats.process_pair_data`, passing the path to the R1 and
R2 input files:

.. code-block:: python

    spats.process_pair_data(path_to_R1_fastq, path_to_R2_fastq)

Finally, call :meth:`.spats.Spats.compute_profiles` to compute
beta/theta values, and :meth:`.spats.Spats.write_reactivities` to
output the results to a file:

.. code-block:: python

    spats.compute_profiles()
    spats.write_reactivities(path_to_reactivities_out)

For a case that requires no non-default configuration or processing,
you can use the :func:`run_spats` convenience function, which will
perform all of the above steps.

For saving data for later analysis / output / visualization, use the
:meth:`.spats.Spats.store` method:

.. code-block:: python

    spats.store('my_run.spats')

This can be later loaded using the :meth:`.spats.Spats.load` method.


"""


from spats import Spats


_VERSION = "1.9.34"
_PRODUCTION = True
_PUBLIC_RELEASE = True
version = _VERSION + (" Public" if _PUBLIC_RELEASE else " Private") + (" Production" if _PRODUCTION else " Beta")


def run_spats(target_path, r1_path, r2_path, output_path, cotrans = False):
    """Convenience function for a common-case SPATS run that doesn't
       require any non-default configuration.

    :param target_path: path to the targets FASTA file
    :param r1_path: path to the R1 input data FASTQ file
    :param r2_path: path to the R2 input data FASTQ file
    :param output_path: path to write resulting reactivities
    :param cotrans: pass `True` for cotrans-style experiments.

    """

    spats = Spats(cotrans = cotrans)
    spats.addTargets(target_path)
    spats.process_pair_data(r1_path, r2_path)
    spats.compute_profiles()
    spats.write_reactivities(output_path)



# jupyter extension support

def _jupyter_server_extension_paths():
    return [{
        "module": "spats_shape_seq"
    }]

def _jupyter_nbextension_paths():
    return [ { "section" : "notebook",
               "src" : "static",
               "dest" : "spats_shape_seq",
               "require" : "spats_shape_seq/main" } ]

def load_jupyter_server_extension(nbapp):
    pass

"""Typically, you'll create a :class:`.spats.Spats` instance:

.. code-block:: python

    from spats_shape_seq import Spats
    spats = Spats()

Then, configure, if desired, based on the :class:`.config.SpatsConfig`:

.. code-block:: python

    spats.config.minimum_target_length = 10

Add the masks used for treated and untreated samples (in that order):

.. code-block:: python

    spats.addMasks('RRRY', 'YYYR')

And the .fa file with the target(s):

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

"""


from config import spats_config
from spats import Spats
from util import reverse_complement


def run_spats(target_path, r1_path, r2_path, output_path, masks = ['RRRY', 'YYYR' ]):
    """Convenience function for a common-case SPATS run that doesn't
       require any non-default configuration. Uses `10` for
       :attr:`config.SpatsConfig.minimum_target_match_length`.

    :param target_path: path to the targets FASTA file
    :param r1_path: path to the R1 input data FASTQ file
    :param r2_path: path to the R2 input data FASTQ file
    :param output_path: path to write resulting reactivities
    :param masks: list of strings (two expected) to use as masks for treated and untreated handles.

    """

    spats = Spats()
    spats.config.minimum_target_match_length = 10
    spats.addMasks(*masks)
    spats.addTargets(target_path)
    spats.process_pair_data(r1_path, r2_path)
    spats.compute_profiles()
    spats.write_reactivities(output_path)

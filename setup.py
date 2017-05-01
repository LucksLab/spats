from setuptools import setup


DESCRIPTION = 'The Spats package implements a read mapping and reactivity analysis pipeline for calculating SHAPE-Seq reactivities from an input set of next-generation reads.'

README = '''
Spats processes reads and calculates SHAPE reactivities for SHAPE-Seq experiments on multiple RNAs. It accepts raw paired-end sequencing reads in fastq format, and a target sequence file containing the sequences of RNAs present in the experimental pool. Spats then performs read alignment to calculate distributions of read ends in the SHAPE (+) and (-) channel for each nucleotide in each RNA. Spats then estimates nucleotide resolution SHAPE reactivities for each RNA, using a model-driven maximum likelihood procedure based on a model of the reverse transcriptase process used in the SHAPE-Seq experiment.	Spats is a collaborative effort between the Aviran Lab at UC Davis, the Pachter Lab at UC Berkeley, the Trapnell Lab at the University of Washington, and the Lucks Lab at Northwestern University. 

Spats is provided under the OSI-approved Boost License.

http://luckslab.github.io/spats/
'''

setup(name='spats_shape_seq',
      version='1.9.0',
      description=DESCRIPTION,
      long_description=README,
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: Boost Software License 1.0 (BSL-1.0)',
          'Programming Language :: Python :: 2.7',
      ],
      keywords='spats',
      url='http://luckslab.github.io/spats',
      author='LucksLab',
      author_email='spats.shape@gmail.com',
      license='MIT',
      packages=['spats_shape_seq'],
      install_requires=[],
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True,
      zip_safe=False)

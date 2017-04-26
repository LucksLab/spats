from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='spats_shape',
      version='1.9',
      description='The Spats package implements a read mapping and reactivity analysis pipeline for calculating SHAPE-Seq reactivities from an input set of next-generation reads.',
      long_description=readme(),
      classifiers=[
          'Development Status :: 4 - Beta',
          'Boost Software License 1.0 (BSL-1.0)',
          'Programming Language :: Python :: 2.7',
      ],
      keywords='spats',
      url='http://luckslab.github.io/spats',
      author='LucksLab',
      author_email='spats.shape@gmail.com',
      license='MIT',
      packages=['spats_shape'],
      install_requires=[],
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True,
      zip_safe=False)


Advanced Usage
==============

- Creating on-disk DB (using PairDB)

  - ``$ spats_tool makedb [db_path] [targets_path] [r1_path] [r2_path]``

  - ``$ spats_tool dbrun [db_path] [run_name]``

  - for comparing v102:

    - ``$ spats_tool addv102 [db_path] [targets_path] [v102_spats_out_path]``

    - ``$ spats_tool rdiff [db_path] [first_run_name] v102``

  - note ``spats_tool rdiff`` can be used for comparing runs with different branches/config/etc

  - TODO: change this from cmdl to pythonic

- Stats on input data

- Analyzing specific pairs

- Diagrams


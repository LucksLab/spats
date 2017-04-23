import pstats
p = pstats.Stats('tmp/runprof.out')
p.strip_dirs().sort_stats('cumulative').print_stats(100)

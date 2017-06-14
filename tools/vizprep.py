
import os
import sys

CFG_TEMPLATE = '''
[server]
host = 0.0.0.0
port = 7865
portfile = /tmp/uif.port

[data]
dbfile = {}
result_set_name = tags
'''

def run():

    if len(sys.argv) != 5:
        print "Usage: PYTHONPATH=. python tools/vizprep.py [targets_path] [data_r1_path] [data_r2_path] [db_out_path]"
        exit(1)

    target_path = sys.argv[1]
    r1_path = sys.argv[2]
    r2_path = sys.argv[3]
    db_path = sys.argv[4]

    from spats_shape_seq import Spats
    from spats_shape_seq.db import PairDB
    from spats_shape_seq.tag import TagProcessor
    from spats_shape_seq.util import reverse_complement

    pair_db = PairDB(db_path)
    print "Parsing to db..."
    pair_db.wipe()
    pair_db.add_targets_table(target_path)
    pair_db.parse(r1_path, r2_path, sample_size = 0)

    s = Spats()
    s.run._processor_class = TagProcessor
    s.run.writeback_results = True
    s.run.result_set_name = "tags"
    #s.run.num_workers = 1
    s.loadTargets(pair_db)

    s.run.allow_indeterminate = True
    s.run.allowed_target_errors = 2
    s.run.allowed_adapter_errors = 2

    p = s._processor
    for target in pair_db.targets():
        p.addTagTarget(target[0], target[1])
        p.addTagTarget(target[0] + "_rc", reverse_complement(str(target[1])))
        p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
        p.addTagTarget("adapter_b", s.run.adapter_b)

    s.process_pair_db(pair_db)
    rsid = pair_db.result_set_id_for_name(s.run.result_set_name)
    pair_db.count_tags(rsid)
    counts = pair_db.tag_counts(rsid)
    for k in counts.keys():
        print "  {}: {}".format(k, counts[k])

    cfg_path = os.path.expanduser("~/.spats_viz")
    cfg_data = CFG_TEMPLATE.format(os.path.abspath(db_path))
    open(cfg_path, 'wb').write(cfg_data)

    print "Data processed, and default config updated (~/.spats_viz). Update config if desired, and then use viz tool: `make viz`"

if __name__ == '__main__':
    run()

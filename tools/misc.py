
def id_to_site():
    from spats_common import id_to_site
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/t7/"
    id_to_site(bp + "NOMASK_1.fq", bp + "RRRY.sam", bp + "YYYR.sam", 143)


def make_subset():
    from spats_common import make_subset
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    make_subset(bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                bp + "t11/5s_missing.ids",
                bp + "t11/x")

def clean_sites():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    spats(bp + "5s/5S.fa",
          bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
          bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
          bp + "t11")

def d5s_writeback_run():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"

    from spats_shape_seq.db import PairDB
    pair_db = PairDB(bp + "dev_out/pairs.db")
    pair_db.add_targets_table(bp + "5s/5S.fa")

    from spats_shape_seq import Spats
    s = Spats()
    s.addTargets(bp + "5s/5S.fa")
    s.writeback_results = True
    s.result_set_name = "pure_python"
    s.process_pair_db(pair_db)

def d5s_run():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"

    from spats_shape_seq.db import PairDB
    pair_db = PairDB(bp + "dev_out/pairs.db")
    if False:
        pair_db.add_targets_table(bp + "5s/5S.fa")
        pair_db.parse(bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                      bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq")

    from spats_shape_seq import Spats
    from spats_shape_seq.partial import PartialFindProcessor
    s = Spats()
    #s.run._processor_class = PartialFindProcessor
    s.run.writeback_results = True
    #s.run.resume_processing = True
    s.run.result_set_name = "lookup"
    s.addTargets(bp + "5s/5S.fa")
    s.process_pair_db(pair_db)

def ligation_run():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/Shape_Seq_ligation/"
    from spats_shape_seq import Spats
    s = Spats()
    #s.config.debug = True
    #s.run.minimum_target_match_length = 10
    #s.run.num_workers = 1
    #from spats_shape_seq.partial import PartialFindProcessor
    #s.run._processor_class = PartialFindProcessor
    s.run.skip_database = True
    s.addTargets(bp + "panel_RNAs_complete.fa")
    s.process_pair_data(bp + "data/KEW1_S1_L001_R1_001.fastq",
                        bp + "data/KEW1_S1_L001_R2_001.fastq")

def cotrans_run():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/"
    from spats_shape_seq import Spats
    s = Spats()
    #from spats_shape_seq.partial import PartialFindProcessor
    #s.run._processor_class = PartialFindProcessor
    s.run.skip_database = True
    #s.run.writeback_results = True
    #s.run.resume_processing = True
    #s.run.result_set_name = "lookup"
    s.addTargets(bp + "F_wt.fa")
    s.process_pair_data(bp + "data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq",
                        bp + "data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq")
    s.compute_profiles()
    s.write_reactivities(bp + "dev_out/rx.out")

def cotrans_test():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/"
    from spats_shape_seq import Spats
    s = Spats()
    from spats_shape_seq.partial import PartialFindProcessor
    #s.run._processor_class = PartialFindProcessor
    s.addTargets(bp + "F_wt.fa")
    from spats_shape_seq.pair import Pair
    pair = Pair()
    pair.set_from_data('x', 'GAGCGTCCTTGGTGCCCGAGTCAGAAATAGACTCCT', 'TATCACTACTGGTAGGAGTCTATTTCTGACTCGGGC')
    s.process_pair(pair)
    print "{}: {}".format(pair.target.name, pair.site)

def run_5sq():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
    spats(bp + "5S.fa",
          bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
          bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq", 
          bp + "t4",
          show_sites = False)

def misc():
    if False:
        bp = "/Users/jbrink/mos/tasks/1RwIBa/refactor/spats/test/Read_Mapping/"
        spats(bp + "SRP_All_Stops.fa", bp + "SRP_All_Stops_R1.fq", bp + "SRP_All_Stops_R2.fq", bp + "t5")
    elif False:
        bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
        spats(bp + "5s/5S.fa",
              bp + "5s_sample/filtered_R1.fq",
              bp + "5s_sample/filtered_R2.fq",
              bp + "t11")
    else:
        bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
        spats(bp + "5S.fa",
              bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
              bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq", 
              bp + "t4")

def spats(target, r1, r2, out, show_sites = True):
    from spats_shape_seq import Spats, spats_config
    s = Spats()
    s.addTargets(target)
    s.addMasks("RRRY", "YYYR")
    if show_sites:
        spats_config.show_id_to_site = True
    s.process_pair_data(r1, r2)
    if not show_sites:
        s.compute_profiles()
        s.write_reactivities(out + "/rx.out")
              
def test_refactor():
    from spats_clean import Spats
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
    out = bp + "t3/"
    s = Spats(bp + "5S.fa", out)
    s.setup()
    s.process_pair_data(bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
                        bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq")
    s.compute_profiles()
    s.write_reactivities()
    import subprocess
    subprocess.check_call(["diff", bp + "t2/rx.out", out + "/rx.out"])
    print "Diff OK"

def make_id_case(arg = None, show = True, skip_site = False):
    import subprocess
    ident = arg or sys.argv[2]
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    R1 = subprocess.check_output([ "awk", "/" + ident + "/{getline; print; exit}", bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq" ]).strip('\n')
    R2 = subprocess.check_output([ "awk", "/" + ident + "/{getline; print; exit}", bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq" ]).strip('\n')
    if skip_site:
        site = None
    else:
        site = subprocess.check_output([ "awk", "/" + ident + "/{print $3; exit}", bp + "t11/5s-counts_sorted.out"]).strip('\n')
    case = "    [ '*{}', '{}', '{}', {} ],".format(ident, R1, R2, site or None)
    if show:
        print case
    return [ ident, R1, R2, site or None ]

def show_failure_types():
    from spats_clean import Spats, Pair, FastqRecord
    spats = Spats("test/5s/5s.fa", "test/5s")
    spats.setup()
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"

    with open(bp + "t11/x/filtered_R1.fq", 'rb') as r1_in:
        with open(bp + "t11/x/filtered_R2.fq", 'rb') as r2_in:
            r1_record = FastqRecord()
            r2_record = FastqRecord()
            pair = Pair()
            while True:
                r1_record.read(r1_in)
                if not r1_record.identifier:
                    break
                r2_record.read(r2_in)
                pair.set_from_records(r1_record, r2_record)

                spats.process_pair(pair)

                summary = "{} :: {}".format(pair.identifier, pair.site if pair.has_site else pair.failure)
                if pair.r1.match_errors:
                    summary += " R1!: {}".format(pair.r1.match_errors)
                if pair.r1.adapter_errors:
                    summary += " R1A!: {}, adapter_len={}".format(pair.r1.adapter_errors, pair.r1._rtrim)
                if pair.r2.match_errors:
                    summary += " R2!: {}".format(pair.r2.match_errors)
                if pair.r2.adapter_errors:
                    summary += " R2A!: {}, adapter_len={}".format(pair.r2.adapter_errors, pair.r2._rtrim - 4)
                print summary

def diag_case():
    from spats_shape_seq import Spats, spats_config
    from spats_shape_seq.pair import Pair
    from spats_shape_seq.tests.test_pairs import cases
    from diagram import diagram
    spats_config.minimum_target_match_length = 8
    spats = Spats()
    spats.addMasks('RRRY', 'YYYR')
    spats.addTargets("test/5s/5s.fa")
    spats_config.debug = True
    spats._case_errors = False
    def run_case(case):
        pair = Pair()
        pair.set_from_data(case[0], case[1], case[2])
        spats.process_pair(pair)
        print diagram(pair)
        if case[3] != pair.site:
            spats._case_errors = True
            print "******* mismatch: {} != {}".format(case[3], pair.site)
    for case in cases:
        if case[0].startswith("*"):
            run_case(case)
    spats_config.debug = False
    if spats._case_errors:
        raise Exception("Case failed")

def test_compare_v102():
    from spats_shape_seq.v102 import compare_v102
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    compare_v102(bp + "s3",
                 bp + "5s/5S.fa",
                 'RRRY', 'YYYR',
                 bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                 bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                 bp + "5S-2p1-18x/spats_out")

def test_compare_v102():
    from spats_shape_seq.v102 import compare_v102
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    compare_v102(bp + "s5",
                 bp + "datasets/Shape_Seq_ligation/panel_RNAs_complete.fa",
                 'RRRY', 'YYYR',
                 bp + "datasets/Shape_Seq_ligation/data/KEW1_S1_L001_R1_001.fastq",
                 bp + "datasets/Shape_Seq_ligation/data/KEW1_S1_L001_R2_001.fastq",
                 bp + "datasets/Shape_Seq_ligation/spats_out")

def diagram_v102():
    from spats_shape_seq.v102 import diagram_case
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    diagram_case(bp + "s4/pairs.db",
                 14,
                 bp + "datasets/Shape_Seq_ligation/panel_RNAs_complete.fa",
                 'RRRY', 'YYYR')

def rc():
    from spats_shape_seq.util import reverse_complement
    print reverse_complement(sys.argv[2])

def makedb():
    db_path = sys.argv[2]
    targets_path = sys.argv[3]
    r1_path = sys.argv[4]
    r2_path = sys.argv[5]
    from spats_shape_seq.db import PairDB
    db = PairDB(db_path)
    db.show_progress_every = 200000
    db.load_and_index(targets_path, r1_path, r2_path)

def dbrun():
    db_path = sys.argv[2]
    run_name = sys.argv[3]
    from spats_shape_seq import Spats
    from spats_shape_seq.db import PairDB
    db = PairDB(db_path)
    s = Spats()
    s.run.writeback_results = True
    s.run.result_set_name = run_name
    #s.run.resume_processing = True
    s.process_pair_db(db)

def addv102():
    db_path = sys.argv[2]
    targets_path = sys.argv[3]
    out_path = sys.argv[4]
    from spats_shape_seq.db import PairDB
    db = PairDB(db_path)
    db.add_v102_comparison(targets_path, out_path)

def rdiff():
    db_path = sys.argv[2]
    rs1_name = sys.argv[3]
    rs2_name = sys.argv[4]
    from spats_shape_seq.db import PairDB
    db = PairDB(db_path)
    n1 = db.num_results(rs1_name)
    n2 = db.num_results(rs2_name)
    print "{}: {} results  /  {}: {} results".format(rs1_name, n1, rs2_name, n2)
    if not n1 or not n2:
        print "** Abort."
        exit(1)
    print "Diffs:"
    ours_only = []
    theirs_only = []
    differences = []
    for r in db.differing_results(rs1_name, rs2_name):
        if r[4] == -1:
            assert(r[7] != -1)
            theirs_only.append(r)
        elif r[7] == -1:
            ours_only.append(r)
        else:
            differences.append(r)
    all_lists = [ ours_only, theirs_only, differences ]
    for l in all_lists:
        reasons = {}
        for r in l:
            key = r[5] or r[8]
            assert(key)
            rlist = reasons.get(key)
            if not rlist:
                rlist = []
                reasons[key] = rlist
            rlist.append(r)
        for reason, rlist in reasons.iteritems():
            for r in rlist[:min(len(rlist), 10)]:
                print "  {}:{} ({}) -- {}:{} ({})   ({}: {} / {})".format(r[3] or 'x', r[4], r[5] or "OK",
                                                                          r[6] or 'x', r[7], r[8] or "OK",
                                                                          r[0], r[1], r[2])
            if len(rlist) > 0:
                print "... {} total.".format(len(rlist))
    print "{} total diffs.".format(sum(map(len, all_lists)))

def test_tags():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
    from spats_shape_seq import Spats
    s = Spats()
    from spats_shape_seq.tag import TagProcessor
    s.run._processor_class = TagProcessor

    #from spats_shape_seq.target import Targets
    #s.addTargets(bp + "5S.fa")
    s.addTarget("5s", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC")
    #s.addTarget("rc(5s)", "GTCCTTGGTGCCCGAGTCAGATGCCTGGCAGTTCCCTACTCTCGCATGGGGAGACCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTACGGCCGCCAGGCATCC")
    #s.addTarget("adapter_t", s.run.adapter_t)
    #s.addTarget("adapter_b", s.run.adapter_b)
    #s._targets._index_word_length = 8
    #s._targets._minimum_length = 8
    #s.addTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
    #s.addTarget("adapter_b_rc", reverse_complement(s.run.adapter_b))

    p = s._processor
    p.addTagTarget("5s", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC")
    p.addTagTarget("5s_rc", "GTCCTTGGTGCCCGAGTCAGATGCCTGGCAGTTCCCTACTCTCGCATGGGGAGACCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTACGGCCGCCAGGCATCC")
    from spats_shape_seq.util import reverse_complement
    p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
    p.addTagTarget("adapter_b", s.run.adapter_b)

    from spats_shape_seq.pair import Pair
    cases = [
        [ "1101:20069:1063", "TTTAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "TCCCACCTGACCCCATGCCGAACTCAGAAGTGAAA" ],
        [ "1101:11562:1050", "AAACGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG" ],
        [ "21189", "TTTGGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACCAAAAGATCGGAAGA" ],
        [ "1101:12888:8140", "GGATGTCCTTGGTGCCCGAGTCAGATGCCAGATCG", "GGCATCTGACTCGGGCACCAAGGACATACAGATCG" ],
        [ "18333", "GAGTGTCCTTGGTGCCCGAGTCAGTGGTAGATCGG", "ACCACTGACTCGGGCACCAAGGACACTCAGATCGG" ],
    ]

    pair = Pair()
    for case in cases:
        pair.set_from_data(case[0], case[1], case[2])
        s.process_pair(pair)

        print pair.r1.original_seq
        print pair.r1.tags
        print pair.r2.original_seq
        print pair.r2.tags
        print "-----------------------------"

    #s.run.skip_database = True
    #s.run.writeback_results = True
    #s.run.resume_processing = True
    #s.run.result_set_name = "lookup"
    #s.process_pair_data(bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
    #                    bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq")


def tags():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/"

    from spats_shape_seq.db import PairDB
    pair_db = PairDB(bp + "db/pairs.db")
    if True:
        print "Parsing to db..."
        pair_db.wipe()
        pair_db.add_targets_table(bp + "cotrans_single.fa")
        pair_db.parse(bp + "data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq",
                      bp + "data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq",
                      sample_size = 100000)

    from spats_shape_seq import Spats
    from spats_shape_seq.tag import TagProcessor
    from spats_shape_seq.util import reverse_complement
    s = Spats()
    s.run._processor_class = TagProcessor
    s.run.writeback_results = True
    s.run.result_set_name = "tags"
    s.run.num_workers = 1
    s.run.cotrans = True
    s.run.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'
    s.loadTargets(pair_db)

    s.run.allow_indeterminate = True
    s.run.allowed_target_errors = 2
    s.run.allowed_adapter_errors = 2

    p = s._processor
    for target in pair_db.targets():
        p.addTagTarget(target[0], target[1])
        p.addTagTarget(target[0] + "_rc", reverse_complement(target[1]))
    p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
    p.addTagTarget("adapter_b", s.run.adapter_b)
    if s.run.cotrans:
        p.addTagTarget("linker_cotrans", s.run.cotrans_linker)
        p.addTagTarget("linker_cotrans_rc", reverse_complement(s.run.cotrans_linker))

    s.process_pair_db(pair_db)
    rsid = pair_db.result_set_id_for_name(s.run.result_set_name)
    pair_db.count_tags(rsid)
    print pair_db.tag_counts(rsid)


def tquery():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/"

    from spats_shape_seq.db import PairDB
    pair_db = PairDB(bp + "db/pairs.db")
    print pair_db.results_matching(1, [ "linker_cotrans", "adapter" ], [ "match" ])

def tag_test():
    from spats_shape_seq import Spats
    s = Spats()
    s.run.cotrans = True
    s.run.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'
    from spats_shape_seq.partial import PartialFindProcessor
    s.run._processor_class = PartialFindProcessor

    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/"
    s.addTargets(bp + "cotrans_single.fa")

    from spats_shape_seq.pair import Pair
    pair = Pair()
    import cjb.util
    d = cjb.util.jsonAtPath("/tmp/spats_test.json")
    pair.set_from_data(str(d['id']), str(d['r1']), str(d['r2']))
    print "{}\n{} / {}".format(pair.identifier, pair.r1.original_seq, pair.r2.original_seq)
    s.process_pair(pair)
    if pair.has_site:
        print "{}: {} / {}".format(pair.target.name, pair.site, pair.right)
    else:
        print "FAIL: {}".format(pair.failure)

def cotrans_debug():
    from spats_shape_seq import Spats
    s = Spats()
    s.run.cotrans = True
    s.run.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'
    s.run._v102_compat = True
    s.run.minimum_target_match_length = 10

    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/"
    s.addTargets(bp + "cotrans_single.fa")

    from spats_shape_seq.pair import Pair
    pair = Pair()

    import cjb.util
    d = cjb.util.jsonAtPath("/tmp/spats_test.json")
    pair.set_from_data(str(d['id']), str(d['r1']), str(d['r2']))
    #c = ['683779', 'TCCGGTCCTTGGTGCCCGAGTCAGAAAAAAATAGAA', 'TCTATTTTTTTCTGACTCGGGCACCAAGGACCGGAA', 82, 71]
    #c = [ "1116:19486:8968", "TCCGGTCCTTGGTGCCCGAGTCAGTCCTTCCTCCTA", "GAGTCTATTTTTTTAGGAGGAAGGACTGACTCGGGC", 93, 68 ]
    #c = [ "301028", "AAGTGTCCTTGGTGCCCGAGTCAGAGATAGATCGGA", "ATCTCTGACTCGGGCACCAAGGACACTTAGATCGGA", 96, 92 ]
    #c = [ "31631284", "TTCAGTCCTTGGTGCCCGAGTCAGAGATAGATCGGA", "ATCTCTGACTCGGGCACCAATGACCGGAAGATCGGA", 96, 92 ]
    #c = [ "7232", "AGGTGTCCTTGGTGCCCGAGTCAGTAGCTAAGAAAT", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAG", -1, -1 ]
    #c = [ "16845404", "AAATGTCCTTGGTGCCCGAGTCAGACTGGTAGGAGT", "TCTTATAGGCGATGGAGTTCGCCATAAACGCTGCTT", -1, -1 ]
    #c = [ "24102328", "AAGCGTCCTTGGTGCCCGAGTCAGGAGTCATAGATC", "ATGACTCCTGACTCGGGCACCAAGGACGCTTAGATC", 46, 39 ]
    c = [ "51216106", "GGGTGTCCTTGGTGCCCGAGTCAGATTAGCTAAGCA", "AGCTAATCTGACTCGGGCACCAAGGACGCTGCTTAG", 41, 34 ]
    pair.set_from_data(c[0], c[1], c[2])
    print "{}\n{} / {}".format(pair.identifier, pair.r1.original_seq, pair.r2.original_seq)
    s.process_pair(pair)
    if pair.has_site:
        print "{}: {} / {}".format(pair.target.name, pair.site, pair.right)
    else:
        print "FAIL: {}".format(pair.failure)
    
if __name__ == "__main__":
    import sys
    globals()[sys.argv[1]]()

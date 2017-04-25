
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

def profile_run():
    bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
    spats(bp + "5s/5S.fa",
          bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
          bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
          bp + "t11",
          show_sites = False,
          max_pairs = 0)
    #Processed 64052 properly paired fragments, kept 4936/32073 (15.4%) treated, 8900/31979 (27.830764%) untreated

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

def spats(target, r1, r2, out, show_sites = True, max_pairs = 0):
    from spats_clean import Spats, spats_config
    s = Spats(target, out)
    s.setup()
    if show_sites:
        spats_config.show_id_to_site = True
    s.process_pair_data(r1, r2, max_pairs)
    if not show_sites:
        s.compute_profiles()
        s.write_reactivities()
              
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

if __name__ == "__main__":
    import sys
    globals()[sys.argv[1]]()

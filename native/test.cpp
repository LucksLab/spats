

#include "ats.hpp"
#include "lookup.hpp"
#include "pair.hpp"
#include "parse.hpp"
#include "r1lookup.hpp"
#include "r1tree.hpp"
#include "seq.hpp"
#include <map>
#include "db.hpp"



void
test_parse(const char * ftext)
{
    size_t len = strlen(ftext);
    Fragment f(ftext, len);
    std::string res = f.string();
    ATS_DEBUG("Res is: %s", res.c_str());
    ATS_DEBUG(" ft is: %s", ftext);
    ATS_ASSERT(res == std::string(ftext));
}

void
test_parse_fns(void)
{
    test_parse("ACGTGTGC");
    test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCC");
    test_parse("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG");
    test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCG");
    test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGCTCTGCTAGCTAGCATCGACGACGA");
    //test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGGGGGTTTTTTTAAAAAAACCCCCC"); // should assert for length
    test_parse("AANT");
    // test_parse("AANTY"); // should assert for invalid char
    ATS_DEBUG("s: %d", (int)sizeof(Fragment));
}

void
lins(FragmentStorage * t, const char * s)
{
    FragmentResult * fr = new FragmentResult();
    fr->m_fragment.parse(s, strlen(s));
    t->insert(fr);
    /* leak fr */
}

void
lfind(FragmentStorage * t, const char * s)
{
    Fragment f(s, 10);
    FragmentResult * r = t->find(&f);
    if (NULL != r)
        printf("Found: %s\n", r->m_fragment.string().c_str());
    else
        printf("No dice for: %s\n", s);
}

void
tfs(FragmentStorage * fs)
{
    lins(fs, "AAAACCCGGT");
    lins(fs, "CAAACCCGGT");
    lins(fs, "AGAACCCGGT");
    fs->dump();
    lfind(fs, "AGAACCCGGT");
    lfind(fs, "AGAACCCGAT");
    lfind(fs, "CAAACCCGGT");
}

void
tlookup()
{
    Lookup t(14);
    tfs(&t);
}

void
ttree()
{
    R1Tree t(10);
    tfs(&t);
}

#if 0
void
tr1l()
{
    Target t(1, "5s", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC");
    R1Lookup l(&t, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 31, 1);
    Fragment f("CAGGAACCACGGGCTCAGTAGATCGGAAGAG");
    FragmentResult * fr = l.find(&f);
    if (fr)
        printf("FR site: %d [%d]\n", fr->m_site, fr->m_errors);
    else
        printf("not found\n");
}

void
tf()
{
    Fragment f("CAGGAACCACGGGCTCAGTAGATCGGAAGAG", 31);
    printf("%s\n", f.string().c_str());
    f.del(0);
    printf("%s\n", f.string().c_str());
    f.insert(1, T_bits);
    printf("%s\n", f.string().c_str());
}
#endif

#define MASK_NO_MATCH 0
#define MASK_TREATED 1
#define MASK_UNTREATED 2

// RRRY = treated
// YYYR = untreated
int
mask_match(const char * handle)
{
    bool treated = true;
    bool untreated = true;
    char ch;
    for (int i = 0; i < 3; ++i) {
        ch = handle[i];
        if (ch == 'A' || ch == 'G')
            untreated = false;
        else if (ch == 'C' || ch == 'T')
            treated = false;
        else
            untreated = treated = false;
    }

    ch = handle[3];
    if (ch == 'A' || ch == 'G')
        treated = false;
    else if (ch == 'C' || ch == 'T')
        untreated = false;
    else
        untreated = treated = false;

    assert(!(treated  &&  untreated));
    return (treated ? MASK_TREATED : (untreated ? MASK_UNTREATED : MASK_NO_MATCH));
}

#define MAX_SITES 256
int g_total = 0;
int g_indeterminate = 0;
int g_matched = 0;
int g_matched_with_errors = 0;
int g_mask_failure = 0;
int g_treated_sites[MAX_SITES] = { 0 };
int g_untreated_sites[MAX_SITES] = { 0 };
R1Lookup * g_r1l = NULL;
R2Lookup * g_r2l = NULL;
Fragment * g_linker = NULL;
Target * g_target = NULL;
std::map< int, std::map< int, std::map < int, int > > > g_counters;
std::string g_adapter_t_rc;
int g_cotrans_minimum_length = 20;
pthread_mutex_t g_mutex;


bool
lookup_handler(Fragment * r1, Fragment * r2, const char * handle)
{
    //printf("R1: %s\nR2: %s\nH: %s\n", r1->string().c_str(), r2->string().c_str(), handle);
    ++g_total;
    //if (g_total > 20)
    //    exit(-1);
    //if (r1->has_errors() || r2->has_errors()) {
    // ++g_indeterminate;
    //    return true;
    //}
    FragmentResult * res = g_r1l->find(r1);
    if (res) {
        ++(res->m_count);
        if (res->m_errors)
            ++g_matched_with_errors;
        else
            ++g_matched;
        int site = res->m_site;
        if (-1 == site) {
            FragmentResult * res2 = g_r2l->find(r2, res->m_target);
            if (res2)
                site = res2->m_site;
        }
        else {
            // TODO: need to build & verify R2 frag
            site = res->m_site;
        }
        if (-1 != site) {
            switch (mask_match(handle)) {
            case MASK_TREATED:
                ++g_treated_sites[site]; break;
            case MASK_UNTREATED:
                ++g_untreated_sites[site]; break;
            default:
                ++g_mask_failure;
            }
        }
    }
    return true;// false;
}

void
t5s()
{
    Targets t;
    t.parse("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/5S.fa");
    g_r1l = new R1Lookup(&t, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 31, 0);
    g_r2l = new R2Lookup(35, 0);
    g_r2l->addTarget(t.target(0));
    fastq_parse_handler("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                        &lookup_handler);
    for (int i = 0; i < MAX_SITES; ++i) {
        if (g_treated_sites[i] > 0 || g_untreated_sites[i] > 0)
            printf("  %d: %d / %d\n", i, g_treated_sites[i], g_untreated_sites[i]);
    }
    printf("%d matched, %d w/errors, %d indeterminate, %d mask failure, %d total\n", g_matched, g_matched_with_errors, g_indeterminate, g_mask_failure, g_total);
    //g_r1l->dump();
}

void
tfa()
{
    Targets t;
    t.parse("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/Shape_Seq_ligation/panel_RNAs_complete.fa");
    printf("%d\n", t.size());
}

void
tpanel()
{
    Targets t;
    t.parse("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/Shape_Seq_ligation/panel_RNAs_complete.fa");
    g_r1l = new R1Lookup(&t, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 31, 0);
    g_r2l = new R2Lookup(35, 0);
    for (int i = 0; i < t.size(); ++i)
        g_r2l->addTarget(t.target(i));
    fastq_parse_handler("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/Shape_Seq_ligation/data/KEW1_S1_L001_R1_001.fastq",
                        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/Shape_Seq_ligation/data/KEW1_S1_L001_R2_001.fastq",
                        &lookup_handler);
    printf("%d matched, %d w/errors, %d indeterminate, %d mask failure, %d total\n", g_matched, g_matched_with_errors, g_indeterminate, g_mask_failure, g_total);
    for (int i = 0; i < MAX_SITES; ++i) {
        if (g_treated_sites[i] > 0 || g_untreated_sites[i] > 0)
            printf("  %d: %d / %d\n", i, g_treated_sites[i], g_untreated_sites[i]);
    }
    //g_r1l->dump();
}

bool
try_lookup_hit(FragmentResult * res, Fragment * r1, Fragment * r2, const char * handle)
{
    int linker_len = g_linker->len();
    int pair_len = r2->len();
    int L = res->m_site;
    int trim = res->m_trim;
    int site = -1;
    std::string r2_seq(r2->string());
    std::string tseq(g_target->seq());

    if (0 == trim) {
        int r2_match_len = g_r2l->r2_length();
        Fragment r2f(r2_seq.c_str(), r2_match_len);
        FragmentResult * r2_res = g_r2l->find(&r2f, g_target);
        if (NULL != r2_res) {
            site = r2_res->m_site;
        }
        else {
            //pair.failure = Failures.nomatch
            return false;
        }
    }
    else {
        site = L - (pair_len - linker_len - 4) + trim;
    }

    // now need to verify R2 contents
    // R2: [target][linker][handle][adapter]

    int target_match_len = std::min(pair_len, L - site);
    if (target_match_len <= 0  ||  r2_seq.substr(0, target_match_len) != tseq.substr(site, target_match_len)) {
        //pair.failure = Failures.match_errors
        return false;
    }

    if (target_match_len < pair_len) {
        int linker_match_len = std::min(linker_len, pair_len - target_match_len);
        if (r2_seq.substr(target_match_len, linker_match_len) != g_linker->string().substr(0, linker_match_len)) {
            //pair.failure = Failures.linker
            return false;
        }

        if (trim > 0  &&  r2_seq.substr(r2_seq.length() - trim, trim) != g_adapter_t_rc.substr(0, trim)) {
            //pair.failure = Failures.adapter_trim
            return false;
        }
        if (pair_len - target_match_len - linker_len - 4 > trim) {
            //pair.failure = Failures.adapter_trim
            return false;
        }
    }

    if (r1->has_errors() || r2->has_errors()) {
        ++g_indeterminate;
        return false;
    }

    int mask = mask_match(handle);
    if (MASK_NO_MATCH == mask) {
        ++g_mask_failure;
        return false;
    }

    pthread_mutex_lock(&g_mutex);
    {
        g_counters[L][site][mask] += 1;
    }
    pthread_mutex_unlock(&g_mutex);

    return true;
}

bool
cotrans_lookup_handler(Fragment * r1, Fragment * r2, const char * handle)
{
    //printf("R1: %s\nR2: %s\nH: %s\n", r1->string().c_str(), r2->string().c_str(), handle);
    ++g_total;
    if (0 == g_total % 1000000) {
        printf(".");
        fflush(stdout);
    }
    //if (g_total > 500) {
    //    printf("%s.%s / %s...\n", handle, r1->string().c_str(), r2->string().c_str());
    //}
    //    exit(-1);
    //if (r1->has_errors() || r2->has_errors()) {
    // ++g_indeterminate;
    //    return true;
    //}
    FragmentResult * res = g_r1l->find(r1);
    if (!res) {
        //pair.failure = Failures.nomatch
        return true;
    }
    if (res->m_site < g_cotrans_minimum_length) {
        //pair.failure = Failures.cotrans_min
        return true;
    }
    if (NULL != res->m_next) {
        FragmentResult * cur = res;
        while (cur != NULL) {
            FragmentResult * other = res;
            while (other != NULL) {
                if (other != cur  &&  other->m_trim == cur->m_trim) {
                    //pair.failure = Failures.multiple_R1
                    return true;
                }
                other = other->m_next;
            }
            cur = cur->m_next;
        }
    }

    while (NULL != res) {
        if (try_lookup_hit(res, r1, r2, handle)) {
            ++g_matched;
            return true;
        }
        res = res->m_next;
    }
    return true;
}


void
tcotrans()
{
    pthread_mutex_init(&g_mutex, NULL);
    g_adapter_t_rc = reverse_complement("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
    Targets t;
    t.parse("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/cotrans_single.fa");
    g_target = t.target(0);
    g_linker = new Fragment("CTGACTCGGGCACCAAGGAC", 20);
    g_r1l = new R1Lookup(&t, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 36, 0, g_linker);
    g_r2l = new R2Lookup(14, 0); // 14 = conservative guess for target->self_match_len; should not exceed pair_len - linker_len
    for (int i = 0; i < t.size(); ++i)
        g_r2l->addTarget(t.target(i));
#if 1
    fastq_parse_handler(
#if 1
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq",
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq",
#else
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/small_R1.fastq",
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/small_R2.fastq",
#endif
                &cotrans_lookup_handler);
    printf("%d matched, %d w/errors, %d indeterminate, %d mask failure, %d total\n", g_matched, g_matched_with_errors, g_indeterminate, g_mask_failure, g_total);
    for (int i = 0; i < MAX_SITES; ++i) {
        if (g_treated_sites[i] > 0 || g_untreated_sites[i] > 0)
            printf("  %d: %d / %d\n", i, g_treated_sites[i], g_untreated_sites[i]);
    }
    printf("{\n");
    for (int mask = MASK_TREATED; mask <= MASK_UNTREATED; ++mask) {
        printf("\"%s\" : {\n", mask == MASK_TREATED ? "RRRY" : "YYYR");
        for (int L = g_cotrans_minimum_length; L < g_target->n() + 1; ++L) {
            if (L > g_cotrans_minimum_length)
                printf(",\n");
            printf(" %d : [", L);
            for (int site = 0; site <= L; ++site) {
                printf(" %d%s", g_counters[L][site][mask], (L==site?"":","));
            }
            printf(" ]");
        }
        printf("}%s\n", mask == MASK_TREATED ? "," :"");
    }
    printf("}\n");
#else
    Case c("1116:19486:8968", "TCCGGTCCTTGGTGCCCGAGTCAGTCCTTCCTCCTA", "GAGTCTATTTTTTTAGGAGGAAGGACTGACTCGGGC");
    cotrans_lookup_handler(&c.r1, &c.r2, c.handle.c_str());

    for (int mask = MASK_TREATED; mask <= MASK_UNTREATED; ++mask) {
        for (int L = g_cotrans_minimum_length; L < g_target->n() + 1; ++L) {
            for (int site = 0; site <= L; ++site) {
                if (g_counters[L][site][mask]) {
                    printf("%s / %d / %d \n", (MASK_TREATED == mask ? "RRRY" : "YYYR"), L, site);
                    return;
                }
            }
        }
    }
    printf("NO MATCH\n");
#endif
}

void
tcotrans2()
{
    Spats s(true);
    s.addTargets("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/cotrans_single.fa");
    //PairDB * pdb = new PairDB("/Users/jbrink/tmp/bar.spats");
    //s.m_writeback = true;
    //s.m_writeback_db = pdb;
    //pdb->start_worker();
    s.run_fastq(
#if 1
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq",
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq"
#else
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/med_R1.fq",
        "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/med_R2.fq"
#endif
        );
    //pdb->commit_results();
    printf("Counts: %s\n\n%s\n", s.counters()->count_json().c_str(), s.counters()->site_json(s.cotrans_minimum_length()).c_str());
    //printf("Counts: %s\n", s.counters()->count_json().c_str());
    //printf("NW: %d\n", pdb->num_written());
}

void
tcotrans3()
{
    Spats s(true);
    s.addTargets("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/cotrans_single.fa");
    s.run_db("/Users/jbrink/tmp/foo.spats");
    s.m_writeback = true;
    printf("Counts: %s\n\n%s\n", s.counters()->count_json().c_str(), s.counters()->site_json(s.cotrans_minimum_length()).c_str());
}

int
tcase(int argc, char ** argv)
{
    Spats s(true);
    s.addTargets("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/cotrans_single.fa");
    Case c("", argv[1], argv[2]);
    s.run_case(&c);
    if (MASK_NO_MATCH == c.mask  ||  c.L <= 0  ||  c.site < 0) {
        printf("[ None, None, None ]\n");
        return 1;
    }
    else {
        printf("[ \"%s\", %d, %d ]\n", (MASK_TREATED == c.mask ? "RRRY" : "YYYR"), c.L, c.site);
        return 0;
    }
}

void
f(char ch)
{
    int r = ((ch >> 1) & 3);
    //if (3 == r  &&  ch == 'G')
    //    r = 2;
    printf("%c: 0x%x, %d\n", ch, ch, r);
}

extern void thash();

void
tdb()
{
    PairDB pdb("/Users/jbrink/tmp/foo.spats");
    pdb.test();
}

int
main(int argc, char ** argv)
{
    if (5 != argc) {
        fprintf(stderr, "Usage: cotrans [target_path] [r1_path] [r2_path] [spats_out_path]\n");
        return -1;
    }
    Spats s(true);
    s.addTargets(argv[1]);
    s.process_pair_data(argv[2], argv[3]);
    s.store(argv[4]);
    return 0;

    //tcotrans2();
    //tdb();
    //return tcase(argc, argv);
    //tcotrans2();
    //thash();
    //f('A'); f('C'); f('G'); f('T'); f('\n'); f('\r');
    //tpanel();
    //t5s();
    //tlookup();
    //ttree();
    //tr1l();
    //tf();
}

/*

ok, what's shortest way to test this, say on 5s:
* create Target class with name/seq
  * this should just have a string seq (for now)
  * include a function that pull any subsequence
* create R1Lookup class which:
  * takes a target and an adapter_b
  * creates the lookup table
* make it handle 1bp toggles and indels
  * this is not complicated...
* parse an R1 fastq file and find matches
  * for now, just report how many
* R2 lookup
  * and test
* handles
* cotrans basic test
- track down all discrepancies with .py...

*/

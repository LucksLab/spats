
#include "ats.hpp"
#include "lookup.hpp"
#include "pair.hpp"
#include "parse.hpp"
#include "r1lookup.hpp"
#include "r1tree.hpp"
#include "seq.hpp"


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
p2s(void)
{
    pairs_to_spin("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                  "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                  "/Users/jbrink/mos/tasks/1RwIBa/tmp/5s.spin");
    //pairs_to_spin("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R1.fastq",
    //              "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/data/EJS_6_F_10mM_NaF_Rep1_GCCAAT_R2.fastq",
    //              "/Users/jbrink/mos/tasks/1RwIBa/tmp/cotrans.spin");
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

int g_total = 0;
int g_indeterminate = 0;
int g_matched = 0;
int g_matched_with_errors = 0;
int g_sites[256] = { 0 };
R1Lookup * g_r1l = NULL;

bool
lookup_handler(Fragment * r1, Fragment * r2, const char * handle)
{
    //printf("R1: %s\nR2: %s\nH: %s\n", r1->string().c_str(), r2->string().c_str(), handle);
    ++g_total;
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
    }
    return true;// false;
}

void
tlp()
{
    Target t(1, "5s", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC");
    g_r1l = new R1Lookup(&t, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 31, 0);
    fastq_parse("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                "/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/5S-2p1-18x/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                &lookup_handler);
    printf("%d matched, %d w/errors, %d indeterminate, %d total\n", g_matched, g_matched_with_errors, g_indeterminate, g_total);
    //g_r1l->dump();
}

int
main(void)
{
    tlp();
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

*/

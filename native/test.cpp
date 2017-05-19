
#include "ats.hpp"
#include "pair.hpp"
#include "parse.hpp"
#include "r1tree.hpp"
#include "seq.hpp"


void
test_parse(const char * ftext)
{
    size_t len = strlen(ftext);
    Fragment f(ftext, len);
    std::string res = f.string(len);
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
tins(R1Tree * t, const char * s)
{
    FragmentResult * fr = new FragmentResult();
    fr->m_fragment.parse(s, strlen(s));
    t->insert(fr);
    /* leak fr */
}

void
tfind(R1Tree * t, const char * s)
{
    Fragment f(s, 10);
    FragmentResult * r = t->find(f);
    if (NULL != r)
        printf("Found: %s\n", r->m_fragment.string(10).c_str());
    else
        printf("No dice for: %s\n", s);
}

int
main(void)
{
    R1Tree t(10);
    tins(&t, "AAAACCCGGT");
    tins(&t, "CAAACCCGGT");
    tins(&t, "AGAACCCGGT");
    t.dump();
    tfind(&t, "AGAACCCGGT");
    tfind(&t, "AGAACCCGAT");
    tfind(&t, "CAAACCCGGT");
    return 0;
}

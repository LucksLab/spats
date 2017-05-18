
#include "ats.hpp"
#include "pair.hpp"
#include "seq.hpp"

void
test_parse(const char * ftext)
{
    int len = strlen(ftext);
    Fragment f;
    parse_to_fragment(ftext, f, len);
    std::string res = fragment_to_string(f, len);
    ATS_DEBUG("Res is: %s", res.c_str());
    ATS_DEBUG(" ft is: %s", ftext);
    ATS_ASSERT(res == std::string(ftext));
}

int
main(void)
{
    test_parse("ACGTGTGC");
    test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCC");
    test_parse("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG");
    test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCG");
    test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGCTCTGCTAGCTAGCATCGACGACGA");
    //test_parse("AAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGGGGGTTTTTTTAAAAAAACCCCCC"); // should assert for length
    test_parse("AANT");
    // test_parse("AANTY"); // should assert for invalid char
    return 0;
}

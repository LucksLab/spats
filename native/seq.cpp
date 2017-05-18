
#include <string.h>

#include "ats.hpp"
#include "seq.hpp"

#define SEQ_DEBUG ATS_VERBOSE

void
parse_to_fragment(const char * text, Fragment &f, int length)
{
    int word_size = sizeof(uint64_t);
    SEQ_DEBUG("WL: %d, MTL: %d", length >> 2, length >> 3);
    ATS_ASSERT(length <= (sizeof(Fragment) - word_size) << 2);  /* 2 bits / nt  ==>  4nt / byte */
    ATS_ASSERT(length <= word_size << 3);                       /* a bit per nt to flag as mistranscribed */

    memset(&f, 0, sizeof(Fragment));

    int frag_sel = 0;
    int frag_idx = 0;
    uint64_t bits;
    uint64_t mistranscribed = 0;

    /* it'd almost certainly be better to do this by parsing a byte at
     * a time based on a dictionary of string -> byte mapping */
    for (int idx = 0; idx < length; ++idx)
    {
        bits = 0;
        switch (text[idx])
        {
        case 'A': bits = A_bits; break;
        case 'C': bits = C_bits; break;
        case 'G': bits = G_bits; break;
        case 'T': bits = T_bits; break;

        case 'N':
            mistranscribed += (0x1 << idx);
            break;

        default:
            ATS_ASSERT_NOT_REACHED();
            break;
        }
        frag_idx = idx << 1;
        if (frag_idx < (word_size << 3))
        {
            frag_sel = 0;
        }
        else
        {
            frag_sel = 1;
            frag_idx -= (word_size << 3);
        }
        ATS_ASSERT(frag_idx < (word_size <<3));
        f[frag_sel] += (bits << frag_idx);
        SEQ_DEBUG("FS/FI: %d[%d] <- %c", frag_sel, frag_idx, text[idx]);
    }

    f[SPATS_FRAG_WORDS - 1] = mistranscribed;
}


std::string
fragment_to_string(Fragment &f, int length)
{
    int word_size = sizeof(uint64_t);
    std::string res;
    int frag_sel = 0;
    int frag_idx = 0;
    uint64_t mistranscribed = f[SPATS_FRAG_WORDS - 1];
    char nt;
    for (int idx = 0; idx < length; ++idx)
    {
        if ((mistranscribed >> idx) & 0x1) {
            res.append(1, 'N');
        }
        else {
            frag_idx = idx << 1;
            if (frag_idx < (word_size << 3))
            {
                frag_sel = 0;
            }
            else
            {
                frag_sel = 1;
                frag_idx -= (word_size << 3);
            }
            switch ((f[frag_sel] >> frag_idx) & 0x3)
            {
            case A_bits: nt = 'A'; break;
            case C_bits: nt = 'C'; break;
            case G_bits: nt = 'G'; break;
            case T_bits: nt = 'T'; break;
            }
            SEQ_DEBUG("FS/FI: %d[%d] -> %c", frag_sel, frag_idx, nt);
            res.append(1, nt);
        }
    }
    return res;
}

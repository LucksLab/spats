
#include <string.h>

#include "ats.hpp"
#include "seq.hpp"


#define SEQ_DEBUG ATS_VERBOSE

const uint64_t A_bits = 0x0;
const uint64_t C_bits = 0x1;
const uint64_t G_bits = 0x2;
const uint64_t T_bits = 0x3;

void
Fragment::parse(const char * text, int length)
{
    SEQ_DEBUG("WL: %d, MTL: %d", length >> 2, length >> 3);
    ATS_ASSERT(length <= (sizeof(m_words) << 2));  /* 2 bits / nt  ==>  4nt / byte */
    ATS_ASSERT(length <= (sizeof(m_errors) << 3));        /* a bit per nt to flag as mistranscribed */

    memset(&m_words, 0, sizeof(m_words));
    memset(&m_errors, 0, sizeof(m_errors));

    int frag_sel = 0;
    int frag_idx = 0;
    int word_bits = (sizeof(uint64_t) << 3);
    uint64_t bits;

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
            m_errors += (0x1 << idx);
            break;

        default:
            ATS_ASSERT_NOT_REACHED();
            break;
        }
        frag_idx += 2;
        if (frag_idx >= word_bits)
        {
            ++frag_sel;
            frag_idx -= word_bits;
        }
        ATS_ASSERT(frag_sel < sizeof(m_words));
        m_words[frag_sel] += (bits << frag_idx);
        SEQ_DEBUG("FS/FI: %d[%d] <- %c", frag_sel, frag_idx, text[idx]);
    }
}


std::string
Fragment::string(int length) const
{
    std::string res;
    int word_bits = (sizeof(uint64_t) << 3);
    int frag_sel = 0;
    int frag_idx = 0;
    char nt;
    for (int idx = 0; idx < length; ++idx)
    {
        frag_idx += 2;
        if ((m_errors >> idx) & 0x1) {
            res.append(1, 'N');
        }
        else {
            if (frag_idx >= word_bits)
            {
                ++frag_sel;
                frag_idx -= word_bits;
            }
            switch ((m_words[frag_sel] >> frag_idx) & 0x3)
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

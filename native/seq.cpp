
#include <string.h>
#include <map>

#include "ats.hpp"
#include "seq.hpp"


#define SEQ_DEBUG ATS_VERBOSE


void
Fragment::reset()
{
    memset(&m_words, 0, sizeof(m_words));
    memset(&m_errors, 0, sizeof(m_errors));
}

void
Fragment::parse(const char * text, size_t length)
{
    SEQ_DEBUG("WL: %d, MTL: %d", length >> 2, length >> 3);
    ATS_ASSERT(length <= (sizeof(m_words) << 2));     /* 2 bits / nt  ==>  4nt / byte */
    ATS_ASSERT(length <= (sizeof(m_errors) << 3));    /* a bit per nt to flag as mistranscribed */

    this->reset();

    int frag_sel = 0;
    int frag_idx = 0;
    int word_bits = (sizeof(uint64_t) << 3);
    uint64_t bits;

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
        if (frag_idx >= word_bits)
        {
            ++frag_sel;
            frag_idx -= word_bits;
        }
        ATS_ASSERT(frag_sel < sizeof(m_words));
        m_words[frag_sel] += (bits << frag_idx);
        frag_idx += 2;
        SEQ_DEBUG("FS/FI: %d[%d] <- %c", frag_sel, frag_idx, text[idx]);
    }
}


std::string
Fragment::string(size_t length) const
{
    std::string res;
    int word_bits = (sizeof(uint64_t) << 3);
    int frag_sel = 0;
    int frag_idx = 0;
    char nt = 0;
    for (int idx = 0; idx < length; ++idx)
    {
        if ((m_errors >> idx) & 0x1) {
            res.append(1, 'N');
        }
        else {
            if (frag_idx >= word_bits)
            {
                ++frag_sel;
                frag_idx -= word_bits;
            }
            nt = nt_bits_to_ch((m_words[frag_sel] >> frag_idx) & 0x3)
            SEQ_DEBUG("FS/FI: %d[%d] -> %c", frag_sel, frag_idx, nt);
            res.append(1, nt);
        }
        frag_idx += 2;
    }
    return res;
}

uint64_t
Fragment::at(int index) const
{
    int frag_sel = 0;
    int frag_idx = index << 1;
    int word_bits = (sizeof(uint64_t) << 3);
    while (frag_idx >= word_bits)
    {
        ++frag_sel;
        frag_idx -= word_bits;
    }
    return ((m_words[frag_sel] >> frag_idx) & 0x3);
}

bool
Fragment::equals(Fragment * other, size_t length) const
{
    int frag_sel = 0;
    size_t frag_idx = length << 1;
    int word_bits = (sizeof(uint64_t) << 3);
    while (frag_idx >= word_bits)
    {
        if (m_words[frag_sel] != other->m_words[frag_sel])
            return false;
        ++frag_sel;
        frag_idx -= word_bits;
    }
    for (int idx = 0; idx < frag_idx; idx += 2) {
        if ( ((m_words[frag_sel] >> idx) & 0x3) != ((other->m_words[frag_sel] >> idx) & 0x3) )
            return false;
    }
    if (m_errors != other->m_errors)
        return false;
    return true;
}

void
Fragment::clone(Fragment * other)
{
    memcpy(&m_words, other->m_words, sizeof(m_words));
    m_errors = other->m_errors;
}

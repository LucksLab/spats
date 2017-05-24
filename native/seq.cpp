

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
    m_length = 0;
}

void
Fragment::parse(const char * text, size_t in_length)
{
    int length = (-1 == in_length ? strlen(text) : in_length);
    SEQ_DEBUG("WL: %d, MTL: %d", length >> 2, length >> 3);
    ATS_ASSERT(length <= (sizeof(m_words) << 2));     /* 2 bits / nt  ==>  4nt / byte */
    ATS_ASSERT(length <= (sizeof(m_errors) << 3));    /* a bit per nt to flag as mistranscribed */

    this->reset();
    m_length = length;

    int frag_sel = 0;
    int frag_idx = 0;
    int word_bits = (sizeof(uint64_t) << 3);
    uint64_t bits;
    uint64_t curWord = 0LL;
    uint64_t errors = 0LL;
    char ch = 0;

    for (int idx = 0; idx < length; ++idx)
    {
        ch = text[idx];
        if (ch == 'N') {
            errors += (0x1LL << idx);
        }
        else {
            // fast translation of character byte to *_bits
            bits = ((ch >> 1) & 3);
        }
        if (frag_idx >= word_bits)
        {
            m_words[frag_sel] = curWord;
            ++frag_sel;
            frag_idx -= word_bits;
        }
        ATS_ASSERT(frag_sel < sizeof(m_words));
        curWord += (bits << frag_idx);
        frag_idx += 2;
        SEQ_DEBUG("FS/FI: %d[%d] <- %c", frag_sel, frag_idx, text[idx]);
    }
    m_words[frag_sel] = curWord;
}


std::string
Fragment::string() const
{
    std::string res;
    int word_bits = (sizeof(uint64_t) << 3);
    int frag_sel = 0;
    int frag_idx = 0;
    char nt = 0;
    for (int idx = 0; idx < m_length; ++idx)
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

void
Fragment::set(int index, uint64_t nt)
{
    int frag_sel = 0;
    int frag_idx = index << 1;
    int word_bits = (sizeof(uint64_t) << 3);
    while (frag_idx >= word_bits)
    {
        ++frag_sel;
        frag_idx -= word_bits;
    }
    uint64_t mask = (0x3LL << frag_idx);
    m_words[frag_sel] = (m_words[frag_sel] & (~mask)) | (nt << frag_idx);
}

void
Fragment::insert(int index, uint64_t nt)
{
    // doesn't work with more than one word yet
    assert(index <= 31);
    int frag_sel = 0;
    int frag_idx = index << 1;
    int word_bits = (sizeof(uint64_t) << 3);
    while (frag_idx >= word_bits)
    {
        ++frag_sel;
        frag_idx -= word_bits;
    }
    uint64_t mask = ((0x1LL << frag_idx) - 0x1LL);
    m_words[frag_sel] = ( (m_words[frag_sel] & mask)  |
                          (nt << frag_idx)            |
                          ((m_words[frag_sel] & (~mask)) << 2) );
    ++m_length;
}

void
Fragment::del(int index)
{
    // doesn't work with more than one word yet
    assert(index <= 31);
    int frag_sel = 0;
    int frag_idx = index << 1;
    int word_bits = (sizeof(uint64_t) << 3);
    while (frag_idx >= word_bits)
    {
        ++frag_sel;
        frag_idx -= word_bits;
    }
    uint64_t mask = ((0x1LL << frag_idx) - 0x1LL);
    m_words[frag_sel] = ( (m_words[frag_sel] & mask)  |
                          ((m_words[frag_sel] & (~mask)) >> 2) );
    --m_length;
}

bool
Fragment::equals(Fragment * other, int length, int start_index) const
{
    int use_length = (-1 == length ? m_length :length );
    if (other->m_length < use_length)
        return false;
    int frag_sel = 0;
    int start_frag_idx = start_index << 1;
    size_t frag_idx = use_length << 1;
    int word_bits = (sizeof(uint64_t) << 3);
    while (frag_idx >= word_bits)
    {
        if (m_words[frag_sel] != other->m_words[frag_sel])
            return false;
        ++frag_sel;
        frag_idx -= word_bits;
        start_frag_idx -= word_bits;
    }
    uint64_t mask = (1LL << frag_idx) - 1LL;
    if ((m_words[frag_sel] & mask) != (other->m_words[frag_sel] & mask))
        return false;
    if (m_errors != other->m_errors)
        return false;
    return true;
}

#define m1 0x5555555555555555LL
#define m2 0x6666666666666666LL

int
Fragment::hamming_distance(Fragment * other) const
{
    // TODO: deal with length > one word..
    register uint64_t a = m_words[0] ^ other->m_words[0];
    register uint64_t b = ((a & m1) | ((a & m2) >> 1));
    return __builtin_popcountll(b);
}

void
Fragment::clone(Fragment * other)
{
    memcpy(&m_words, other->m_words, sizeof(m_words));
    m_errors = other->m_errors;
    m_length = other->m_length;
}

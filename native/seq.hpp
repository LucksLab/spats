
#ifndef __SPATS_SEQ_HPP_INCLUDED__
#define __SPATS_SEQ_HPP_INCLUDED__


#include <string>
#include <inttypes.h>


#define SPATS_FRAG_WORDS 2


#define A_bits  0x0LL
#define C_bits  0x1LL
#define G_bits  0x3LL
#define T_bits  0x2LL


class Fragment
{
private:
    int m_length;
    uint64_t m_words[SPATS_FRAG_WORDS];
    uint64_t m_errors;
    void reset();
public:
    Fragment() { reset(); }
    Fragment(const char * text, size_t length = -1) { this->parse(text, length = -1); }
    bool has_errors() const { return (0LL != m_errors); }
    int len() const { return m_length; }
    void parse(const char * text, size_t length = -1);
    std::string string() const;
    const char * str() const { return string().c_str(); }
    uint64_t at(int index) const;
    void set(int index, uint64_t nt);
    void insert(int index, uint64_t nt);
    void del(int index);
    uint64_t key() const { return m_words[0]; }
    bool equals(Fragment * other, int length = -1, int start_index = 0) const;
    int hamming_distance(Fragment * other) const;
    void clone(Fragment * other);
};

inline
char
nt_bits_to_ch(uint64_t bits)
{
    switch (bits)
    {
    case A_bits: return 'A';
    case C_bits: return 'C';
    case G_bits: return 'G';
    case T_bits: return 'T';
    default: return '?';
    }
}

inline
uint64_t
nt_bits_from_ch(char ch)
{
    switch (ch)
    {
    case 'A': return A_bits;
    case 'C': return C_bits;
    case 'G': return G_bits;
    case 'T': return T_bits;
    }
    return -1;
}


#endif // __SPATS_SEQ_HPP_INCLUDED__

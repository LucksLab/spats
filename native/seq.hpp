
#ifndef __SPATS_SEQ_HPP_INCLUDED__
#define __SPATS_SEQ_HPP_INCLUDED__


#include <string>
#include <inttypes.h>

#define SPATS_FRAG_WORDS 2

class Fragment
{
private:
    uint64_t m_words[SPATS_FRAG_WORDS];
    uint64_t m_errors;
public:
    Fragment() {}
    Fragment(const char * text, int length) { this->parse(text, length); }
    void parse(const char * text, int length);
    std::string string(int length) const;
};


#endif // __SPATS_SEQ_HPP_INCLUDED__


#include "target.hpp"
#include <assert.h>

#define BUFFER_SIZE 2048

std::string
reverse_complement(std::string seq)
{
    char buffer[BUFFER_SIZE];
    int len = seq.size();
    assert(len < BUFFER_SIZE);
    const char * cstr = seq.c_str();
    char ch = 0;
    for (int i = 0; i < len; ++i) {
        ch = cstr[len - (i + 1)];
        switch (ch) {
        case 'A': ch = 'T'; break;
        case 'C': ch = 'G'; break;
        case 'G': ch = 'C'; break;
        case 'T': ch = 'A'; break;
        default: break;
        }
        buffer[i] = ch;
    }
    return std::string(buffer, len);
}

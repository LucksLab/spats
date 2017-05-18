
#ifndef __SPATS_SEQ_HPP_INCLUDED__
#define __SPATS_SEQ_HPP_INCLUDED__


#include <string>
#include <inttypes.h>


/*
 * first (N-1) words are encodings of nt's ('ACGT')
 * last word flags any transcription errors (eg, 'N')
 */
#define SPATS_FRAG_WORDS 3
typedef uint64_t Fragment[SPATS_FRAG_WORDS];

const uint64_t A_bits = 0x0;
const uint64_t C_bits = 0x1;
const uint64_t G_bits = 0x2;
const uint64_t T_bits = 0x3;

void
parse_to_fragment(const char * text, Fragment &f, int length);

std::string
fragment_to_string(Fragment &f, int length);



#endif // __SPATS_SEQ_HPP_INCLUDED__


#ifndef __SPATS_MASK_HPP_INCLUDED__
#define __SPATS_MASK_HPP_INCLUDED__


#include "ats.hpp"


#define MASK_NO_MATCH 0
#define MASK_TREATED 1
#define MASK_UNTREATED 2


// RRRY = treated
// YYYR = untreated
inline int
match_mask(const char * handle)
{
    bool treated = true;
    bool untreated = true;
    char ch;
    for (int i = 0; i < 3; ++i) {
        ch = handle[i];
        if (ch == 'A' || ch == 'G')
            untreated = false;
        else if (ch == 'C' || ch == 'T')
            treated = false;
        else
            untreated = treated = false;
    }

    ch = handle[3];
    if (ch == 'A' || ch == 'G')
        treated = false;
    else if (ch == 'C' || ch == 'T')
        untreated = false;
    else
        untreated = treated = false;

    ATS_ASSERT(!(treated  &&  untreated));
    return (treated ? MASK_TREATED : (untreated ? MASK_UNTREATED : MASK_NO_MATCH));
}


#endif // __SPATS_MASK_HPP_INCLUDED__

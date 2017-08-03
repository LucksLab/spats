
#ifndef __SPATS_COUNTERS_HPP_INCLUDED__
#define __SPATS_COUNTERS_HPP_INCLUDED__


#include <map>
#include <string.h>


#define COUNT_INDEX_TOTAL         0
#define COUNT_INDEX_MATCHED       1
#define COUNT_INDEX_MASK_FAILURE  2
#define COUNT_INDEX_INDETERMINATE 3
#define NUM_FIXED_COUNTS          4

#define COUNTER(cstruct,counterName) (cstruct->counts[COUNT_INDEX_##counterName])
#define INCR_COUNTER(cstruct,counterName) ++(cstruct->counts[COUNT_INDEX_##counterName])

class Counters
{
public:
    int counts[NUM_FIXED_COUNTS];
    int n;
    int n2;
    int * counters;

    Counters(int length) : n(length + 1), n2((length + 1) * (length + 1))
    {
        for (int i = 0; i < NUM_FIXED_COUNTS; ++i)
            counts[i] = 0;
        counters = new int[n2 << 1];
        memset(counters, 0, n2 << 1);
    }
    ~Counters() { delete [] counters; }

    inline void register_site(int mask, int L, int site) { ++counters[(mask == 2 ? n2 : 0) + (n * L) + site]; }

    void aggregate(Counters * other);
    std::string count_json();
    std::string site_json(int cotrans_min_length);
};



#endif // __SPATS_COUNTERS_HPP_INCLUDED__

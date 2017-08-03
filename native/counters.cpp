
#include <sstream>

#include "counters.hpp"
#include "mask.hpp"


void
Counters::aggregate(Counters * other)
{
    for (int i = 0; i < NUM_FIXED_COUNTS; ++i)
        counts[i] += other->counts[i];
    for (int i = 0; i < (n2 << 1); ++i)
        counters[i] += other->counters[i];
}

std::string
Counters::count_json()
{
    std::stringstream json;
    json << "[";
    for (int i = 0; i < NUM_FIXED_COUNTS; ++i) {
        if (i > 0)
            json << ",";
        json << " " << counts[i];
    }
    json << " ]";
    return json.str();
}

std::string
Counters::site_json(int cotrans_min_length)
{
    std::stringstream json;
    json << "{";
    for (int mask = MASK_TREATED; mask <= MASK_UNTREATED; ++mask) {
        json << "\"" << (mask == MASK_TREATED ? "RRRY" : "YYYR") << "\":{\n";
        for (int L = cotrans_min_length; L < n; ++L) {
            if (L > cotrans_min_length)
                json << ",";
            json << "\n" << L << ": [";
            for (int site = 0; site <= L; ++site) {
                json << counters[(mask == 2 ? n2 : 0) + (n * L) + site];
                if (L != site)
                    json << ",";
            }
            json << "]";
        }
        json << "}";
        if (mask == MASK_TREATED)
            json << ",";
    }
    json << "}";
    return json.str();
}

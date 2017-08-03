
#ifndef __SPATS_SPATS_HPP_INCLUDED__
#define __SPATS_SPATS_HPP_INCLUDED__


#include "r1lookup.hpp"
#include "seq.hpp"
#include "target.hpp"
#include "counters.hpp"
#include "mask.hpp"


class PairDB;


#define COUNT_INDEX_TOTAL         0
#define COUNT_INDEX_MATCHED       1
#define COUNT_INDEX_MASK_FAILURE  2
#define COUNT_INDEX_INDETERMINATE 3
#define NUM_COUNTS                4

class Case
{
public:
    int pair_id;
    std::string id;
    std::string handle;
    Fragment r1;
    Fragment r2;
    int mask;
    int L;
    int site;
    Case() : pair_id(0), mask(MASK_NO_MATCH), L(-1), site(-1) {}
    Case(const char * _id, const char * _r1, const char * _r2) : id(_id), pair_id(0), r1(&_r1[4]), r2(_r2), handle(_r1, 4), L(-1), site(-1), mask(MASK_NO_MATCH)  { }
    Case(int _id, const char * _r1, const char * _r2) : pair_id(_id), r1(&_r1[4]), r2(_r2), handle(_r1, 4), L(-1), site(-1), mask(MASK_NO_MATCH)  { }
    bool valid() const { return (L > 0) && (site >= 0) && (mask != MASK_NO_MATCH); }
};


class Spats
{
public:
    bool m_cotrans;

    int m_pair_len;
    int m_cotrans_minimum_length;

    Targets m_targets;
    R1Lookup * m_r1l;
    R2Lookup * m_r2l;
    Fragment * m_linker;
    Target * m_cotrans_target;
    std::string m_adapter_t_rc;
    std::string m_adapter_b;

    Counters * m_counters;

    bool m_writeback;
    PairDB * m_writeback_db;

    void setup();
    bool try_lookup_hit(FragmentResult * res, Fragment * r1, Fragment * r2, const char * handle, Counters * counters, Case * caseinfo);

public: // internal use only
    bool spats_handler(Fragment * r1, Fragment * r2, const char * handle, Counters * counters, Case * caseinfo = NULL);

public:

    Spats(bool cotrans = false) : m_cotrans(cotrans),  m_r1l(NULL), m_r2l(NULL), m_cotrans_target(NULL),
          m_pair_len(36), m_cotrans_minimum_length(20), m_counters(NULL), m_writeback(false), m_writeback_db(NULL)
    {
        m_adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        m_adapter_t_rc = reverse_complement("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
        m_linker = new Fragment("CTGACTCGGGCACCAAGGAC", 20);
    }

    ~Spats() { }

    Counters * counters() const { return m_counters; };

    void addTargets(const char * path);

    void run_fastq(const char * r1_path, const char * r2_path);
    void run_case(Case * c);
    void run_db(const char * db_path);
    PairDB * writeback_db() const { return (m_writeback ? m_writeback_db : NULL); }

    int cotrans_minimum_length() const { return m_cotrans_minimum_length; }

};



#endif // __SPATS_SPATS_HPP_INCLUDED__

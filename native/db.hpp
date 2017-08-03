
#ifndef __SPATS_DB_HPP_INCLUDED__
#define __SPATS_DB_HPP_INCLUDED__


#include <sqlite3.h>
#include "ats.hpp"
#include "spats.hpp"

#define CASE_QUEUE_LEN 0x10000

class PairDB
{
    const char * m_path;
    sqlite3 * m_handle;

    pthread_t m_worker_thread;
    pthread_mutex_t m_mutex;
    bool m_working;
    Case m_cases[CASE_QUEUE_LEN];
    int m_head;
    int m_tail;
    int m_num_written;

public: // internal use only
    void worker_fn();

public:
    PairDB(const char * path);
    void test();
    void run_cases(Spats * s);
    void write_result(Case * c);
    void start_worker();
    void submit_result(Case * c);
    void commit_results();
    int num_written() const { return m_num_written; }
};

#endif // __SPATS_DB_HPP_INCLUDED__

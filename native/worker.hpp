
#ifndef __SPATS_WORKER_HPP_INCLUDED__
#define __SPATS_WORKER_HPP_INCLUDED__

#include <pthread.h>
#include "parse.hpp"

#define FRAG_BUFFER_SIZE 40
#define CQ_SIZE 512
#define NUM_WORKERS 6

#if 0
# define WORKER_TRACE(theStr) printf(theStr);
#else
# define WORKER_TRACE(theStr)
#endif


struct WorkItem
{
    bool ready;
    char r1chars[FRAG_BUFFER_SIZE];
    char r2chars[FRAG_BUFFER_SIZE];
    int pair_id;
};


struct WorkContext
{
    pair_handler handler;
    Spats * spats;
    size_t fragment_len;
};

struct Worker
{
    WorkContext * context;
    int id;
    bool done;
    pthread_t thread;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    WorkItem items[CQ_SIZE];
    int start;
    int end;
    int empty_worker;
    int count;
    Counters * counters;
};

void *
worker_fn(void * arg);

#endif // __SPATS_WORKER_HPP_INCLUDED__

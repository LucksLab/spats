
#include <unistd.h>
#include <sys/time.h>

#include "worker.hpp"


void *
worker_fn(void * arg)
{
    Worker * w = (Worker *)arg;
    Fragment r1_frag;
    Fragment r2_frag;
    size_t fragment_len = 0;
    Spats * spats = w->context->spats;
    pair_handler handler = w->context->handler;

    for (int i = 0; i < CQ_SIZE; ++i) {
        w->items[i].r1chars[fragment_len] = 0;
        w->items[i].r2chars[fragment_len] = 0;
    }

    while (0 == fragment_len) {
        fragment_len = w->context->fragment_len;
        usleep(10);
    }

    while (true) {
        while (true) {
            WorkItem * wi = &w->items[w->start];
            if (!wi->ready)
                break;

            /* handle it */
            //printf("%s  --  %s\n", wi->r1chars, wi->r2chars);
            r1_frag.parse(&wi->r1chars[4], fragment_len - 4);
            r2_frag.parse(wi->r2chars, fragment_len);
            //printf("%s  --  %s\n", r1_frag.string().c_str(), r2_frag.string().c_str());
            wi->r1chars[4] = 0;
            if (spats)
                spats->spats_handler(&r1_frag, &r2_frag, wi->r1chars, w->counters);
            else
                handler(&r1_frag, &r2_frag, wi->r1chars);

            WORKER_TRACE("^");
            wi->ready = false;
            w->start = (w->start + 1) % CQ_SIZE;
        }
        if (w->done)
            break;
        //sleep(1);
        WORKER_TRACE("z");
        ++w->empty_worker;

        usleep(50); // this is highly tune-able based on the work done in handler...
    }
    return w;
}

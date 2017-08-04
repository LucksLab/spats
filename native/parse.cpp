
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <unistd.h>

#include "ats.hpp"
#include "seq.hpp"
#include "parse.hpp"
#include "worker.hpp"


/* as of 8/2017, when turning this on:
 *  - timing is nearly identical (13.5s, whether parsing only or doing full data analysis)
 *  - ~90% of the time in the parsing-only case is spent in fread
 * so, i don't think it's likely to optimize much further... 
 */
#define BENCHMARK_PARSING_ONLY 0




#define BUFFER_SIZE (1024 << 8)
#define MAX_LINE_SIZE 1024

#define FIND_EOL(cptr) do { if (*cptr == '\n' || *cptr == '\0') { break; } ++cptr; } while (1)


void
fastq_parse_driver(const char * r1_path, const char * r2_path, pair_handler handler, Spats * spats)
{
    FILE * r1 = fopen(r1_path, "rb");
    FILE * r2 = fopen(r2_path, "rb");

    char r1buf[BUFFER_SIZE + 1];
    char r2buf[BUFFER_SIZE + 1];
    r1buf[BUFFER_SIZE] = r2buf[BUFFER_SIZE] = 0;
    size_t r1r2read = 0;
    int buf_idx = 0;
    char r1line[MAX_LINE_SIZE];
    char r2line[MAX_LINE_SIZE];
    char * r1start;
    char * r1end;
    char * r2start;
    char * r2end;

    size_t skip = 0;
    size_t full_skip = 0;
    size_t skip_padding = 0;
    WorkContext context;
    context.handler = handler;
    context.spats = spats;
    context.fragment_len = 0;

    Worker workers[NUM_WORKERS];
    int num_workers = (spats ? NUM_WORKERS : 1);
    for (int i = 0; i < num_workers; ++i) {
        Worker * w = &workers[i];
        w->id = i;
        w->context = &context;
        w->start = w->end = 0;
        w->items[0].ready = false;
        if (spats)
            w->counters = new Counters(spats->counters()->n - 1);
        else
            w->counters = NULL;
#if !(BENCHMARK_PARSING_ONLY)
        pthread_mutex_init(&(w->mutex), NULL);
        pthread_cond_init(&(w->cond), NULL);
        pthread_create(&(w->thread), NULL, &worker_fn, (void *)w);
#endif
    }
    Worker * cur_worker = NULL;
    int cur_worker_idx = 0;
    WorkItem * cur_work_item = NULL;
    int worker_full = 0;
    int fragment_len = 0;
    int pair_idx = 0;
    int wrap_len = 0;

    /* first-pass set up */
    {
        r1r2read = fread(r1buf, 1, BUFFER_SIZE, r1);
        fread(r2buf, 1, BUFFER_SIZE, r2);
        if (r1r2read == 0)
            goto fastq_parse_done;
        buf_idx = 0;

        char * lineid_end = strchr(r1buf, '\n');
        context.fragment_len = fragment_len = (int)(strchr(&lineid_end[1], '\n') - lineid_end - 1);
        full_skip = (2) + (fragment_len + 1) + (lineid_end - r1buf - 2); // comment line, quality line, ID line
        skip_padding = full_skip + 2 * fragment_len;
        buf_idx = (int)(lineid_end - r1buf + 1);
        r1start = &r1buf[buf_idx];
        r2start = &r2buf[buf_idx];
    }

    while (true)
    {
        /* at this point, r1start and r2start should be pointed at the beginning of the next fragment */
        skip = fragment_len;
        ATS_ASSERT((int)(buf_idx + skip + 8) < (int)r1r2read);
        r1end = r1start + skip; // hotpath opt: since we know structure of fastq, we know when there won't be \n
        r2end = r2start + skip;
        while (*r1end >= 0x10) ++r1end; // hotpath opt: \0, \r, \n are all < 0x10, all ascii are > 0x10
        buf_idx += (r1end - r1start + 1); // include the +1 at the end
        r2end = r2start + (r1end - r1start); // match r2
        *r1end = 0; // null-terminate
        *r2end = 0;

        ATS_ASSERT((int)buf_idx < (int)r1r2read);

#if !(BENCHMARK_PARSING_ONLY)
        /* r1start, r2start now have null-terminated strings; process the line */
        while (true) {
            cur_worker = &workers[cur_worker_idx];
            if (++cur_worker_idx >= num_workers)
                cur_worker_idx = 0;
            if (cur_worker->done)
                goto fastq_parse_done;
            cur_work_item = &cur_worker->items[cur_worker->end];
            if (cur_work_item->ready) {
                ++worker_full;
                WORKER_TRACE("!");
            }
            else {
                break;
            }
        }

        cur_worker->end = (cur_worker->end + 1) % CQ_SIZE;

        memcpy(cur_work_item->r1chars, r1start, fragment_len);
        memcpy(cur_work_item->r2chars, r2start, fragment_len);
        cur_work_item->pair_id = ++pair_idx;
        ATS_VERBOSE("P: %d (%p) %20s / %20s\n", cur_work_item->pair_id, cur_work_item, cur_work_item->r1chars, cur_work_item->r2chars);
        cur_work_item->ready = true;
        WORKER_TRACE("v");
#endif

        /* now try to skip to the beginning of the next fragment, 4 lines down */
        skip = full_skip;

        if (buf_idx + skip_padding >= r1r2read) {
            /* special-case handling: out of data in the buffer, need to read in more and wrap the contents */
            wrap_len = (int)(r1r2read - buf_idx);
            memcpy(r1line, &r1buf[buf_idx], wrap_len);
            memcpy(r2line, &r2buf[buf_idx], wrap_len);

            r1r2read = fread(r1buf, 1, BUFFER_SIZE, r1);
            fread(r2buf, 1, BUFFER_SIZE, r2);
            if (r1r2read == 0)
                goto fastq_parse_done;

            buf_idx = (int)(skip - wrap_len);

            // copy over plenty to be sure the next fragment is in the line buffer
            memcpy(&r1line[wrap_len], r1buf, std::max(0, buf_idx) + (fragment_len << 1));
            memcpy(&r2line[wrap_len], r2buf, std::max(0, buf_idx) + (fragment_len << 1));

            r1end = r1start = &r1line[skip];
            r2end = r2start = &r2line[skip];
        }
        else {
            /* normal case, just continue on in the buffer */
            buf_idx += skip;
            ATS_ASSERT(buf_idx >= 0);
            r1end = r1start = &r1buf[buf_idx];
            r2end = r2start = &r2buf[buf_idx];
        }

        /* now, find the next r1/r2 start... */
        while (*r1end >= 0x10) ++r1end;
        skip = (r1end - r1start + 1); // include the +1 at the end
        buf_idx += skip;
        r1start += skip;
        r2start += skip;
        /* ...before continuing with the next line.*/

    } // while (true)


fastq_parse_done:
    fclose(r1);
    fclose(r2);

    int wempty = 0;
#if !(BENCHMARK_PARSING_ONLY)
    for (int i = 0; i < num_workers; ++i) {
        Worker * w = &workers[i];
        w->done = true;
        pthread_join(w->thread, NULL);
        wempty += w->empty_worker;
        if (spats)
            spats->counters()->aggregate(w->counters);
    }
#endif
    ATS_DEBUG("\n%d wfull, %d empty\n", worker_full, wempty);
}

void
fastq_parse_handler(const char * r1_path, const char * r2_path, pair_handler handler)
{
    fastq_parse_driver(r1_path, r2_path, handler, NULL);
}

void
fastq_parse_spats(const char * r1_path, const char * r2_path, Spats * spats)
{
    fastq_parse_driver(r1_path, r2_path, NULL, spats);
}

int
appx_number_of_fastq_pairs(const char * r1_path)
{
    FILE * r1 = fopen(r1_path, "rb");
    char r1buf[MAX_LINE_SIZE + 1];
    int nread = (int)fread(r1buf, 1, MAX_LINE_SIZE, r1);
    fseek(r1, 0, SEEK_END);
    long sz = ftell(r1);
    fclose(r1);
    if (0 == nread)
        return 0;
    char * fourthnl = strchr(strchr(strchr(strchr(r1buf, '\n') + 1, '\n') + 1, '\n') + 1, '\n') + 1;
    // the +1 is since first records tend to be short, and we'd rather underestimate than overestimate
    int frag_len = 1 + (int)(fourthnl - r1buf);
    return (int)((float)(sz) / (float)(frag_len));
}


void
fasta_parse(const char * fasta_path, target_handler handler)
{
    std::ifstream infile(fasta_path);
    std::string line;
    std::string name;
    while (true) {
        if (!getline(infile, line))
            break;
        if (0 == line.size())
            continue;
        if ('>' == line.at(0))
            name = line.substr(1);
        else
            handler(name.c_str(), line.c_str());
    }
}


void
portable_srandomdev()
{
    FILE * random_file = fopen("/dev/random", "r");
    char random_seed = getc(random_file);
    unsigned int seed = (getc(random_file) << 24) + (getc(random_file) << 16) + (getc(random_file) << 8) + getc(random_file);
    srandom(random_seed);
    fclose(random_file);
}

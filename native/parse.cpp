
#include <stdio.h>
#include <fstream>

#include "ats.hpp"
#include "seq.hpp"
#include "parse.hpp"

#define BUFFER_SIZE (1024 << 8)
#define MAX_LINE_SIZE 1024

#define FIND_EOL(cptr) do { if (*cptr == '\n' || *cptr == '\0') { break; } ++cptr; } while (1)


void
fastq_parse(const char * r1_path, const char * r2_path, pair_handler handler)
{
    FILE * r1 = fopen(r1_path, "rb");
    FILE * r2 = fopen(r2_path, "rb");

    char r1buf[BUFFER_SIZE + 1];
    char r2buf[BUFFER_SIZE + 1];
    r1buf[BUFFER_SIZE] = r2buf[BUFFER_SIZE] = 0;
    size_t r1r2read;
    int buf_idx = 0;
    char r1line[MAX_LINE_SIZE];
    char r2line[MAX_LINE_SIZE];
    char * r1start;
    char * r1end;
    char * r2start;
    char * r2end;
    int wrap = 0;
    size_t fragment_len = 0;

    int line_index = 0;
    size_t skip = 0;
    Fragment r1_frag;
    Fragment r2_frag;
    size_t skips[4];
    bool first_pass = true;

fastq_parse_read_chunk:
    {
        r1r2read = fread(r1buf, 1, BUFFER_SIZE, r1);
        fread(r2buf, 1, BUFFER_SIZE, r2);
        if (r1r2read == 0)
            goto fastq_parse_done;

        if (first_pass)
        {
            char * lineid_end = strchr(r1buf, '\n');
            fragment_len = strchr(&lineid_end[1], '\n') - lineid_end - 1;
            skips[0] = lineid_end - r1buf - 2;
            skips[1] = fragment_len - 1;
            skips[2] = 0;
            first_pass = false;
        }

        buf_idx = 0;
        r1start = &r1buf[0];
        r2start = &r2buf[0];

    fastq_parse_parse_lines:
        {
            skip = skips[line_index];
            if (wrap  ||  buf_idx + skip >= r1r2read)
                skip = 0;
            r1end = r1start + skip; // hotpath opt: since we know structure of fastq, we know when there won't be \n
            r2end = r2start + skip;
            while (*r1end >= 0x10) ++r1end; // hotpath opt: \0, \r, \n are all < 0x10, all ascii are > 0x10
            while (*r2end >= 0x10) ++r2end;
            buf_idx += (r1end - r1start + 1); // include the +1 at the end

            if (wrap)
            {
                size_t r1x = strlen(r1line);
                size_t r2x = strlen(r2line);
                memcpy(&r1line[r1x], r1start, r1end - r1start);
                memcpy(&r2line[r2x], r2start, r2end - r2start);
                r1line[r1x + r1end - r1start] = 0;
                r2line[r2x + r2end - r2start] = 0;
                r1start = &r1line[0];
                r2start = &r2line[0];
                wrap = 0;
            }
            else if (buf_idx >= r1r2read)
            {
                /* need to wrap; +1 to get the null terminator */
                memcpy(r1line, r1start, r1end - r1start + 1);
                memcpy(r2line, r2start, r2end - r2start + 1);
                wrap = 1;
                goto fastq_parse_read_chunk;
            }
            else
            {
                *r1end = 0;
                *r2end = 0;
            }

            /* r1start, r2start now have null-terminated strings; process the line */
            if (1 == line_index)
            {
                r1_frag.parse(r1start + 4, fragment_len - 4);
                r2_frag.parse(r2start, fragment_len);
                r1start[4] = 0;
                if (!handler(&r1_frag, &r2_frag, r1start))
                    goto fastq_parse_done;
            }
            line_index = ((line_index + 1) & 0x3);

            r1start = r1end + 1;
            r2start = r2end + 1;
            goto fastq_parse_parse_lines;
        }

    }

fastq_parse_done:
    fclose(r1);
    fclose(r2);
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
pairs_to_spin(const char * r1_path, const char * r2_path, const char * spin_path)
{
    FILE * r1 = fopen(r1_path, "rb");
    FILE * r2 = fopen(r2_path, "rb");
    FILE * spin = fopen(spin_path, "wb");

    char r1buf[BUFFER_SIZE + 1];
    char r2buf[BUFFER_SIZE + 1];
    r1buf[BUFFER_SIZE] = r2buf[BUFFER_SIZE] = 0;
    size_t r1read, r2read;
    char r1line[MAX_LINE_SIZE];
    char r2line[MAX_LINE_SIZE];
    char * r1start;
    char * r1end;
    char * r2start;
    char * r2end;
    int wrap = 0;
    int num_pairs = 0;
    size_t fragment_len = 0;

    int line_index = 0;
    Fragment r1_frag;
    Fragment r2_frag;

    while (1)
    {
        r1read = fread(r1buf, 1, BUFFER_SIZE, r1);
        r2read = fread(r2buf, 1, BUFFER_SIZE, r2);
        if (r1read == 0  || r2read == 0)
            break;
        r1start = &r1buf[0];
        r2start = &r2buf[0];
        while (1)
        {
            r1end = r1start;
            r2end = r2start;
            FIND_EOL(r1end);
            FIND_EOL(r2end);
            if (wrap)
            {
                size_t r1x = strlen(r1line);
                size_t r2x = strlen(r2line);
                memcpy(&r1line[r1x], r1start, r1end - r1start);
                memcpy(&r2line[r2x], r2start, r2end - r2start);
                r1line[r1x + r1end - r1start] = 0;
                r2line[r2x + r2end - r2start] = 0;
                r1start = &r1line[0];
                r2start = &r2line[0];
                wrap = 0;
            }
            else if (*r1end == '\0'  ||  *r2end == '\0')
            {
                /* need to wrap; +1 to get the null terminator */
                memcpy(r1line, r1start, r1end - r1start + 1);
                memcpy(r2line, r2start, r2end - r2start + 1);
                wrap = 1;
                break;
            }
            else
            {
                *r1end = 0;
                *r2end = 0;
            }
            {
                /* r1start, r2start now have null-terminated strings; process the line */
                if (1 == (line_index & 0x3))
                {
                    if (0 == fragment_len)
                        fragment_len = strlen(r1start);
                    r1_frag.parse(r1start, fragment_len);
                    r2_frag.parse(r2start, fragment_len);
                    fwrite(&r1_frag, 1, sizeof(r1_frag), spin);
                    fwrite(&r2_frag, 1, sizeof(r2_frag), spin);
                    ++num_pairs;
                }
                ++line_index;
            }
            r1start = r1end + 1;
            r2start = r2end + 1;
        } // while 1 (parsing lines from this chunked read)
    } // while 1 (reading from the input files)
    printf("\n%d pairs parsed.\n", num_pairs);
    fclose(spin);
    fclose(r1);
    fclose(r2);
}

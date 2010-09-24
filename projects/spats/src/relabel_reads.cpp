/*
 *  relabel_reads.cpp
 *  spats
 *
 *  Created by Cole Trapnell on 4/15/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#endif

#include <stdio.h>
#include <cassert>
#include <vector>
#include <cstring>
#include <cstdlib>

#include "common.h"
#include "reads.h"
#include "tokenize.h"
#include "qual.h"

bool fastq_db = true;

void format_qual_string(const string& orig_qual_str,
						string& out_qual_str)
{
	out_qual_str = orig_qual_str;
	for (size_t i = 0; i < orig_qual_str.size(); ++i)
	{
		out_qual_str[i] = charToPhred33(orig_qual_str[i], 
										solexa_quals, 
										phred64_quals);
	}
}

static const int nuc_A = (1);
static const int nuc_C = (1<<1);
static const int nuc_G = (1<<2);
static const int nuc_T = (1<<3);

char get_encoding(char seq_char)
{
    char mask_char = -1;
    
    switch (seq_char)
    {
        case 'A':
            mask_char = nuc_A;
            break;
        case 'C':
            mask_char = nuc_C;
            break;
        case 'G':
            mask_char = nuc_G;
            break;
        case 'T':
        case 'U':
            mask_char = nuc_T;
            break;
        case 'R': /* A or G */
            mask_char = nuc_A | nuc_G;
            break;
        case 'Y': /* C or T */
            mask_char = nuc_C | nuc_T;
            break;
        case 'S': /* G or C */
            mask_char = nuc_G | nuc_C;
            break;
        case 'W': /* A or T */
            mask_char = nuc_A | nuc_T;
            break;
        case 'K': /* G or T */
            mask_char = nuc_G | nuc_T;
            break;
        case 'M': /* A or C */
            mask_char = nuc_A | nuc_C;
            break;
        case 'B': /* C or G or T */
            mask_char = nuc_C | nuc_G | nuc_T;
            break;
        case 'D': /* A or G or T */
            mask_char = nuc_A | nuc_G | nuc_T;
            break;
        case 'H': /* A or C or T */
            mask_char = nuc_A | nuc_C | nuc_T;
            break;
        case 'V': /* A or C or G */
            mask_char = nuc_A | nuc_C | nuc_G;
            break;
        case 'N': /* any */
            mask_char = nuc_A | nuc_C | nuc_G | nuc_T;
            break;
        default:
            break;
    }
    
    return mask_char;
}


bool mask_matches(const string& read_seq, 
                  const vector<char>& mask)
{
    size_t min_len = min(read_seq.length(), mask.size());
    for (size_t i = 0; i < min_len; ++i)
    {
        char encoding = get_encoding(read_seq[i]);
        if (encoding == -1)
            return false;
        if (!(mask[i] & encoding))
            return false;
    }
    return true;
}

void relabel_reads(vector<FILE*> reads_files, 
                   vector<pair<vector<char>, FILE*> >& masks_files)
{	
	int num_reads_chucked = 0, num_reads = 0;
	int next_id = 0;
	for (size_t fi = 0; fi < reads_files.size(); ++fi)
	{
		Read read;
		FILE* fa = reads_files[fi];
		while(!feof(fa))
		{
			read.clear();
			
			// Get the next read from the file
			if (reads_format == FASTA)
			{
				if (!next_fasta_record(fa, read.name, read.seq))
					break;
			}
			else if (reads_format == FASTQ)
			{
				string orig_qual;
				if (!next_fastq_record(fa, read.name, read.seq, read.alt_name, orig_qual))
					break;
				format_qual_string(orig_qual, read.qual);
			}
			
			++num_reads;
			++next_id;
			
            vector<bool> mask_status(masks_files.size(), false);
            size_t matched = masks_files.size();
            for (size_t i = 0; i < masks_files.size(); ++i)
            {
                mask_status[i] = mask_matches(read.seq, masks_files[i].first);
                if (mask_status[i])
                {
                    matched = i;
                    break;
                }
            }
            
            //bool matched_treated = mask_matches(read.seq, treated_mask);
            //bool matched_untreated = mask_matches(read.seq, untreated_mask);
            
            if (matched != masks_files.size())
            {
                fprintf(masks_files[matched].second,
                        "%d\n",
                        next_id); 
            }
            else 
            {
                num_reads_chucked++;
                continue;
            }

            if (!fastq_db)
            {
                if (reads_format == FASTA)
                    printf(">%s\n%s\n", read.name.c_str(), read.seq.c_str());
                else if (reads_format == FASTQ)
                    printf("@%s\n%s\n+\n%s\n", 
                           read.name.c_str(), read.seq.c_str(),read.qual.c_str());
            }
            else
            {
                if (reads_format == FASTA)
                {
                    printf("@%d\n%s\n+%s\n%s\n",
                           next_id,
                           read.seq.c_str(),
                           read.name.c_str(),
                           string(read.seq.length(), 'I').c_str());
                }
                else if (reads_format == FASTQ)
                {
                    printf("@%d\n%s\n+%s\n%s\n",
                           next_id,
                           read.seq.c_str(),
                           read.name.c_str(),
                           read.qual.c_str());
                }
            }
		}
	}
    fprintf(stderr, "Kept %d of %d reads\n", next_id - num_reads_chucked, next_id);
}

void init_handle_mask(const string& handle_iupac_seq, vector<char>& mask)
{
    for (size_t i = 0; i < handle_iupac_seq.length(); ++i)
    {
        char mask_char = get_encoding(handle_iupac_seq[i]);
        if (mask_char == -1)
        {
            fprintf(stderr, "Error: found non-IUPAC symbol in handle string\n");
            exit(1);
            break;
        }
        
        mask.push_back(mask_char); 
    }
}


void print_usage()
{
    fprintf(stderr, "Usage:   relabel_reads <reads1.fa/fq,...,readsN.fa/fq> [handle1] [handle2] .. [handleN]\n");
}


int main(int argc, char *argv[])
{
	//fprintf(stderr, "relabel_reads v%s\n", PACKAGE_VERSION); 
	//fprintf(stderr, "---------------------------\n");
	
	int parse_ret = parse_options(argc, argv, print_usage);
	if (parse_ret)
		return parse_ret;
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string reads_file_list = argv[optind++];
    
    vector<pair<vector<char>, FILE*> > masks_files;
    
    while(optind < argc)
    {
        string handle_seq = argv[optind++];
        vector<char> mask;
        
        init_handle_mask(handle_seq, mask);
        string read_id_filename = handle_seq + ".ids";
        FILE* read_id_file = fopen(read_id_filename.c_str(), "w");
        if (read_id_file == NULL)
        {
            fprintf(stderr, 
                    "Error could not open read id file %s for writing\n", 
                    read_id_filename.c_str()); 
            exit(1);
        }
        masks_files.push_back(make_pair(mask, read_id_file));
    }
    
    
	vector<string> reads_file_names;
    vector<FILE*> reads_files;
    tokenize(reads_file_list, ",",reads_file_names);
    for (size_t i = 0; i < reads_file_names.size(); ++i)
    {
        FILE* seg_file = fopen(reads_file_names[i].c_str(), "r");
        if (seg_file == NULL)
        {
            fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                    reads_file_names[i].c_str());
            exit(1);
        }
        reads_files.push_back(seg_file);
    }
	
    for (size_t i = 0; i < masks_files.size(); i++)
    {
        for (size_t j = 0; j < masks_files.size(); ++j)
        {
            if (i == j)
                continue;
            
            const vector<char>& mask_i = masks_files[i].first;
            const vector<char>& mask_j = masks_files[j].first;
            
            size_t smaller_mask_idx = min(mask_i.size(), mask_j.size());
            bool compatible_masks = false;
            for (size_t k = 0; k < smaller_mask_idx; ++k)
            {
                if (mask_i[k] & mask_j[k])
                {
                    compatible_masks = true;
                    break;
                }
            }
            if (compatible_masks)
            {
                fprintf(stderr, "Error: some reads might satisfy both handle %lu and handle %lu sequences\n", i, j);
                exit(1);
            }
        }
    }
    
	// Only print to standard out the good reads
	relabel_reads(reads_files, masks_files);
	
	return 0;
}
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

struct FragmentMaskSet
{
    FragmentMaskSet(const vector<char>& m, FILE* h, FILE* n)
    : mask(m), handle_reads(h), nonhandle_reads(n) {}
    vector<char> mask;
    FILE* handle_reads;
    FILE* nonhandle_reads;
};

void relabel_reads(vector<FILE*> handle_reads_files, 
                   vector<FILE*> nonhandle_reads_files,
                   vector<FragmentMaskSet>& masks)
{	
    
    if (handle_reads_files.size() != nonhandle_reads_files.size())
    {
        fprintf (stderr, "Error: some read files appear to be missing!\n");
    }   
	int num_fragments_chucked = 0, num_fragments = 0;
	int next_id = 0;
	for (size_t fi = 0; fi < handle_reads_files.size(); ++fi)
	{
		Read handle_read;
        Read nonhandle_read;
		FILE* handle_fa = handle_reads_files[fi];
        FILE* nonhandle_fa = nonhandle_reads_files[fi];
        
		while(!feof(handle_fa))
		{
			handle_read.clear();
			
			// Get the next read from the file
			if (reads_format == FASTA)
			{
				if (!next_fasta_record(handle_fa, handle_read.name, handle_read.seq))
					break;
                if (!next_fasta_record(nonhandle_fa, nonhandle_read.name, nonhandle_read.seq))
					break;
			}
			else if (reads_format == FASTQ)
			{
				string orig_qual;
				if (!next_fastq_record(handle_fa, handle_read.name, handle_read.seq, handle_read.alt_name, orig_qual))
					break;
				format_qual_string(orig_qual, handle_read.qual);
                
				if (!next_fastq_record(nonhandle_fa, nonhandle_read.name, nonhandle_read.seq, nonhandle_read.alt_name, orig_qual))
					break;
				format_qual_string(orig_qual, nonhandle_read.qual);
			}
			
			++num_fragments;
			++next_id;
			
            vector<bool> mask_status(masks.size(), false);
            size_t matched = masks.size();
            for (size_t i = 0; i < masks.size(); ++i)
            {
                mask_status[i] = mask_matches(handle_read.seq, masks[i].mask);
                if (mask_status[i])
                {
                    matched = i;
                    break;
                }
            }
            
            //bool matched_treated = mask_matches(read.seq, treated_mask);
            //bool matched_untreated = mask_matches(read.seq, untreated_mask);
            
            FILE* handle_f_out = NULL;
            FILE* nonhandle_f_out = NULL;
            
            if (matched != masks.size() && 
                handle_read.seq.length() > (masks[matched].mask.size() + 4) &&
                nonhandle_read.seq.length() > (masks[matched].mask.size() + 4))
            {
                handle_f_out = masks[matched].handle_reads;
                nonhandle_f_out = masks[matched].nonhandle_reads;
            }
            else 
            {
                num_fragments_chucked++;
                continue;
            }

            handle_read.seq = handle_read.seq.substr(masks[matched].mask.size());
            if (!handle_read.qual.empty())
            {
                handle_read.qual = handle_read.qual.substr(masks[matched].mask.size());   
            }
            
            if (nonhandle_read.seq.length() > masks[matched].mask.size())
            {
                nonhandle_read.seq = nonhandle_read.seq.substr(0, nonhandle_read.seq.length() - masks[matched].mask.size());
                if (!nonhandle_read.qual.empty())
                {
                    nonhandle_read.qual = nonhandle_read.qual.substr(0, nonhandle_read.qual.length() - masks[matched].mask.size());   
                } 
            }
            
            if (!fastq_db)
            {
                if (reads_format == FASTA)
                {
                    fprintf(handle_f_out, ">%s\n%s\n", handle_read.name.c_str(), handle_read.seq.c_str());
                    fprintf(nonhandle_f_out, ">%s\n%s\n", nonhandle_read.name.c_str(), nonhandle_read.seq.c_str());
                }
                else if (reads_format == FASTQ)
                {
                    fprintf(handle_f_out, "@%s\n%s\n+\n%s\n", 
                            handle_read.name.c_str(), handle_read.seq.c_str(),handle_read.qual.c_str());
                    fprintf(nonhandle_f_out, "@%s\n%s\n+\n%s\n", 
                            nonhandle_read.name.c_str(), nonhandle_read.seq.c_str(), nonhandle_read.qual.c_str());
                }
            }
            else
            {
                if (reads_format == FASTA)
                {
                    fprintf(handle_f_out,
                            "@%d\n%s\n+%s\n%s\n",
                            next_id,
                            handle_read.seq.c_str(),
                            handle_read.name.c_str(),
                            string(handle_read.seq.length(), 'I').c_str());
                    fprintf(nonhandle_f_out,
                            "@%d\n%s\n+%s\n%s\n",
                            next_id,
                            nonhandle_read.seq.c_str(),
                            nonhandle_read.name.c_str(),
                            string(nonhandle_read.seq.length(), 'I').c_str());
                }
                else if (reads_format == FASTQ)
                {
                    fprintf(handle_f_out,
                            "@%d\n%s\n+%s\n%s\n",
                            next_id,
                            handle_read.seq.c_str(),
                            handle_read.name.c_str(),
                            handle_read.qual.c_str());
                    fprintf(nonhandle_f_out,
                            "@%d\n%s\n+%s\n%s\n",
                            next_id,
                            nonhandle_read.seq.c_str(),
                            nonhandle_read.name.c_str(),
                            nonhandle_read.qual.c_str());
                }
            }
		}
	}
    fprintf(stderr, "Kept %d of %d reads\n", next_id - num_fragments_chucked, next_id);
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
    fprintf(stderr, "Usage:   relabel_reads <reads1_1.fa/fq,...,readsN_1.fa/fq> <reads1_2.fa/fq,...,readsN_2.fa/fq> [handle1] [handle2] .. [handleN]\n");
}


int main(int argc, char *argv[])
{
	//fprintf(stderr, "relabel_reads v%s\n", PACKAGE_VERSION); 
	//fprintf(stderr, "---------------------------\n");
	
    fastq_db = true;
    
	int parse_ret = parse_options(argc, argv, print_usage);
	if (parse_ret)
		return parse_ret;
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string handle_reads_file_list = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string nonhandle_reads_file_list = argv[optind++];
    
    vector<pair<vector<char>, string> > masks;
    
    while(optind < argc)
    {
        string handle_seq = argv[optind++];
        vector<char> mask;
        
        init_handle_mask(handle_seq, mask);
        masks.push_back(make_pair(mask, handle_seq));
    }
    
    
	vector<string> handle_reads_file_names;
    vector<FILE*> handle_reads_files;
    tokenize(handle_reads_file_list, ",",handle_reads_file_names);
    for (size_t i = 0; i < handle_reads_file_names.size(); ++i)
    {
        FILE* seg_file = fopen(handle_reads_file_names[i].c_str(), "r");
        if (seg_file == NULL)
        {
            fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                    handle_reads_file_names[i].c_str());
            exit(1);
        }
        handle_reads_files.push_back(seg_file);
    }
    
    vector<string> nonhandle_reads_file_names;
    vector<FILE*> nonhandle_reads_files;
    tokenize(nonhandle_reads_file_list, ",",nonhandle_reads_file_names);
    for (size_t i = 0; i < nonhandle_reads_file_names.size(); ++i)
    {
        FILE* seg_file = fopen(nonhandle_reads_file_names[i].c_str(), "r");
        if (seg_file == NULL)
        {
            fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                    nonhandle_reads_file_names[i].c_str());
            exit(1);
        }
        nonhandle_reads_files.push_back(seg_file);
    }
    
    vector<FragmentMaskSet> masks_files;
	
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < masks.size(); ++j)
        {
            if (i == j)
                continue;
            
            const vector<char>& mask_i = masks[i].first;
            const vector<char>& mask_j = masks[j].first;
            
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
        
        string handle_out_name = output_dir + "/" + masks[i].second + "_1.fq";
        string nonhandle_out_name = output_dir + "/" + masks[i].second + "_2.fq";
        
        FILE* handle_reads_out = fopen(handle_out_name.c_str(), "w");
        FILE* nonhandle_reads_out = fopen(nonhandle_out_name.c_str(), "w");
        
        if (!handle_reads_out)
        {
            fprintf(stderr, "Error: cannot open reads file %s for writing\n",
                    handle_out_name.c_str());
            exit(1);
        }
        
        if (!nonhandle_reads_out)
        {
            fprintf(stderr, "Error: cannot open reads file %s for writing\n",
                    nonhandle_out_name.c_str());
            exit(1);
        }
        
        FragmentMaskSet m(masks[i].first, handle_reads_out, nonhandle_reads_out);
        masks_files.push_back(m);
    }
    
    if (masks.empty())
    {
        string handle_out_name = output_dir + "/" + "NOMASK_1.fq";
        string nonhandle_out_name = output_dir + "/" + "NOMASK_2.fq";
        
        FILE* handle_reads_out = fopen(handle_out_name.c_str(), "w");
        FILE* nonhandle_reads_out = fopen(nonhandle_out_name.c_str(), "w");
        
        if (!handle_reads_out)
        {
            fprintf(stderr, "Error: cannot open reads file %s for writing\n",
                    handle_out_name.c_str());
            exit(1);
        }
        
        if (!nonhandle_reads_out)
        {
            fprintf(stderr, "Error: cannot open reads file %s for writing\n",
                    nonhandle_out_name.c_str());
            exit(1);
        }
        
        vector<char> mask;
        init_handle_mask("", mask);
        masks.push_back(make_pair(mask, ""));
        
        FragmentMaskSet m(masks[0].first, handle_reads_out, nonhandle_reads_out);
        masks_files.push_back(m);
    }
    
	// Only print to standard out the good reads
	relabel_reads(handle_reads_files, nonhandle_reads_files, masks_files);
	
	return 0;
}
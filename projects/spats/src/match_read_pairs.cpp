/*
 *  match_read_pairs.cpp
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

void print_paired_reads(FILE* left_reads_in_file,
                        FILE* right_reads_in_file,
                        FILE* left_reads_out_file,
                        FILE* right_reads_out_file)
{	
    Read curr_left_read;
    uint32_t curr_left_read_id = 0;
    
    Read curr_right_read;
    uint32_t curr_right_read_id = 0;
    
    while(curr_left_read_id != 0xFFFFFFFF && curr_right_read_id != 0xFFFFFFFF)
    {
        if (curr_left_read_id &&
            curr_left_read_id == curr_right_read_id)
        {
            fprintf(left_reads_out_file,
                    "@%s\n%s\n+%s\n%s\n",
                    curr_left_read.name.c_str(),
                    curr_left_read.seq.c_str(),
                    curr_left_read.alt_name.c_str(),
                    curr_left_read.qual.c_str());
            
            fprintf(right_reads_out_file,
                    "@%s\n%s\n+%s\n%s\n",
                    curr_right_read.name.c_str(),
                    curr_right_read.seq.c_str(),
                    curr_right_read.alt_name.c_str(),
                    curr_right_read.qual.c_str());
        }
        
        if (curr_left_read_id >= curr_right_read_id)
        {
            curr_right_read.clear();
            if (!next_fastq_record(right_reads_in_file, 
                                   curr_right_read.name, 
                                   curr_right_read.seq, 
                                   curr_right_read.alt_name, 
                                   curr_right_read.qual))
                break;
            

        }
        
        if (curr_left_read_id <= curr_right_read_id)
        {
            curr_left_read.clear();
            if (!next_fastq_record(left_reads_in_file, 
                                   curr_left_read.name, 
                                   curr_left_read.seq, 
                                   curr_left_read.alt_name, 
                                   curr_left_read.qual))
                break;
            

        }
        
        curr_left_read_id = atoi(curr_left_read.name.c_str());
        if (curr_left_read_id == 0)
        {
            fprintf(stderr, "Error: bad read IID encountered %s\n", 
                    curr_right_read.name.c_str());
            exit(1);
        }
        curr_right_read_id = atoi(curr_right_read.name.c_str());
        if (curr_left_read_id == 0)
        {
            fprintf(stderr, "Error: bad read IID encountered %s\n", 
                    curr_left_read.name.c_str());
            exit(1);
        }
    }
}

void print_usage()
{
    fprintf(stderr, "Usage:   match_read_pairs <reads_1.fq> <reads_2.fq> <reads_1_matched.fq> <reads_2_matched.fq>\n");
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
    
    string left_reads_in_name = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string right_reads_in_name = argv[optind++];
    
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string left_reads_out_name = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string right_reads_out_name = argv[optind++];
    
    FILE* left_reads_in_file = fopen(left_reads_in_name.c_str(), "r");
    if (left_reads_in_file == NULL)
    {
        fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                left_reads_in_name.c_str());
        exit(1);
    }
    
    FILE* right_reads_in_file = fopen(right_reads_in_name.c_str(), "r");
    if (right_reads_in_file == NULL)
    {
        fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                right_reads_in_name.c_str());
        exit(1);
    }
    
    FILE* left_reads_out_file = fopen(left_reads_out_name.c_str(), "w");
    if (left_reads_out_file == NULL)
    {
        fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                left_reads_out_name.c_str());
        exit(1);
    }
    
    FILE* right_reads_out_file = fopen(right_reads_out_name.c_str(), "w");
    if (right_reads_out_file == NULL)
    {
        fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                right_reads_out_name.c_str());
        exit(1);
    }
    
    print_paired_reads(left_reads_in_file,
                       right_reads_in_file,
                       left_reads_out_file,
                       right_reads_out_file);

	return 0;
}
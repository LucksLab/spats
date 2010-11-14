/*
 *  common.cpp
 *  Spats
 *
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <getopt.h>
#include "common.h"

using namespace std;

int segment_length = 25;
int segment_mismatches = 2;


ReadFormat reads_format = FASTQ;

bool verbose = false;

int max_multihits = 40;

string output_dir = "spats_out";


bool solexa_quals = false;
bool phred64_quals = false;

extern void print_usage();

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */

int parseInt(int lower, const char *errmsg, void (*print_usage)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            cerr << errmsg << endl;
            print_usage();
            exit(1);
        }
        return (int32_t)l;
    }
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static float parseFloat(float lower, float upper, const char *errmsg, void (*print_usage)()) {
    float l;
    l = (float)atof(optarg);
	
    if (l < lower) {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    if (l > upper)
    {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    return l;
	
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}



const char *short_options = "";

#define OPT_FASTA						33
#define OPT_FASTQ						34
#define OPT_VERBOSE						37
#define OPT_OUTPUT_DIR					42
#define OPT_MAX_MULTIHITS				46
#define OPT_SEGMENT_LENGTH				50
#define OPT_SEGMENT_MISMATCHES			51
#define OPT_SOLEXA_QUALS			    63
#define OPT_PHRED64_QUALS				64
#define OPT_NO_RELABEL                  65

static struct option long_options[] = {
{"fasta",				no_argument,		0,	OPT_FASTA},
{"fastq",				no_argument,		0,	OPT_FASTQ},
{"verbose",				no_argument,		0,	OPT_VERBOSE},
{"output-dir",			required_argument,	0,	OPT_OUTPUT_DIR},
{"segment-length",		required_argument,	0,  OPT_SEGMENT_LENGTH},
{"segment-mismatches",	required_argument,	0,  OPT_SEGMENT_MISMATCHES},
{"solexa-quals",		no_argument,		0,	OPT_SOLEXA_QUALS},
{"phred64-quals",		no_argument,		0,	OPT_PHRED64_QUALS},
{"no-relabel",          no_argument,		0,	OPT_NO_RELABEL},
{0, 0, 0, 0} // terminator
};

int parse_options(int argc, char** argv, void (*print_usage)())
{
    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:    
				break;
			case OPT_FASTA:
				reads_format = FASTA;
				break;
			case OPT_FASTQ:
				reads_format = FASTQ;
				break;
			case OPT_VERBOSE:
				verbose = true;
				break;
			case OPT_OUTPUT_DIR:
				output_dir = optarg;
				break;
			case OPT_MAX_MULTIHITS:
				max_multihits = parseInt(1, "--max-multihits arg must be at least 1", print_usage);
				break;
			case OPT_SEGMENT_LENGTH:
				segment_length = parseInt(4, "--segment-length arg must be at least 4", print_usage);
				break;
			case OPT_SEGMENT_MISMATCHES:
				segment_mismatches = parseInt(0, "--segment-mismatches arg must be at least 0", print_usage);
				break;
			case OPT_SOLEXA_QUALS:
				solexa_quals = true;
				break;
			case OPT_PHRED64_QUALS:
				phred64_quals = true;
				break;
            case OPT_NO_RELABEL:
				fastq_db = false;
				break;
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
    return 0;
}

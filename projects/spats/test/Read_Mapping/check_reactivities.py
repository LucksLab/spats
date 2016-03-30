#!/usr/bin/env python

"""
Summarizes resulting reactivities.out of spats based on expected results. Part of Read_Mapping tests.
Usage: python check_reactivities.py --input <> --output <> --linker <> --sequence <>
Options: 
 --input	reactivities.out file from spats
 --output	File to output result summary
 --linker	Linker sequence
 --sequence	Sequence of RNA target
Version: 0.1
Date: March 29, 2016
Author: Angela M Yu
Copyright (c) 2016 Lucks Laboratory - all rights reserved.
"""

import getopt
import sys
from itertools import repeat

def getopts(short_arg_string, long_arg_list):
    """
    Returns a dictionary of command line arguments as defined by short_arg_string and long_arg_list
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:],short_arg_string,long_arg_list)
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(2)
    return dict(opts)

opts = getopts("", ["input=", "output=", "linker=", "sequence="])
spats_reactivities_out = opts["--input"]
output = opts["--output"]
linker_len = len(opts["--linker"])
sequence_len = len(opts["--sequence"])

# Read in reactivities.out
reads = []
with open(spats_reactivities_out, "r") as f:
    header = f.readline()
    for line in f:
        fields = line.split("\t")
        reads += [int(fields[4]), int(fields[5])] # treated_reads is [4], untreated_reads is [5]

# Build expected reads
case_exp = zip(range(1,sequence_len+2), range(1, sequence_len+2)[::-1]) #creates pairs [(1,N), (2,N-1), ..., (N+1,1)]
case_exp = [ele for i in case_exp for ele in i] + [0]*(2*linker_len-4)  #add 0s where linker does not map

# Calculate summary
correct = sum([1 if a[0] == a[1] else 0 for a in zip(case_exp, reads)])
incorrect = len(reads) - correct
expected_read_lines = 2 * (sequence_len + linker_len - 1)

if correct == len(reads) and len(reads) == expected_read_lines:
    result = "OK - %s read positions out of %s expected, %s correct, %s incorrect\n"%(len(reads), exp_read_lines, correct, incorrect)
else:
    result = "FAILED - %s read positions out of %s expected, %s correct, %s incorrect\n"%(len(reads), exp_read_lines, correct, incorrect)

with open(output, "a") as f:
    f.write(result)

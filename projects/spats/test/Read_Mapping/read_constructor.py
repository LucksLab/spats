#!/usr/bin/env python
"""
read_constructor.py builds fastq reads from input RNA sequence fragments to help simulate read mapping by spats.

This script was written in July 2015 by JBL and edited in March 2016 by AMY.
Copyright (c) 2016 Lucks Laboratory - all rights reserved.
"""

import sys
import os
import getopt
import shutil
from datetime import datetime, date, time
import string
import fileinput
import random
from itertools import product

current_dir = os.getcwd()+'/'

help_message = '''
read_constructor.py builds fastq reads from input RNA sequence fragments to help simulate read mapping by spats.

Usage:
   read_constructor.py [options] [-s OR -f]

Options
-h, --help                  opens help message
-v, --version               displays version number
-o, --output                Output file name. Default R1.fq, R2.fq
-l, --length                Read length of output reads. Default is 35.
-s, --sequence              DNA sequence to use to construct reads representing a stop at each position.
-f, --file                  File containing DNA sequences to use.
-a, --adapter <sequence>    Adapter sequence (5'->3') to add at the 3' end of the RNA w/ barcodes (if applicable)
                            as it appears on the RT primer (after handle) (ex CACTCGGGCACCAAGGA)
--permutation               Select which permutation of set of 4 expected reactivity mappings to use. 0-15, default 1
'''

def get_version():
    return "0.0.2"


class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg


class Params:
    
    def __init__(self):
        pass
    
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvo:l:s:f:a:",
                                       ["version",
                                        "help",
                                        "output=",
                                        "length=",
                                        "sequence=",
                                        "file=",
                                        "adapter=",
                                        "permutation="])
        
        #Add linker option
    

        except getopt.error, msg:
            raise Usage(msg)
        
        out_1 = 'R1.fq'
        out_2 = 'R2.fq'
        input_file = None
        sequence = None
        read_length = 35
        adapter = ''
        permutation = 1

        for option, value in opts:
            if option in ("-v", "--version"):
                print "read_constructor.py v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o","--output"):
                out_1 = value+'_R1.fq'
                out_2 = value+'_R2.fq'
            if option in ("-l","--length"):
                read_length = int(value)
            if option in ("-s","--sequence"):
                sequence = value
            if option in ("-f","--file"):
                input_file = value
            if option in ("-a","--adapter"):
                adapter = value
            if option == "--permutation":
                permutation = int(value)
            
        if (sequence is None) and (input_file is None):
            raise Usage('At least one of -s and -f must be specified'+help_message)
        
        return args,out_1,out_2,read_length,input_file,sequence,adapter,permutation
    
    def check(self):
        pass


def fastq_format(r_1,r_2,read_name=""):
    """format a read pair in fastq format"""
    r_1_string = ''
    r_2_string = ''
    qual_string = ''.join('I' for n in range(len(r_1)))

    if read_name is "":
        #choose random id
        rand_id = '@'+''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(6)])
    else:
        rand_id = '@'+read_name

    r_1_string = rand_id+'\n'+r_1+'\n+\n'+qual_string
    r_2_string = rand_id+'\n'+r_2+'\n+\n'+qual_string
    
    return [r_1_string,r_2_string]


def reverse_complement(s):
    #This section was taken from Cole's code
    nuc_table = { 'A' : 'T',
                  'T' : 'A',
                  'C' : 'G',
                  'G' : 'C',
                  'a' : 't',
                  't' : 'a',
                  'c' : 'g',
                  'g' : 'c'  }
    sl = list(s)
    try:
        rsl = [nuc_table[x] for x in sl]
    except KeyError, k:
        print >> sys.stderr, "Error: adapter sequences must contain only A,C,G,T"
        exit(1)
    rsl.reverse()
    return ''.join(rsl)


def write_read_files(f_out_1, f_out_2, read_1_string, read_2_string):
    f_out_1.write(read_1_string)
    f_out_1.write('\n')
    f_out_2.write(read_2_string)
    f_out_2.write('\n')


def main(argv=None,):
    
    #Definitions of A_T and A_B based on SHAPE-Seq 2.0 (Loughrey, Watters NAR Figure S4)
    adapter_T = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    adapter_B = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

    all_permutations = [seq for seq in product([0,1], repeat=4)]
    
    params = Params()
        
    try:
        if argv is None:
            argv = sys.argv
            args,out_1,out_2,read_length,input_file,sequence,adapter,permutation = params.parse_options(argv)
            params.check()
        input_fragments = [] 
        
        #Load sequences to build with
        if input_file is not None:
            for line in fileinput.input(input_file):
                input_fragments.append(line.rstrip('\n'))
        
        if sequence is not None:
            for i in range(len(sequence)):
                input_fragments.append(sequence[i:])
        
        input_fragments.append('') #always include empty fragment
        
        #Open output files
        f_out_1 = file(out_1,'w')
        f_out_2 = file(out_2,'w')
        
        #Loop through input fragments

        read_i = -1
        read_ii = 0
        create_reads_perm = all_permutations[permutation]
        for seq in input_fragments:
            read = adapter+reverse_complement(seq)
            read_i += 1

            #Append +/- barcode
            plus_read = 'GGGC'+read
            minus_read = 'CCCG'+read
            
            for lib_read in (plus_read,minus_read):
                read_ii += 1
                #Construct Read_1 and Read_2 based on length
                temp_read = lib_read+adapter_T
                temp_read = lib_read+adapter_B
                read_1 = temp_read[:read_length]
                
                temp_read = reverse_complement(adapter_B+lib_read)
                temp_read = reverse_complement(adapter_T+lib_read)
                read_2 = temp_read[:read_length]

                #Format read into Fastq
                read_1_string,read_2_string = fastq_format(read_1,read_2,str(read_i)+"_"+str((read_ii+1)%2))
                
                #Print to out_1 and out_2 if defined by permutation
                if read_ii % 4 == 1 and create_reads_perm[0] == 1:
                    write_read_files(f_out_1, f_out_2, read_1_string, read_2_string)
                if read_ii % 4 == 2 and create_reads_perm[1] == 1:
                    write_read_files(f_out_1, f_out_2, read_1_string, read_2_string)
                if read_ii % 4 == 3 and create_reads_perm[2] == 1:
                    write_read_files(f_out_1, f_out_2, read_1_string, read_2_string)
                if read_ii % 4 == 0 and create_reads_perm[3] == 1:
                    write_read_files(f_out_1, f_out_2, read_1_string, read_2_string)

        f_out_1.close()
        f_out_2.close()


    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
"""
read_constructor.py builds fastq reads from input RNA sequence fragments to help simulate read mapping by spats.

This script was written in July 2015 by J. B. Lucks and edited in March 2016 by Angela M.Yu.
Copyright (c) 2016 Lucks Laboratory - all rights reserved.
"""

import sys
import os
import getopt
import string
import fileinput
import random

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
-s, --sequence              DNA sequence to use to construct reads for each stop position. N Plus reads are generated at
                              each position where N increments from 5' to 3' (i.e. there is one full length plus read, two
                              N-1 plus reads, etc.) L-N minus Minus reads are generated at each position where L is the length
                              of the sequence (i.e. there are L full length minus reads, L-1 minus reads at position N-1, etc.)
-f, --file                  File containing specific DNA sequences to use. Only these sequences will be mapped. Use --sequence to automatically
                              construct a comprehensive test case.
-a, --linker <sequence>     Linker sequence (5'->3') to add at the 3' end of the RNA w/ barcodes (if applicable)
                            as it appears on the RT primer (after handle) (ex CACTCGGGCACCAAGGA)
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
                                        "linker="])
        

        except getopt.error, msg:
            raise Usage(msg)
        
        out_1 = 'R1.fq'
        out_2 = 'R2.fq'
        input_file = None
        sequence = None
        read_length = 35
        linker = ''

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
            if option in ("-a","--linker"):
                linker = value
            
        if (sequence is None) and (input_file is None):
            raise Usage('At least one of -s and -f must be specified'+help_message)
        
        return args,out_1,out_2,read_length,input_file,sequence,linker
    
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
        print >> sys.stderr, "Error: linker sequences must contain only A,C,G,T"
        exit(1)
    rsl.reverse()
    return ''.join(rsl)


def write_read_files(f_out_1, f_out_2, read_1_string, read_2_string):
    f_out_1.write(read_1_string)
    f_out_1.write('\n')
    f_out_2.write(read_2_string)
    f_out_2.write('\n')


def main(argv=None,):
    #Definitions of A_T and A_B based on SHAPE-Seq 2.0 (Loughrey, Watters NAR Figure S4, doi: 10.1093/nar/gku909)
    Adapter_T = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    Adapter_B = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

    params = Params()
        
    try:
        if argv is None:
            argv = sys.argv
            args,out_1,out_2,read_length,input_file,sequence,linker = params.parse_options(argv)
            params.check()
        input_fragments = [] 
        
        #Load sequences to build with
        if input_file is not None:
            for line in fileinput.input(input_file):
                input_fragments.append(line.rstrip('\n'))
        
        if sequence is not None:
            for i in range(len(sequence)):
                #build sequence fragments with decreasing length stopping at every position
                input_fragments.append(sequence[i:]) 
        
        input_fragments.append('') #always include empty fragment
        
        #Open output files
        f_out_1 = file(out_1,'w')
        f_out_2 = file(out_2,'w')
        
        #Loop through input fragments

        read_i = -1
        read_ii = 0
        for seq in input_fragments:
            if True: #JBL - need this?
                # Constructing reads. See Loughrey, Watters NAR Figure S4, doi: 10.1093/nar/gku909
                read = linker+reverse_complement(seq)
                read_i += 1

                #Prepend +/- handle index
                plus_read = 'GGGC'+read # Plus is RRRY
                minus_read = 'CCCG'+read # Minus is YYYR
            
                for lib_read in (plus_read,minus_read):
                    read_ii += 1
                    #Construct Read_1 and Read_2 based on length
                    temp_read = lib_read+Adapter_B #temp_read contains full top strand
                    read_1 = temp_read[:read_length] #cut down by read_length to get read_1
                
                    temp_read = reverse_complement(Adapter_T+lib_read) #temp_read contains full bottow strand, revcomp based on sequencing process
                    read_2 = temp_read[:read_length] #cut down by read_length to get read_2

                    #Format read into Fastq
                    read_1_string,read_2_string = fastq_format(read_1,read_2,str(read_i)+"_"+str((read_ii+1)%2))

                    if read_ii % 2 == 0:
                        #For every minus read, repeat this mate pair as many times as the length of the sequence
                        # Since fragments are in a list in descending length, this creates a pattern of descending 
                        # number of reads at each postion. i.e. Full length minus reads equal the sequence length,
                        # full length -1 minus reads equal sequence length -1 etc.
                        for _ in range(len(seq)+1):
                            write_read_files(f_out_1, f_out_2, read_1_string, read_2_string)
                    else:
                        #For every plus read, repeat this mate pair in ascending length. i.e. 1 full length plus read,
                        # 2 full length -1 plus reads, etc.
                        for _ in range(len(sequence)-len(seq)+1):
                            write_read_files(f_out_1, f_out_2, read_1_string, read_2_string)

        f_out_1.close()
        f_out_2.close()


    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())

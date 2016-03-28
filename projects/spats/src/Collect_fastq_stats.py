#!/usr/bin/env python

#Search_fastq_R1.py takes the output of ReadsAnalyzer and searches for Spats targets
#to determine the number of reads obtained per target before alignment with Spats

#This script was written in August 2014 by Kyle Watters
#Copyright (c) 2014 Kyle Watters. All rights reserved.

import sys
import os
import getopt
from collections import OrderedDict

current_dir = os.getcwd()+'/'

help_message = '''
Search_fastq_R1.py uses a targets file for spats to search the R1 output of ReadsAnalyzer
to determine how many reads can aligned to each target.  Can be used for %Aligned calculations

Usage: Search_fastq_R1.py [options] <ReadsAnalyzer Output> <targets filename>

Options
-h, --help                  opens help message
-v, --version               displays version number
-o, --output                changes the output text file name (default: <targets_filename>_Counts-per-target.txt)
'''

def get_version():
    return "0.0.1"

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    
    def __init__(self):
        pass
    
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvo:",
                                       ["version",
                                        "help",
                                        "output="])
        
        except getopt.error, msg:
            raise Usage(msg)
        
        for option, value in opts:
            if option in ("-v", "--version"):
                print "Search_fastq_R1.py v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o","--output"):
                output_filename = value
      
        if len(args) != 2:
            raise Usage(help_message)    
    
        output_filename = args[0].split(".")[0] + "_Counts-per-target.txt"
 
        for option, value in opts:
            if option in ("-o","--output"):
                output_filename = value
        
        return args,output_filename
    
    def check(self):
        pass

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

def discover_read_len(input_file):
        #Detect the read length in the collapsed file
    with open(input_file,"r") as input_read:
        first_line = input_read.readline()
        read_len = len(first_line.split("\t")[1].strip())
    return read_len        
            
def main(argv=None,):
    
    params = Params()
        
    try:
        if argv is None:
            argv = sys.argv
            args,output_filename = params.parse_options(argv)
            params.check()
        
        input_file = args[0]
        targets_file = args[1]
            
        read_len = discover_read_len(input_file)     
        #Load the targets file and determine what to search for from the 3' end of each
        targets = file(targets_file,"r")
        lines = []
        for line in targets:
            read_line = line.strip()
            if read_line != "":
                if read_line[0] == ">":
                    read_line = read_line + " "
                lines.append(read_line)
        split_lines = ''.join(lines).split(">")[1:]
        all_targets = OrderedDict()
        for target in split_lines:
            #Takes the target name and sequence pairs, and splits them into a dictionary, using the rev-comp of the 3' end of the target (same length as reads)
            all_targets[target.split(" ")[0]] = reverse_complement(target.split(" ")[1].strip()[-read_len:])
        
        #Initialize results txt file
        with open(output_filename, "w") as output:
            firstline = "\t".join(['Target', 'Counts', 'Sequence in Read 1\n'])
            output.write(firstline)
            pass
        
        #Search for each target sequentially, if found, stop searching, and add to the results txt file
        num_lines = sum(1 for line in open(input_file))
        
        iteration = 0
        for key in all_targets:
            inputfile = file(input_file, "r")
            for line in inputfile:
                iteration += 1
                temp_counts = line.strip().split("\t")[0].split("-")[1]
                temp_seq = line.strip().split("\t")[1]
                if temp_seq == all_targets[key]:
                    with open(output_filename, "a") as outputfile:
                        outputfile.write("\t".join([key, temp_counts, temp_seq])+"\n")
                    break
                if iteration == num_lines:
                    with open(output_filename, "a") as outputfile:
                        outputfile.write("\t".join([key, '0', temp_seq])+"\n")
            iteration = 0
        #Search line by line of the collapsed file, stop on first hit and write to file with target name
        
                
                
                
                
                
        #Reverse complement of adapter is added to RNA sequence, so rev-comp:
       # adapter = reverse_complement(adapter)
           
       # with open(output,"w") as output_text:
       #     for pos in range(min_len, max_len+1):
       #         if adapter_name != '':
       #             name_label = ">" + RNA_name + "_" + adapter_name + "_" + str(pos) + "nt" +  "\n"
       #         else:
       #             name_label = ">" + RNA_name + "_" + str(pos) + "nt" +  "\n"
       #         RNA_seq_processed = RNA_seq[:pos].upper() + adapter.upper() + "\n\n"
       #         output_text.write(name_label) 
       #         output_text.write(RNA_seq_processed)   
       #         if pos == max_len:
       #              print >> sys.stderr, ""
       #              print >> sys.stderr, "Longest target created: "
       #              print >> sys.stderr, "--------------------------"
       #              print >> sys.stderr, name_label,RNA_seq_processed
           
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())
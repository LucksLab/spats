#!/usr/bin/env python

#Cotrans_targets builds a targets file for spats automatically for all lengths of an RNA
#using all of the intermediate lengths - that is, for a 100 nt RNA, it will automatically make 
#the targets file for all the intermediate lengths from 5' -> 3' with the appropriate 3' adapter

#This script was written in February 2014 by Kyle Watters
#Copyright (c) 2014 Kyle Watters. All rights reserved.

import sys
import os
import getopt
import shutil
from datetime import datetime, date, time
import string

current_dir = os.getcwd()+'/'

help_message = '''
Cotrans_targets.py generates a targets file for spats when analyzing intermediate lengths 
of an RNA with an attached adapter for RT priming.

For an RNA of N nt, produces N-y (y=user chosen length) intermediate lengths in a
targets file with the supplied adapter sequence at the 3' end.

Usage:
   Cotrans_targets.py [options] <RNA name> <RNA sequence> <targets filename(.fa recommended)>

Options
-h, --help                  opens help message
-v, --version               displays version number
-a, --adapter <sequence>    Adapter sequence (5'->3') to add at the 3' end of the RNA w/ barcodes (if applicable)
                            as it appears on the RT primer (after handle) (ex for IDT 2: gtccttggtgcccgagtg; for IDT2_mod: gtccttggtgcccgagtcag)
-n, --adapter-name <string> Name for adapter to include in names in targets file (ex: IDT2)                                                                      
-i, --min-len <N>           Minimum length of RNA from 5' end to make target with (default = 20)
-m, --max-len <N>           Maximum length of RNA from 5' end to make target with (default is entire RNA)
-x, --entire-only           Only produce the target for the max length
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
            opts, args = getopt.getopt(argv[1:], "hva:i:m:xn:",
                                       ["version",
                                        "help",
                                        "adapter=",
                                        "min-len=",
                                        "max-len=",
                                        "adapter-name=",
                                        "entire-only"])
        
        except getopt.error, msg:
            raise Usage(msg)
        
        adapter = ''
        adapter_name = ''
        min_len = 20
        max_len = 0
        
        for option, value in opts:
            if option in ("-v", "--version"):
                print "Cotrans_targets.py v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-a","--adapter"):
                adapter = value
            if option in ("-n","--adapter-name"):
                adapter_name = value
            if option in ("-i","--min-len"):
                min_len = int(value)
            if option in ("-m","--max-len"):
                max_len = int(value)
            
        if len(args) != 3:
            raise Usage(help_message)    
        
        #Set the maximum length to the RNA sequence length if no maximum length provided
        if max_len == 0: 
            max_len = len(args[1])
    
        for option, value in opts:    
            if option in ("-x","--entire-only"):
                min_len = max_len
        
        return args,adapter,min_len,max_len,adapter_name
    
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

def main(argv=None,):
    
    params = Params()
        
    try:
        if argv is None:
            argv = sys.argv
            args,adapter,min_len,max_len,adapter_name = params.parse_options(argv)
            params.check()
        
        RNA_name = args[0]
        RNA_seq = args[1]
        output = args[2]
                
        #Reverse complement of adapter is added to RNA sequence, so rev-comp:
        adapter = reverse_complement(adapter)
           
        with open(output,"w") as output_text:
            for pos in range(min_len, max_len+1):
                if adapter_name != '':
                    name_label = ">" + RNA_name + "_" + adapter_name + "_" + str(pos) + "nt" +  "\n"
                else:
                    name_label = ">" + RNA_name + "_" + str(pos) + "nt" +  "\n"
                RNA_seq_processed = RNA_seq[:pos].upper() + adapter.upper() + "\n\n"
                output_text.write(name_label) 
                output_text.write(RNA_seq_processed)   
                if pos == max_len:
                     print >> sys.stderr, ""
                     print >> sys.stderr, "Longest target created: "
                     print >> sys.stderr, "--------------------------"
                     print >> sys.stderr, name_label,RNA_seq_processed
           
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    sys.exit(main())

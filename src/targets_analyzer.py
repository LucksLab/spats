#!/usr/bin/env python

"""
targets_analyzer.py

Created by Ray Zhuang, Kyle E. Watters and Julius B. Lucks 2013-2014.
Copyright (c) 2014 Ray Zhuang, Kyle E. Watters and Julius B. Lucks. All rights reserved.
"""

import sys
import os
import getopt

from sys import stdout

present_dir = os.getcwd()+"/"

help_message = '''
targets_analyzer takes a targets file for alignment in spats and determines how many nucleotides
at the 3' end are required for unique alignment as well as to search for A_adapter_b similarity to set 
a proper trim-match value for adapter clipping.

Usage: 
    targets_analyzer <targets.fa> <A_b_sequence>

Options:
-h,--help           brings up help
-v,--version        displays version number
'''

class Usage(Exception):
   def __init__(self,msg):
       self.msg = msg

def get_version():
    return "0.0.1"

class Params:
    def __init__(self):
       pass

    def parse_options(self, argv):
        
        try:
            opts, args = getopt.getopt(argv[1:],"hv",["help","version"])            

        except getopt.error, msg:
            raise Usage(msg)
 
        for option, arg in opts:
            if option in ("-v", "--version"):
                print "targets_analysis v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
                
        if len(args) != 1:
            raise Usage(help_message)
    
        return args
        
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

def readFASTA(input_file,A_b_sequence):
 
    reader = open(input_file)
    thisFile = reader.readlines()    
    
	# Lower bound on trim-match
	# Will search for better minCount that will ultimately be used for trim-match
    minCount = 5

	#The number of nucleotides on the 3' end needed to provide complete uniqueness among targets
    uniqueNuc = 0

    #Create a local array to hold everything
    tempArray = []

    #Create array that will contain all unique sequences seen so far
    uniqArray = []

    #Here, we just get the smaller set of sequences, without the titles
    for arrayEl in thisFile:
        realEl = arrayEl.strip() #last array element is the newline, so we'll skip it. Python's :blah will not include the blah index
        #KEW ^replaced with strip to remove any trailing newlines, etc.(in cases spaces, etc. are left too)
        if (arrayEl != "\n") & (arrayEl[0] != ">"):
            tempArray.append(realEl)
            #print realEl[(len(realEl)-1):] #This is the last nucleotide
    
    seqArray = tempArray
    shortest = False
    
    #Scan through all possible sequences on the 3' end, getting successively larger, until find uniqueNuc
	## Starting search at 1 nt from the 3' end
	## If uniqArray contains fewer elements than number of seqs (len(tempArray)), then not unique
	## Try a bigger sequence chunk
    while (len(uniqArray) != len(tempArray)):
        for el in tempArray:
            if (not uniqArray.__contains__(el[(len(el)-1-uniqueNuc):])): #If we haven't seen it, then append it
                uniqArray.append(el[(len(el)-1-uniqueNuc):])
        #print uniqArray
        if (len(uniqArray) != len(tempArray)):
            uniqArray = []
            uniqueNuc = uniqueNuc + 1
    
	#A_b_sequence comes in reads as a revcomp
    #First reverse complement of the file, 
    #then add AGAT at the front (to mimic it coming up randomly in the handle sequence)
	#this is a conservative case of some of the 3' end matching A_b_sequence and the handle completing the match
    for el in range(0,len(seqArray)):
        # seqArray[el] = "AGAT" + reverse_complement(seqArray[el]) #JBL - should change to adding 1st 4 of A_b_sequence[0:3]
		seqArray[el] = A_b_sequence[0:3] + reverse_complement(seqArray[el]) 
    
	# JBL REFACTOR - this code should be replaced with re.find implementation
    while ((shortest == False) & (minCount < 35)):
        uniqArray = [] #KEW - needed to flush string to 'reset' length search with each loop
        A_b_match = A_b_sequence[0:minCount]  #A_b_sequence fragment we're trying to match
        #Now check if we already have that sequence
        for seq in seqArray:
            for char in range(0,len(seq)-minCount+1):
				# Scanning through sequence chuncks of length minCount to search for A_b_match
                if ((seq[char:(char+minCount)] == A_b_match)):
                    uniqArray.append(seq)
        #print(uniqArray)
        if(len(uniqArray) != 0): #Should be empty when A_b cannot be found at the given length
            shortest = False
            minCount += 1
        else:
            shortest = True

    return uniqueNuc, minCount

def main(argv=None):

    params = Params()
    
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()     

        inputfile = args[0]
        
        try:
            A_b_sequence = args[1]
        except IndexError:
			A_b_sequence = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG" #default
   
        #stdout.write("\rRunning the Parser for Script........")
        #stdout.flush()       
        
	#Get number of nucleotides to search for adapter	
	unique_num,A_b_search_len = readFASTA(inputfile,A_b_sequence)
        #print("done\n")
        print("Number of nucleotides from 3' end needed to be unique: " + str(unique_num))     
        print("AGAT search length is: " + str(A_b_search_len))              

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2     

if __name__ == "__main__":
    sys.exit(main())


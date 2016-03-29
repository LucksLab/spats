#!/usr/bin/env python
"""
spats.py

Created by Kyle Watters and Julius Lucks 2014
Copyright (c) 2016 Lucks Laboratory. All rights reserved.

Revision by Kyle Watters on 2014-9-20.
change notes:
    -switched fastx_toolkit to cutadapt script (v1.5 at time of change)
"""

import analyze_spats_targets
import fastq_revcomp
import sys
import os
import getopt
import shutil
from datetime import datetime, date, time
import string
import random
import re

current_dir = os.getcwd()+'/'

name = "adapter_trimmer"

help_message = '''
adapter_trimmer trims adapter sequences from fastq files in preparatino for processing by spats.

takes fastq files, searches for short sequences that are reverse complimentary,
trims off adapter sequences, and prepares reads for processing directly with spats

Usage:
   adapter_trimmer [options] <R1_seq.fastq> <R2_seq.fastq> <targets.fa>

Options:
-h, --help                  opens help message
-v, --version               displays version number
-e <N>, --error-rate <N>    allowed error rate when clipping adapters (0.1 is 1/10 bases, etc.)
-q <N>, --quality-min <N>   minimum quality score to consider, 
-p <N>, --num-threads <N>   Number of processors to use (default = 1)
--allow-mismatches          Allow downstream mismatches for Spats (keep sequences that don't perfectly align)
--set-quality <N>           sets all read qualities to the given quality where ascii is (N+33)               
--trim-match <N>            number of nucleotides of adapter to search for with clipper (default = 6) *Adjust trim length this way.
--read-len <N>              Specify a length of reads (ex. 35 for 2x35 bp paired end reads) (autodetected)
--min-read-len <N>          Number of nt on the 3' end to leave after (w/ adapter) (default = 6)
--A-b-sequence <string>     A_adapter_b sequence (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
--A-t-sequence <string>     A_adapter_t sequence (default: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)
--max-handle-len <N>        Number of nucleotides in the longest (+/-) handle (default = 4 for RRRY/YYYR)
-o, --output <string>       If running many simultaneous runs in the same folder, use this option to name the treated reads
                            folder to avoid overwriting during processing

JBL - if trim-match or min-read-len are not set, they will automatically be determined

The algorithm employed is as follows:
 
1.) Filter reads into two sets: Case_I = do not have any adapters at all - i.e. long insert sizes, Case_II = have some adapter.
2.) For Case II, trim all possible with cutadapt (Case_II_clipped).
3.) For reads that aren't trimmed from Case II (Case_II_unclipped):
3.1) First filter these reads out with cutadapt (i.e. collect just Case_II_unclipped)
3.2) Execute trim_search algorithm
4.) Combine all reads into one file to send to spats

'''

def get_version():
    return "0.2.0"

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    
    def __init__(self):
        pass
    
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvo:e:q:p:",
                                       ["version",
                                        "help",
                                        "trim-match=",
                                        "read-len=",
                                        "A-b-sequence=",
                                        "A-t-sequence=",
                                        "min-read-len=",
                                        "output=",
                                        "max-handle-len=",
                                        "error-rate=",
                                        "quality-min=",
                                        "set-quality=",
                                        "num-threads=",
                                        "allow-mismatches"])
        
        except getopt.error, msg:
            raise Usage(msg)
        
        for option, value in opts:
            if option in ("-v", "--version"):
                print "%s v%s" % (name,get_version())
                exit(0)
            elif option in ("-h", "--help"):
                raise Usage(help_message)
        
        
        trim_match = 0          #default length of adapter to search
        trim_auto = True
        read_len = -1         
        A_b_sequence = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        A_t_sequence = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        final_dir = current_dir + "combined_reads_" + id_generator(10) + "/"
        max_handle_len = 4
        min_read_auto = True
        min_read_len = 6
        error_rate = 0
        quality_min = 0
        quality_change = None
        mismatches = False
        threads = 1
        for option, value in opts:
            if option in ("--trim-match"):
                trim_auto = False
                trim_match = int(value)
            elif option in ("--read-len"):
                read_len = int(value)
            elif option in ("--A-b-sequence"):
                A_b_sequence = value
            elif option in ("--A-t-sequence"):
                A_t_sequence = value
            elif option in ("--min-read-len"):
                min_read_len = int(value)
                min_read_auto = False
            elif option in ("-o","--output"):
                final_dir = value
            elif option in ("-e","--error-rate"):
                error_rate = value
            elif option in ("-q","--quality-min"):
                quality_min = value
            elif option in ("-p", "--num-threads"):
                threads = int(value)
            elif option in ("--allow-mismatches"):
                mismatches = True
            elif option in ("--max-handle-len"):
                max_handle_len = int(value)
            elif option in ("--set-quality"):
                if int(value) >= 0 and int(value) <= 93:
                    quality_change = int(value)
                    print >> sys.stderr, "Manually changing all phred scores to '{0}'.".format(chr(quality_change))
                else:
                    print >> sys.stderr, "Quality to set {0} is out of range".format(value)   
        # JBL - A_t shows up as a revcomp in read2 because of library format
        A_t_sequence = reverse_complement(A_t_sequence)
        if read_len == -1:
            read_len = find_read_len(args[0])
        
        if len(args) != 3:
            raise Usage(help_message)
        return [args,trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,
                    max_handle_len,error_rate,quality_min,quality_change,threads,mismatches]
    
    def check(self):
        pass

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    '''From http://stackoverflow.com/questions/2257441/python-random-string-generation-with-upper-case-letters-and-digits'''
    return ''.join(random.choice(chars) for x in range(size))

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)

def find_read_len(fastq):
    #Check that the file is fasta
    extension = fastq.split(".")[-1]
    if extension != "fq" and extension != "fastq":
        print >> sys.stderr, "File is not a fasta extension (fq,fastq)"
    
    #Check file exists
    elif not os.path.isfile(fastq):
        print >> sys.stderr, "No such file as {0}".format(fastq)
        
    else:
        #Detect the read length in the collapsed file
        fp = open(fastq)
        for i, line in enumerate(fp):
            if i == 1:
                read_len = len(line.strip())
            elif i > 1:
                break
        fp.close()
        #print("Read length is {0} nt.".format(read_len))
        
        return read_len

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
 
def trim_min_read_calculator(trim_auto,min_len_auto,input_targets,A_b_sequence,max_handle_len,trim_match,min_read_len,read_len):
    
    #Find the minimal lengths of A_b to search and number of nucleotides to trim from 3' end
    uniqueness = analyze_spats_targets.analyze_unique(input_targets,A_b_sequence)
    A_b_search_len,error_rate = analyze_spats_targets.analyze_clip_min(input_targets,A_b_sequence)
    
    if min_len_auto == True:
        # Add handle length to uniqueness length to find min_read_length considering
        if uniqueness < 6:
            min_read_len = max_handle_len + uniqueness   #
        else:
            min_read_len = 10  #still leaves at least 6 nt of shared sequence to align to the target
        print >> sys.stderr, "[%s] Automatically dropping %s nt from end of RNAs (from 3' end analysis)" % (right_now(),uniqueness)
    else:
        print >> sys.stderr, "[%s] Manually dropping %s nt from end of reads" % (right_now(),min_read_len)
    
    if trim_auto == True:
        if A_b_search_len >= read_len:
            trim_match = read_len - min_read_len 
        else:
            trim_match = A_b_search_len
        print >> sys.stderr, "[%s] Automatically set trim-match to %s based on targets file" % (right_now(),trim_match)
    else:
         print >> sys.stderr, "[%s] Trim-match manually set to %s" % (right_now(),trim_match)
    
    return trim_match,min_read_len,error_rate

def quality_mod(inputfile,read_len,quality_change):
    #Modifies all of the quality reads to to the user's request, usually so that phred scores are ignored
    
    qualmodfile = os.path.splitext(inputfile)[0] + '_qualmod.fq'
    qualstring = chr(quality_change) * read_len
    qualcommand = "awk \'BEGIN{s=0} {s += 1; if (s % 4 == 0) print \"" + qualstring + "\"; else print }\'"
    #this is needed to except the curly brackets in the below command
    os.system(qualcommand +" {0} > {1}".format(inputfile, qualmodfile));
    #print("done")
    return qualmodfile

def clip_search(base,read_len,max_handle_len,output_dir,targets_dir,input_R1,input_R2,threads):
    """Search for revcomp reads using bowtie after trimming to length base
    
    Assuming input_R1 and input_R2 are filenames of files that have been quality_mod()'d
    
    """
    print >> sys.stderr, "[%s] Finding reads that are %s nt long" % (right_now(), (base-max_handle_len)) #JBLQ - what does this mean?
    
    #Make a temporary folder for intermediate files to delete after each trim
    search_dir = current_dir + "search_" + id_generator(10) + "/"
    try:
        os.mkdir(search_dir)
    except:
        exit # should terminate here since need this search_dir to proceed
    
    #Implement Adapter Clipping Step 3: Find reads that the clipper will miss
    ## Trim base-by-base and look for revcomp sequences. 3 cases:
        ## 1) Reads contain no adapter, thus will not revcomp and will not keep these mate pairs.
        ## 2) Reads contain some adapter, will not revcomp and will not keep these mate pairs.
        ## 3) Adapter has been manually trimmed out, reads will revcomp, keep this pair.
    ## Accumulate kept reads in files - will combine all reads at the end of analysis.
    
    #Block of filename constructions for outputfiles
    if base == read_len - 1:
        output_R1 = os.path.splitext(input_R1)[0] + "_" + str(base) + ".fq"
        output_R2 = os.path.splitext(input_R2)[0] + "_" + str(base) + ".fq"
        output_revcomp_R1 = os.path.splitext(input_R1)[0] + "_" + str(base) + "_rc" + ".fq"
    else:
        output_R1 =  output_dir + "Case_II_1_unclipped_{0}.fq".format(base)
        output_R2 =  output_dir + "Case_II_2_unclipped_{0}.fq".format(base)
        output_revcomp_R1 = output_dir + "Case_II_2_unclipped_{0}_rc.fq".format(base)
    output_bowtie = output_dir + "bowtie_output_" + str(base) + ".sam"
    bowtie_results = output_dir + "bowtie_results{0}.fq".format(base)
    bowtie_results_1 = output_dir + "bowtie_results{0}_1.fq".format(base)
    bowtie_results_1_rc = output_dir + "bowtie_results{0}_1_rc.fq".format(base)
    bowtie_unmatched = output_dir + "bowtie_results{0}_unmatched.fq".format(base)
    bowtie_unmatched_1 = output_dir + "bowtie_results{0}_unmatched_1.fq".format(base)
    bowtie_unmatched_1_rc = output_dir + "bowtie_results{0}_unmatched_1_rc.fq".format(base)
    bowtie_unmatched_2 = output_dir + "bowtie_results{0}_unmatched_2.fq".format(base)
    
	#Trim down each file to base
	## Construct cutadapt trimming command to remove nucleotides unconditionally
	## -u N = Number of nucleotides to remove from the end of a read (when N < 0)
    trim_num = read_len - base
    trim_num = 1 #JBL FLAG
    os.system("cutadapt --quiet -u -{0} -o {1} {2}".format(trim_num,output_R1,input_R1))
    os.system("cutadapt --quiet -u -{0} -o {1} {2}".format(trim_num,output_R2,input_R2))
    
    #Use a loaded function to do the reverse comps
    fastq_revcomp.rev_comp(output_R1,output_revcomp_R1)
    #os.system("fastx_reverse_complement -i {0} -o {1}".format(output_R1,output_revcomp_R1))
	
	# Construct bowtie command
	## -q :: fastq input files
	## -v 0 :: report hits with no mismatches - i.e. only perfect revcomps matched
	## --ff :: aligning forward, forward (since revcomped a read) (see diagram - Read2 naturally aligns in forward with targets)
	## --norc :: do not align to reverse-complement reference strand (forcing reads to align in one direction)
	## -3 4 :: remove 4 nts from 3' end (this gets rid of YYYR and RRRY from revcomped read 1, or end of read 2, otherwise would not align)
	## -X 2000 :: setting huge upper limit on insert size (i.e. no limit)
	## -m 1 :: unique alignments
	## --al {4}:: This will generate two files ({4}_1.fq,{4}_2.fq) containing the reads that aligned in fastq format
	## --allow-contain :: keeps reads that match perfectly
	## -1 {2} :: input file 1
	## -2 {3} :: input file 2
	## --un {5} :: This will generate two files ({5}_1.fq,{5}_2.fq) containing the reads that did not align in fastq format
    os.system("bowtie -q --quiet --sam -p 1 -v 0 --ff --norc -3 4 -X 2000 -m 1 --al {4} --un {5} --sam --allow-contain {0} -1 {1} -2 {2} {3}".format(targets_dir + "targets",
              output_revcomp_R1,output_R2,output_bowtie,bowtie_results,bowtie_unmatched))
          
    # If we observe results (i.e. reads mapped)
    if os.path.isfile(os.path.splitext(bowtie_results)[0]+"_1.fq"):
        #Have to switch rev-comp of R1 back to original
        fastq_revcomp.rev_comp(bowtie_results_1,bowtie_results_1_rc)
        output_file_1 = bowtie_results_1_rc
        output_file_2 = output_dir + "bowtie_results{0}_2.fq".format(base)
    else:
        output_file_1 = None
        output_file_2 = None
    
    #Update the search list so that not trimming the same reads again that already already matched
    #revcomp Read 1 reads that were not aligned 
    if os.path.exists(bowtie_unmatched_1):
        fastq_revcomp.rev_comp(bowtie_unmatched_1,bowtie_unmatched_1_rc)
    #Delete temporary files
    shutil.rmtree(search_dir,ignore_errors=True)
    
    #Return the files containing the matched reads and updated fastq to search
    os.remove(output_R1)
    os.remove(output_R2)
    os.remove(output_revcomp_R1)
    return output_file_1, output_file_2,bowtie_unmatched_1_rc,bowtie_unmatched_2

def full_trim(trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,max_handle_len,input_R1,
                                                                input_R2,input_targets,error_rate,quality_min,quality_change,threads,mismatches):
    """execute full trimming algorithm"""
    
    
    # The algorithm employed is as follows:
    # 
    # 1.) Filter reads into two sets: Case_I = do not have any adapters at all - i.e. long insert sizes, Case_II = have some adapter.
    # 2.) For Case II, trim all possible with cutadapt (Case_II_clipped).
    # 3.) For reads that aren't trimmed from Case II (Case_II_unclipped):
    # 3.1) First filter these reads out with cutadapt (i.e. collect just Case_II_unclipped)
    # 3.2) Execute trim_search algorithm
    # 4.) Combine all reads into one file to send to spats
    
    # Keeping temporary files together, but use random suffix to prevent
    # run results from colliding
    directory_suffix = id_generator(10)
    output_dir = current_dir + "At_output_" + directory_suffix + "/"
    targets_dir = current_dir + 'targets_' + directory_suffix + "/"
    
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    #Build bowtie targets for checking rev-comp pairs
    print >> sys.stderr, "[%s] Building Bowtie targets" % (right_now())
    try:
        os.mkdir(targets_dir)
    except:
        pass
    os.system("bowtie-build -q {0} {1}targets".format(input_targets,targets_dir))
    
    #Decide what length is optimal for adapter clipping (shorter lengths save time/space)
	#Also decide the minimum length to leave at the three prime end for unique alignment to targets (if unspecified)
    print >> sys.stderr, "[%s] Determining optimal trim_match and minimum read length" % (right_now())
    trim_match,min_read_len,error_rate = trim_min_read_calculator(trim_auto,min_read_auto,input_targets,A_b_sequence,max_handle_len,trim_match,min_read_len,read_len)
	
	#trim_len = length of reads that will end up with after processing
	## Set to the read length minus the match length
    trim_len = read_len-trim_match
    
    #Define a list of files for R1 and R2 to store all of subsets of reads want to keep
    combine_files_R1 = []
    combine_files_R2 = []
    
    if quality_change != None:
        #Modify all of of the quality scores as requested by the user
        print >> sys.stderr, "[%s] Adjusting quality scores" % (right_now())
        input_R1 = quality_mod(input_R1,read_len,quality_change)
        input_R2 = quality_mod(input_R2,read_len,quality_change)
	  
    # 1.) Filter reads into two sets: Case_I = do not have any adapters at all - i.e. long insert sizes, Case_II = have some adapter.
    ## Using Bowtie
    ## Will need to use -5 4 option to get rid of 5' handle sequence on Read 1, otherwise won't match targets
    ## Will also use -3 4 option to get rid of 4nts on 3' end - this will remove dangling handle bases on Read 2 for 
    ## the special case where read does not have adapter, but starts to bleed into read 2
    print >> sys.stderr, "[%s] Filtering reads and storing Case I reads" % (right_now())
    
    bowtie_results_I = output_dir + "Case_I.fq"
    bowtie_results_I_1 = output_dir + "Case_I_1.fq"
    bowtie_results_I_2 = output_dir + "Case_I_2.fq"
    
    bowtie_results_II = output_dir + "Case_II.fq"
    bowtie_results_II_1 = output_dir + "Case_II_1.fq"
    bowtie_results_II_2 = output_dir + "Case_II_2.fq"
    
    output_bowtie = output_dir + "bowtie_output_filter_cases.sam"
    
    # Construct bowtie command
	## -q :: fastq input files
	## --quiet :: suppresses verbose output
	## -p {6} :: use {6} processors
	## -v 0 :: report hits with no mismatches - i.e. only perfect revcomps matched
	## -5 4 :: remove 4 nts from 5' end (this gets rid of handle sequence in Read 1, otherwise would not align)
	## -3 4 :: remove 4 nts from 3' end (in special case, read 2 will have a little bit of handle before the adapter - need to trim this)
	## -X 2000 :: setting huge upper limit on insert size (i.e. no limit)
	## -m 1 :: unique alignments
	## --al {4}:: This will generate two files ({4}_1.fq,{4}_2.fq) containing the reads that aligned in fastq format (Case I)
	## --un {5}:: This will generate two files ({5}_1.fq,{5}_2.fq) containing the reads that did NOT align in fastq format (Case II)
	## --allow-contain :: keeps reads that match perfectly
	## -1 {1} :: input file 1
	## -2 {2} :: input file 2
    os.system("bowtie -q --quiet --sam -p {6} -v 0 -5 4 -3 4 -X 2000 -m 1 --al {4} --un {5} --sam --allow-contain {0} -1 {1} -2 {2} {3}".format(targets_dir + "targets",
              input_R1,input_R2,output_bowtie,bowtie_results_I,bowtie_results_II,threads))
        
    # add Case I reads to the combined list
    # If we observe results Case I reads, add to combined list 
    if os.path.isfile(bowtie_results_I_1):
        combine_files_R1.append(bowtie_results_I_1)
        combine_files_R2.append(bowtie_results_I_2)
    
	# 2.) For Case II, trim all possible with cutadapt.
	    
    print >> sys.stderr, "[%s] Clipping Case II reads" % (right_now())
    
    #Clip off adapter from the raw reads
    output_II_1_clipped = os.path.splitext(bowtie_results_II_1)[0] + "_clipped.fq"
    output_II_2_clipped = os.path.splitext(bowtie_results_II_2)[0] + "_clipped.fq"
	    
    adapt_1 = A_b_sequence
    adapt_2 = A_t_sequence
    
	# Construct cutadapt command
	## Options:
	## --quiet :: Suppresses verbose output
	## --discard-untrimmed :: keep only trimmed sequences
	## -O {0} :: minimum overlap length {0} to clip sequence (need to match at least this many nt to adapter sequence)
	## -a {1} :: sequence of 3' end adapter
	## -m {2} :: minimum sequence length cutoff for sequences to keep (discard if under {2} long)
	## -o {} :: output destination  
	## -p {} :: Write reads from the paired-end input to {}
    #First pass clips one of the fastq files, and deletes any removed sequences from matching read to prevent orphan sequences
    os.system("cutadapt --quiet --discard-untrimmed -O {0} -a {1} -m {2} -o tmp{5}.1.fastq -p tmp{5}.2.fastq {3} {4}".format(trim_match,
                                                                                                                             adapt_1,
                                                                                                                             min_read_len,
                                                                                                                             bowtie_results_II_1,
                                                                                                                             bowtie_results_II_2,
                                                                                                                             directory_suffix))
    
    #Second pass clips the other fastq file, and similarly removes any deleted lines from the other paired file
    os.system("cutadapt --quiet --discard-untrimmed -O {0} -a {1} -m {2} -o {3} -p {4} tmp{5}.2.fastq tmp{5}.1.fastq".format(trim_match,
                                                                                                                             adapt_2,
                                                                                                                             min_read_len,
                                                                                                                             output_II_2_clipped,
                                                                                                                             output_II_1_clipped,
                                                                                                                             directory_suffix))
    
    #Remove the temporary files     
    os.system("rm tmp{0}.1.fastq tmp{0}.2.fastq".format(directory_suffix))
     
    # Add these reads to final set
    combine_files_R1.append(output_II_1_clipped)
    combine_files_R2.append(output_II_2_clipped)
	
    print >> sys.stderr, "[%s] Excecuting trim_search on remaining Case II reads" % (right_now())
	
    # 3.) For reads that aren't trimmed from Case II:
    
    # 3.1) First filter these reads out with cutadapt
    output_II_1_unclipped = os.path.splitext(bowtie_results_II_1)[0] + "_unclipped.fq"
    output_II_2_unclipped = os.path.splitext(bowtie_results_II_2)[0] + "_unclipped.fq"
	
	# Construct cutadapt command
	## Options:
	## --quiet :: Suppresses verbose output
	## --discard-trimmed :: keep only untrimmed sequences
	## -O {0} :: minimum overlap length {0} to clip sequence (need to match at least this many nt to adapter sequence)
	## -a {1} :: sequence of 3' end adapter
	## -m {2} :: minimum sequence length cutoff for sequences to keep (discard if under {2} long)
	## -o {} :: output destination  
	## -p {} :: Write reads from the paired-end input to {}

    #First pass clips one of the fastq files, and deletes any removed sequences from other to prevent orphan sequences
    os.system("cutadapt --quiet --discard-trimmed -O {0} -a {1} -m {2} -o tmp{5}_2.1.fastq -p tmp{5}_2.2.fastq {3} {4}".format(trim_match,
                                                                                                                               adapt_1,
                                                                                                                               min_read_len,
                                                                                                                               bowtie_results_II_1,
                                                                                                                               bowtie_results_II_2,
                                                                                                                               directory_suffix)) 
    
    #Second pass clips the other fastq file, and similarly removes any deleted lines from the other paired file
    os.system("cutadapt --quiet --discard-trimmed -O {0} -a {1} -m {2} -o {3} -p {4} tmp{5}_2.2.fastq tmp{5}_2.1.fastq".format(trim_match,
                                                                                                                               adapt_2,
                                                                                                                               min_read_len,
                                                                                                                               output_II_2_unclipped,
                                                                                                                               output_II_1_unclipped,
                                                                                                                               directory_suffix)) 
    
    #Remove the temporary files 
    os.system("rm tmp{0}_2.1.fastq tmp{0}_2.2.fastq".format(directory_suffix))    
            
    # 3.2) Execute trim_search algorithm
	
	#Now do a stepwise manual trim to search for rev comp in the window that the clipper misses
	## Start at read_len - 1 to get reads that have at least one nt of adapter
	## End at trim_len + 1 because clipper will have found things with trim_len nts of adapter in previous steps
    for base in range(read_len-1,trim_len,-1): 
        
        R1_results,R2_results,update_1,update_2 = clip_search(base,
                                                              read_len,
                                                              max_handle_len,
                                                              output_dir,
                                                              targets_dir,
                                                              output_II_1_unclipped,
                                                              output_II_2_unclipped,
                                                              threads)

        # Using unclipped reads as input for next iteration of clip_search
        output_II_1_unclipped = update_1
        output_II_2_unclipped = update_2
        
        #Add the aligned rev-comp files to the list of files to be recombined
        if R1_results and R2_results: 
            combine_files_R1.append(R1_results)
            combine_files_R2.append(R2_results)            
        
    #If allowing mismatches in Spats, include the remaining reads that were never aligned perfectly with bowtie:
    if mismatches == True:
        if os.path.exists(output_II_1_unclipped):
            combine_files_R1.append(output_II_1_unclipped)
        if os.path.exists(output_II_2_unclipped):
            combine_files_R2.append(output_II_2_unclipped)          
        
        
    # 3.) Combine all reads into one file to send to spats
    print >> sys.stderr, "[%s] Combining all processed reads" % (right_now())
    try:
        os.mkdir(final_dir)
    except:
        pass
    
    combined_R1 = final_dir + 'combined_R1.fastq'
    combined_R2 = final_dir + 'combined_R2_pretrim.fastq'
    
    if os.path.exists(combined_R1):
        os.remove(combined_R1)
    if os.path.exists(combined_R2):
        os.remove(combined_R2)
    
    # JBL - define combine_files function
    outfile = file(combined_R1,'w')
    for fname in combine_files_R1:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
        infile.close()
    outfile.close()
    
    outfile = file(combined_R2,'w')
    for fname in combine_files_R2:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
            infile.close()
    outfile.close()
    
    # Final set of reads trimmed down by 4nts on 3' end of R2 to remove handle (YYYR,RRRY) that is in R2 for short inserts
	## Construct cutadapt command
	## -t N = trim N nts from the end of the read
    final_R2 = final_dir + 'combined_R2.fastq'
    os.system("cutadapt --quiet -u -4 {0} > {1}".format(combined_R2,final_R2))
    
    
    #Delete the rest of the temporary files
    shutil.rmtree(output_dir,ignore_errors=True)
    shutil.rmtree(targets_dir,ignore_errors=True)
    os.remove(combined_R2)


def main(argv=None,):
    
    # Gather params and then call full_trim() - separated out for easy functional testing
    params = Params()
        
    try:
        if argv is None:
            argv = sys.argv
            [args,trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,
                    min_read_auto,max_handle_len,error_rate,quality_min,quality_change,threads,mismatches] = params.parse_options(argv)
            params.check()
        
        input_R1 = args[0]
        input_R2 = args[1]
        input_targets = args[2]
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning Adapter_trimmer reads processing (v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------------"
        start_time = datetime.now()
        
        full_trim(trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,
                    max_handle_len,input_R1,input_R2,input_targets,error_rate,quality_min,quality_change,threads,mismatches)
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Run complete [%s elapsed]" % formatTD(duration)
        print >> sys.stderr
    
        print final_dir
    
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    # JBL - need more elegant exit than this
    # print result returned by main() to be able to capture in bash
    main()
    sys.exit(0)

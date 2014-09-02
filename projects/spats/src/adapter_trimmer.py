#!/usr/bin/env python

import targets_analyzer
import sys
import os
import getopt
import shutil
from datetime import datetime, date, time
import string
import random
import re

current_dir = os.getcwd()+'/'

help_message = '''
Adapter_trim trims adapter sequences from fastq files in preparatino for processing by spats.

takes fastq files, searches for short sequences that are reverse complimentary,
trims off adapter sequences, and prepares reads for processing directly with spats

Usage:
   Adapter_trim [options] <R1_seq.fastq> <R2_seq.fastq> <targets.fa>

Options:
-h, --help                  opens help message
-v, --version               displays version number
--trim-match <N>            number of nucleotides of adapter to search for with clipper (default = 6) *Adjust trim length this way.
--read-len <N>              Number of length of reads (ex. 35 for 2x35 bp paired end reads) (default = 35)
--min-read-len <n>          Number of nt on the 3' end to leave after (w/ adapter) (default = 6)
--A-b-sequence <string>     A_adapter_b sequence (default: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG)
--A-t-sequence <string>     A_adapter_t sequence (default: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)
--max-handle-len <N>        Number of nucleotides in the longest (+/-) handle (default = 4 for RRRY/YYYR)
-o, --output <string>       If running many simultaneously runs in the same folder, use this option to name the treated reads
                            folder to avoid overwriting during processing

JBL - if trim-match or min-read-len are not set, they will automatically be determined

The algorithm employed is as follows:
 
1.) Filter reads into two sets: Case_I = do not have any adapters at all - i.e. long insert sizes, Case_II = have some adapter.
2.) For Case II, trim all possible with fastx_spats_clipper (Case_II_clipped).
3.) For reads that aren't trimmed from Case II (Case_II_unclipped):
3.1) First filter these reads out with fastx_spats_clipper (i.e. collect just Case_II_unclipped)
3.2) Execute trim_search algorithm
4.) Combine all reads into one file to send to spats

'''

def get_version():
    return "0.0.3"

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
                                        "trim-match=",
                                        "read-len=",
                                        "A-b-sequence=",
                                        "A-t-sequence=",
                                        "min-read-len=",
                                        "output=",
                                        "max-handle-len="])
        
        except getopt.error, msg:
            raise Usage(msg)
        
        trim_match = 6
        trim_auto = True
        read_len = 35
        A_b_sequence = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
        A_t_sequence = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        final_dir = current_dir + "combined_reads_" + id_generator(10) + "/"
        max_handle_len = 4
        min_read_auto = True
        min_read_len = 6
        for option, value in opts:
            if option in ("-v", "--version"):
                print "Adapter_trimmer v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("--trim-match"):
                trim_auto = False
                trim_match = int(value)
            if option in ("--read-len"):
                read_len = int(value)
            if option in ("--A-b-sequence"):
                A_b_sequence = value
            if option in ("--A-t-sequence"):
                A_t_sequence = value
            if option in ("--min-read-len"):
                min_read_len = int(value)
                min_read_auto = False
            if option in ("-o","--output"):
                final_dir = value
            if option in ("--max-handle-len"):
                max_handle_len = int(value)
        # JBL - A_t shows up as a revcomp in read2 because of library format
        A_t_sequence = reverse_complement(A_t_sequence)
        
        if len(args) != 3:
            raise Usage(help_message)
        return args,trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,max_handle_len
    
    def check(self):
        pass

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    '''From http://stackoverflow.com/questions/2257441/python-random-string-generation-with-upper-case-letters-and-digits'''
    return ''.join(random.choice(chars) for x in range(size))

def trim_min_read_calculator(trim_auto,min_len_auto,input_targets,A_b_sequence,max_handle_len,trim_match,min_read_len,read_len):
    
    #Call Ray's script to find the minimal lengths of A_b to search and number of nucleotides to trim from 3' end
    uniqueness,A_b_search_len = targets_analyzer.readFASTA(input_targets,A_b_sequence)
    
    if min_len_auto == True:
        # Add handle length to uniqueness length to find min_read_length considering
        min_read_len = max_handle_len + uniqueness 
        print >> sys.stderr, "[%s] Automatically dropping %s nt from end of RNAs (from 3' end analysis)" % (right_now(),uniqueness)
    else:
        print >> sys.stderr, "[%s] Manually dropping %s nt from end of reads" % (right_now(),min_read_len)
    
    if trim_auto == True:
        if A_b_search_len >= read_len:
            trim_match = read_len - min_read_len  #
        else:
            trim_match = A_b_search_len + 1  #Adds two greater than the longest A_b sequence found
        print >> sys.stderr, "[%s] Automatically set trim-match to %s based on targets file" % (right_now(),trim_match)
    else:
         print >> sys.stderr, "[%s] Trim-match manually set to %s" % (right_now(),trim_match)
    
    return trim_match,min_read_len

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def quality_mod(inputfile,read_len):
    #Modifies all of the quality reads to 'I' so that fastx ignores qualities
    
    #qualmodfile = inputfile[:-3] + '_Qs.fq'
    qualmodfile = os.path.splitext(inputfile)[0] + '_Qs.fq'
    qualstring = "I"*read_len
    qualcommand = "awk \'BEGIN{s=0} {s += 1; if (s % 4 == 0) print \"" + qualstring + "\"; else print }\'"
    #this is needed to except the curly brackets in the below command
    os.system(qualcommand +" {0} > {1}".format(inputfile, qualmodfile));
    #print("done")
    return qualmodfile

def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)

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

def relabel_reads(input_R1, input_R2, treated_handle, untreated_handle, current_dir,relabel=True):
    #filter_cmd = ["prep_reads"]
	# relabels reads with names of increasing integers
    print >> sys.stderr, "[%s] Relabeling reads in %s and %s" % (right_now(), input_R1, input_R2)
    
    #filter_log = open(logging_dir + "relabel_reads.log", "w")
    
    filter_cmd = "relabel_reads --output-dir {0} --phred33-quals --fastq {1} {2}".format(current_dir,input_R1,input_R2)
    os.system(filter_cmd)
    
    #print "\t executing: `%s'" % " ".join(filter_cmd)

def match_read_pairs(left_in_reads, right_in_reads, left_out_reads, right_out_reads):
    
    #filter_cmd = ["prep_reads"]
    print >> sys.stderr, "[%s] Rematching read pairs" % (right_now())
    
    #filter_log = open(logging_dir + "relabel_reads.log", "w")
    
    cmd = ["match_read_pairs"]
    cmd.extend([left_in_reads,
                right_in_reads,
                left_out_reads,
                right_out_reads])
    
    #print "\t executing: `%s'" % " ".join(cmd)
    ret = os.system(" ".join(cmd))
    
    return [left_out_reads, right_out_reads]


def find_ids(filelist):
    """Return all read id's present in filelist"""
    id_list = []
    for fname in filelist:
        with open(fname) as infile:
            for line in infile:
                if re.match('^@',line):
                    read_id = str(line[1:-1])
                    id_list.append(read_id)
        infile.close()
    return id_list

def remove_reads(read_file,id_list):
    """remove all reads in id_list from read_file
    
    returns name of new file
    """
    new_read_file = read_file[:-3] + "_ids_removed.fastq"

    # Make a full copy of the id_list 
    # (We'll be removing ids from this list as we go, but we don't want to modify the master list)
    ids_to_remove = id_list[:]
    
    read_chunk = [] #Each read consists of 4 lines
    outfile = file(new_read_file,'w')
    with open(read_file) as infile:
        for line in infile:
            
            # Grab read chunks
            read_chunk.append(line)
            # A read_chunk is 4 lines
            if len(read_chunk) == 4:
                
                #Grab id for this read (on first line of chunk)
                read_id = str(read_chunk[0][1:-1])
                
                #Test to see if id is in id_list
                remove = False
                if read_id in ids_to_remove:
                    # If found, set flag so don't print out and remove id from the list
                    remove = True
                    ids_to_remove.remove(read_id) # once removed don't need to search again
                
                if not remove:
                    #If don't remove, then write this read_chunk
                    for el in read_chunk:
                        outfile.write(el)
                    
                #reset read_chunk to go to next chunk
                read_chunk = None # clear memory
                read_chunk = []
    outfile.close()
    return new_read_file


def trim_search(base,max_handle_len,output_dir,targets_dir,input_R1,input_R2):
    """Search for revcomp reads using bowtie after trimming to length base
    
    Assuming input_R1 and input_R2 are filenames of files that have been quality_mod()'d
    
    """
    print >> sys.stderr, "[%s] Finding reads that are %s nt long" % (right_now(), (base-max_handle_len)) #JBLQ - what does this mean?
    
    #Make a temporary folder for intermediate files to delete after each trim
    search_dir = current_dir + "search_" + id_generator(10) + "/"
    try:
        os.mkdir(search_dir)
    except:
        pass # should terminate here since need this search_dir to proceed
    
    #Implement Adapter Clipping Step 1: Find reads that the clipper will miss
    ## Trim base-by-base and look for revcomp sequences. 3 cases:
    ## 1) Reads contain no adapter, thus will not revcomp and will not keep these mate pairs.
    ## 2) Reads contain some adapter, will not revcomp and will not keep these mate pairs.
    ## 3) Adapter has been manually trimmed out, reads will revcomp, keep this pair.
    ## Accumulate kept reads in files - will combine all reads at the end of analysis.
    
    #Block of filename constructions for outputfiles
    output_R1 = os.path.splitext(input_R1)[0] + "_" + str(base) + ".fq"
    output_R2 = os.path.splitext(input_R2)[0] + "_" + str(base) + ".fq"
    output_revcomp_R1 = os.path.splitext(input_R1)[0] + "_" + str(base) + "_revcomp" + ".fq"
    output_bowtie = output_dir + "bowtie_output_" + str(base) + ".sam"
    bowtie_results = output_dir + "bowtie_results{0}.fq".format(base)
    bowtie_results_rc = output_dir + "bowtie_results{0}_1_rc.fq".format(base)
    
	#Trim down each file to base
	## Construct fastx_trimmer command
	## -l N = Last base to keep. Default is entire read.
    os.system("fastx_trimmer -l {0} -i {1} -o {2}".format(base,input_R1,output_R1))
    os.system("fastx_trimmer -l {0} -i {1} -o {2}".format(base,input_R2,output_R2))
	
	# Reverse complement read 1
    os.system("fastx_reverse_complement -i {0} -o {1}".format(output_R1,output_revcomp_R1))
	
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
    os.system("bowtie -q --quiet --sam -p 1 -v 0 --ff --norc -3 4 -X 2000 -m 1 --al {4} --sam --allow-contain {0} -1 {1} -2 {2} {3}".format(targets_dir + "targets",
              output_revcomp_R1,output_R2,output_bowtie,bowtie_results))
    
    # If we observe results (i.e. reads mapped)
    if os.path.isfile(os.path.splitext(bowtie_results)[0]+"_1.fq"):
        #Have to switch rev-comp of R1 back to original
        os.system("fastx_reverse_complement -i {0} -o {1}".format((os.path.splitext(bowtie_results)[0]+"_1.fq"),bowtie_results_rc))    
        output_file_1 = bowtie_results_rc
        output_file_2 = output_dir + "bowtie_results{0}_2.fq".format(base)
    else:
        output_file_1 = None
        output_file_2 = None
    
    #Delete temporary files
    shutil.rmtree(search_dir,ignore_errors=True)
    
    #Return the files containing the matched reads
    os.remove(output_R1)
    os.remove(output_revcomp_R1)
    os.remove(output_R2)
    return output_file_1, output_file_2


def full_trim(trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,max_handle_len,input_R1,input_R2,input_targets):
    """execute full trimming algorithm"""
    
    
    # The algorithm employed is as follows:
    # 
    # 1.) Filter reads into two sets: Case_I = do not have any adapters at all - i.e. long insert sizes, Case_II = have some adapter.
    # 2.) For Case II, trim all possible with fastx_spats_clipper (Case_II_clipped).
    # 3.) For reads that aren't trimmed from Case II (Case_II_unclipped):
    # 3.1) First filter these reads out with fastx_spats_clipper (i.e. collect just Case_II_unclipped)
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
    trim_match,min_read_len = trim_min_read_calculator(trim_auto,min_read_auto,input_targets,A_b_sequence,max_handle_len,trim_match,min_read_len,read_len)
	
	#trim_len = length of reads that will end up with after processing
	## Set to the read length minus the match length
    trim_len = read_len-trim_match
    
    #Define a list of files for R1 and R2 to store all of subsets of reads want to keep
    combine_files_R1 = []
    combine_files_R2 = []
    
    #Modify all of of the quality scores to 'I' so that the fastx tools don't get hung up
    print >> sys.stderr, "[%s] Adjusting quality scores" % (right_now())
    input_R1 = quality_mod(input_R1,read_len)
    input_R2 = quality_mod(input_R2,read_len)
	
	#Relabel reads before clipping
	## This will allow us to detect a mate pair where only one of the pair is
	## clipped or dropped (i.e. other side has a mismatch so is not detected) by fastx_spats_clipper below.
	## match_read_pairs will be used later to throw out these anomolies and pair up reads with same indexes
	
	#Generates files NOMASK_1.fq and NOMASK_2.fq
    relabel_reads(input_R1,input_R2,None,None,output_dir)
    input_R1 = output_dir+"NOMASK_1.fq"
    input_R2 = output_dir+"NOMASK_2.fq"
    
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
    os.system("bowtie -q --quiet --sam -p 1 -v 0 -5 4 -3 4 -X 2000 -m 1 --al {4} --un {5} --sam --allow-contain {0} -1 {1} -2 {2} {3}".format(targets_dir + "targets",
              input_R1,input_R2,output_bowtie,bowtie_results_I,bowtie_results_II))
        
    # add Case I reads to the combined list
    # If we observe results Case I reads, add to combined list 
    if os.path.isfile(bowtie_results_I_1):
        combine_files_R1.append(bowtie_results_I_1)
        combine_files_R2.append(bowtie_results_I_2)
    
	# 2.) For Case II, trim all possible with fastx_spats_clipper.
	    
    print >> sys.stderr, "[%s] Clipping Case II reads" % (right_now())
    
    #Clip off adapter from the raw reads
    output_II_1_clipped = os.path.splitext(bowtie_results_II_1)[0] + "_clipped.fq"
    output_II_2_clipped = os.path.splitext(bowtie_results_II_2)[0] + "_clipped.fq"
	
	# Construct fastx_spats_clipper command
	## Options:
	## -c :: Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter).
	## -M trim_match :: require minimum adapter length of trim_match (if < trim_match align with adapter, don't clip it)
	## -a :: adapter string
	## -l {4} :: discard sequences shorter than {4} nts
	## -n :: keep sequences with N's
    os.system("fastx_spats_clipper -c -M {3} -a {0} -l {4} -n -i {1} -o {2}".format(A_b_sequence,bowtie_results_II_1,output_II_1_clipped,trim_match,min_read_len))
    
    # Actually using revcomp(A_t_sequence) since read 2 reads the revcomp of A_t
    os.system("fastx_spats_clipper -c -M {3} -a {0} -l {4} -n -i {1} -o {2}".format(A_t_sequence,bowtie_results_II_2,output_II_2_clipped,trim_match,min_read_len))
    
    #rematch read pairs immediately
	## match_read_pairs throws out any unique id reads - i.e. throws out an orphaned mate read
	## orphan reads generated if sequencing error in adapter sequence in one of the mates prevents fastx_spats_clipper from matching and clipping
    matched_II_1_clipped = os.path.splitext(bowtie_results_II_1)[0] + "_clipped_kept.fq"
    matched_II_2_clipped = os.path.splitext(bowtie_results_II_2)[0] + "_clipped_kept.fq"
    matched_reads = match_read_pairs(output_II_1_clipped, output_II_2_clipped, matched_II_1_clipped, matched_II_2_clipped)
    
    # Add these reads to final set
    combine_files_R1.append(matched_II_1_clipped)
    combine_files_R2.append(matched_II_2_clipped)
	
    print >> sys.stderr, "[%s] Excecuting trim_search on remaining Case II reads" % (right_now())
	
    # 3.) For reads that aren't trimmed from Case II:
    
    # 3.1) First filter these reads out with fastx_spats_clipper
    output_II_1_unclipped = os.path.splitext(bowtie_results_II_1)[0] + "_unclipped.fq"
    output_II_2_unclipped = os.path.splitext(bowtie_results_II_2)[0] + "_unclipped.fq"
	
	# Construct fastx_spats_clipper command
	## Options:
	## -C :: Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).
	## -M trim_match :: require minimum adapter length of trim_match (if < trim_match align with adapter, don't clip it)
	## -a :: adapter string
	## -l {4} :: discard sequences shorter than {4} nts
	## -n :: keep sequences with N's
    os.system("fastx_spats_clipper -C -M {3} -a {0} -l {4} -n -i {1} -o {2}".format(A_b_sequence,bowtie_results_II_1,output_II_1_unclipped,trim_match,min_read_len))
    
    # Actually using revcomp(A_t_sequence) since read 2 reads the revcomp of A_t
    os.system("fastx_spats_clipper -C -M {3} -a {0} -l {4} -n -i {1} -o {2}".format(A_t_sequence,bowtie_results_II_2,output_II_2_unclipped,trim_match,min_read_len))
    
    #rematch read pairs immediately
	## match_read_pairs throws out any unique id reads - i.e. throws out an orphaned mate read
	## orphan reads generated if sequencing error in adapter sequence in one of the mates prevents fastx_spats_clipper from matching and clipping
    matched_II_1_unclipped = os.path.splitext(bowtie_results_II_1)[0] + "_unclipped_kept.fq"
    matched_II_2_unclipped = os.path.splitext(bowtie_results_II_2)[0] + "_unclipped_kept.fq"
    matched_reads = match_read_pairs(output_II_1_unclipped, output_II_2_unclipped, matched_II_1_unclipped, matched_II_2_unclipped)    
    
    # 3.2) Execute trim_search algorithm
	
	#Now do a stepwise manual trim to search for rev comp in the window that the clipper misses
	## Start at read_len - 1 to get reads that have at least one nt of adapter
	## End at trim_len + 1 because clipper will find things with trim_len nts of adapter in later steps
    for base in range(read_len-1,trim_len,-1): 
        
        R1_results, R2_results = trim_search(base,max_handle_len,output_dir,targets_dir,matched_II_1_unclipped,matched_II_2_unclipped)
        
        #Add the aligned rev-comp files to the list of files to be recombined
        if R1_results and R2_results: 
            combine_files_R1.append(R1_results)
            combine_files_R2.append(R2_results)            
        
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
	## Construct fastx_trimmer command
	## -t N = trim N nts from the end of the read
    final_R1 = combined_R1
    final_R2 = final_dir + 'combined_R2.fastq'
    os.system("fastx_trimmer -t 4 -i {0} -o {1}".format(combined_R2,final_R2))
	
    
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
            args,trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,max_handle_len = params.parse_options(argv)
            params.check()
        
        input_R1 = args[0]
        input_R2 = args[1]
        input_targets = args[2]
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning Adapter_trimmer reads processing (v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------------"
        start_time = datetime.now()
        
        full_trim(trim_match,read_len,A_b_sequence,A_t_sequence,min_read_len,final_dir,trim_auto,min_read_auto,max_handle_len,input_R1,input_R2,input_targets)
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Run complete [%s elapsed]" % formatTD(duration)
        print >> sys.stderr
        
        return final_dir
    
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":
    # JBL - need more elegant exit than this
    # print result returned by main() to be able to capture in bash
    print main()
    sys.exit(0)

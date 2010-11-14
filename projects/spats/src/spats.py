#! /usr/bin/env python
# encoding: utf-8
"""
spats.py

Created by Cole Trapnell on 2010-4-11.
Copyright (c) 2010 Cole Trapnell. All rights reserved.
"""

import sys
try:
    import psyco
    psyco.full()
except ImportError:
    pass

import sys
import getopt
import subprocess
import errno
import os
import tempfile
import warnings
import shutil
import copy
from datetime import datetime, date, time

use_message = '''
 Spats builds reactivity profiles from SHAPE-Seq experiments.
 
 Usage:
     spats [options] <rna.fasta> <treated_handle> <untreated_handle> <reads1[,reads2,...,readsN]> <reads1[,reads2,...,readsN]>
     
 Options:
     -o/--output-dir                <string>    [ default: ./spats_out ] 
     --left-adapter                 <string>    [ default: None ]
     --right-adapter                <string>    [ default: None ]
     
SAM Header Options (for embedding sequencing run metadata in output):
    --rg-id                        <string>    (read group ID)
    --rg-sample                    <string>    (sample ID)
    --rg-library                   <string>    (library ID)
    --rg-description               <string>    (descriptive string, no tabs allowed)
    --rg-platform-unit             <string>    (e.g Illumina lane ID)
    --rg-center                    <string>    (sequencing center name)
    --rg-date                      <string>    (ISO 8601 date of the sequencing run)
    --rg-platform                  <string>    (Sequencing platform descriptor) 
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

output_dir = "./spats_out/"
logging_dir = output_dir + "logs/"
run_log = None
run_cmd = None

tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"

#ok_str = "\t\t\t\t[OK]\n"
fail_str = "\t[FAILED]\n"


class SpatsParams:
        
    class SystemParams:
        def __init__(self,
                     bowtie_threads,
                     keep_tmp):
            self.bowtie_threads = bowtie_threads
            self.keep_tmp = keep_tmp
            
        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-p", "--num-threads"):
                    self.bowtie_threads = int(value)
                if option in ("--keep-tmp"):
                    self.keep_tmp = True
        
        def check(self):
            pass
        
    
    class ReadParams:
        def __init__(self,
                     solexa_quals,
                     phred64_quals,
                     seed_length,
                     left_adapter,
                     right_adapter,
                     reads_format,
                     read_group_id,
                     sample_id,
                     library_id,
                     description,
                     seq_platform_unit,
                     seq_center,
                     seq_run_date,
                     seq_platform):
            self.solexa_quals = solexa_quals
            self.phred64_quals = phred64_quals
            self.seed_length = seed_length
            self.left_adapter = left_adapter
            self.right_adapter = right_adapter
            self.reads_format = reads_format
            self.read_group_id = read_group_id 
            self.sample_id = sample_id
            self.library_id = library_id
            self.description = description
            self.seq_platform_unit = seq_platform_unit
            self.seq_center = seq_center
            self.seq_run_date = seq_run_date
            self.seq_platform = seq_platform
            
        def parse_options(self, opts):
            for option, value in opts:
                if option == "--solexa-quals":
                    self.solexa_quals = True
                if option in ("--solexa1.3-quals", "--phred64-quals"):
                    self.phred64_quals = True    
                if option in ("-s", "--seed-length"):
                    self.seed_length = int(value)
                if option in ("--left-adapter"):
                    self.left_adapter = value
                if option in ("--right-adapter"):
                    self.right_adapter = value
                if option == "--rg-id":
                    self.read_group_id = value
                if option == "--rg-sample":
                    self.sample_id = value
                if option == "--rg-library":
                    self.library_id = value
                if option == "--rg-description":
                    self.description = value
                if option == "--rg-platform-unit":
                    self.seq_platform_unit = value
                if option == "--rg-center":
                    self.seq_center = value
                if option == "--rg-date":
                    self.seq_run_date = value    
                if option == "--rg-platform":
                    self.seq_platform = value            

        def check(self):
            if self.seed_length != None and self.seed_length < 20:
                print >> sys.stderr, "Error: arg to --seed-length must be at least 20"
                exit(1)
            if (not self.read_group_id and self.sample_id) or (self.read_group_id and not self.sample_id):
                print >> sys.stderr, "Error: --rg-id and --rg-sample must be specified or omitted together"
                exit(1)
            
            #TODO validate adapters for nucleotide chars
                
                    
    def __init__(self):        
        
        self.read_params = self.ReadParams(False,               # solexa_scale
                                           False,
                                           None,                # seed_length
                                           None,
                                           None,
                                           "fastq",             # quality_format
                                           None,                # read group id
                                           None,                # sample id
                                           None,                # library id
                                           None,                # description
                                           None,                # platform unit (i.e. lane)
                                           None,                # sequencing center
                                           None,                # run date
                                           None)                # sequencing platform
        
        self.system_params = self.SystemParams(1,               # bowtie_threads
                                               False)           # keep_tmp   
        

        self.skip_check_reads = False
        self.segment_length = 25
        self.segment_mismatches = 3
        
    def check(self):
        self.read_params.check()
        self.system_params.check()
       
        if self.segment_length <= 4:
            print >> sys.stderr, "Error: arg to --segment-length must at least 4"
            exit(1)
        if self.segment_mismatches < 0 or self.segment_mismatches > 3:
            print >> sys.stderr, "Error: arg to --segment-mismatches must in [0, 3]"
            exit(1)
        
    def cmd(self):
        cmd = ["--output-dir", output_dir,
               "--segment-length", str(self.segment_length),
               "--segment-mismatches", str(self.segment_mismatches)]
        
        if self.read_params.solexa_quals == True:
            cmd.append("--solexa-quals")
        if self.read_params.phred64_quals == True:
            cmd.append("--phred64-quals")
        return cmd
        
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvp:m:o:", 
                                        ["version",
                                         "help",  
                                         "output-dir=",
                                         "solexa-quals",
                                         "solexa1.3-quals",
                                         "phred64-quals",
                                         "num-threads=",
                                         "splice-mismatches=",
                                         "max-multihits=",
                                         "skip-check-reads",
                                         "segment-length=",
                                         "segment-mismatches=",
                                         "keep-tmp",
                                         "left-adapter=",
                                         "right-adapter=",
                                         "rg-id=",
                                         "rg-sample=",
                                         "rg-library=",
                                         "rg-description=",
                                         "rg-platform-unit=",
                                         "rg-center=",
                                         "rg-date=",
                                         "rg-platform="])
        except getopt.error, msg:
            raise Usage(msg)
            
        self.system_params.parse_options(opts)
        self.read_params.parse_options(opts)
        
        # option processing
        for option, value in opts:
            if option in ("-v", "--version"):
                print "Spats v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option == "--skip-check-reads":
                self.skip_check_reads = True
            if option == "--segment-length":
                self.segment_length = int(value)
            if option == "--segment-mismatches":
                self.segment_mismatches = int(value)
            if option in ("-o", "--output-dir"):
                global output_dir
                global logging_dir
                global tmp_dir
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
                tmp_dir = output_dir + "tmp/"
            
        if len(args) < 5:
            raise Usage(use_message)
        return args

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def prepare_output_dir():
    
    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)
        
    if os.path.exists(tmp_dir):
        pass
    else:        
        os.mkdir(tmp_dir)

def check_bowtie_index(idx_prefix):
    print >> sys.stderr, "[%s] Checking for Bowtie index files" % right_now()
    
    idx_fwd_1 = idx_prefix + ".1.ebwt"
    idx_fwd_2 = idx_prefix + ".2.ebwt"
    idx_rev_1 = idx_prefix + ".rev.1.ebwt"
    idx_rev_2 = idx_prefix + ".rev.2.ebwt"
    
    if os.path.exists(idx_fwd_1) and \
       os.path.exists(idx_fwd_2) and \
       os.path.exists(idx_rev_1) and \
       os.path.exists(idx_rev_2):
        return 
    else:
        bowtie_idx_env_var = os.environ.get("BOWTIE_INDEXES")
        if bowtie_idx_env_var == None:
            print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
            exit(1)
        idx_prefix = bowtie_idx_env_var + idx_prefix 
        idx_fwd_1 = idx_prefix + ".1.ebwt"
        idx_fwd_2 = idx_prefix + ".2.ebwt"
        idx_rev_1 = idx_prefix + ".rev.1.ebwt"
        idx_rev_2 = idx_prefix + ".rev.2.ebwt"
        
        if os.path.exists(idx_fwd_1) and \
           os.path.exists(idx_fwd_2) and \
           os.path.exists(idx_rev_1) and \
           os.path.exists(idx_rev_2):
            return 
        else:
            print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
            exit(1)

def bowtie_idx_to_fa(idx_prefix):
    idx_name = idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Reconstituting reference FASTA file from Bowtie index" % (right_now())
    
    try:    
        tmp_fasta_file_name = output_dir + idx_name + ".fa"
        tmp_fasta_file = open(tmp_fasta_file_name, "w")

        inspect_log = open(logging_dir + "bowtie_inspect_recons.log", "w")

        inspect_cmd = ["bowtie-inspect",
                       idx_prefix]
        #print >> sys.stderr, "Executing: " + " ".join(inspect_cmd) + " > " + tmp_fasta_file_name   
        ret = subprocess.call(inspect_cmd, 
                              stdout=tmp_fasta_file,
                              stderr=inspect_log)

        # Bowtie reported an error
        if ret != 0:
           print >> sys.stderr, fail_str, "Error: bowtie-inspect returned an error"
           exit(1)
           
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-inspect not found on this system.  Did you forget to include it in your PATH?"
  
    return tmp_fasta_file_name

def check_fasta(idx_prefix):
    print >> sys.stderr, "[%s] Checking for reference FASTA file" % right_now()
    idx_fasta = idx_prefix + ".fa"
    if os.path.exists(idx_fasta):
        return idx_fasta
    else:
        idx_name = idx_prefix.split('/')[-1]
        bowtie_idx_env_var = os.environ.get("BOWTIE_INDEXES")
        if bowtie_idx_env_var != None:
            idx_fasta = bowtie_idx_env_var + idx_prefix + ".fa" 
            if os.path.exists(idx_fasta):
                return idx_fasta
        
        print >> sys.stderr, "\tWarning: Could not find FASTA file " + idx_fasta
        idx_fa = bowtie_idx_to_fa(idx_prefix)
        return idx_fa
        #print >> sys.stderr, "Error: Could not find Maq binary fasta file " + idx_bfa
        #exit(1)
    
def check_index(idx_prefix):
    check_bowtie_index(idx_prefix)
    ref_fasta_file = check_fasta(idx_prefix)
    
    return (ref_fasta_file, None)

def get_bowtie_version():
    try:
        # Launch Bowtie to capture its version info
        proc = subprocess.Popen(['bowtie', '--version'],stdout=subprocess.PIPE)
        stdout_value = proc.communicate()[0]
        bowtie_version = None
        bowtie_out = repr(stdout_value)

        # Find the version identifier
        version_str = "bowtie version "
        ver_str_idx = bowtie_out.find(version_str)
        if ver_str_idx != -1:
            nl = bowtie_out.find("\\n", ver_str_idx)
            version_val = bowtie_out[ver_str_idx + len(version_str):nl]
            dash = version_val.find("-")
            if dash != -1:
                version_val = version_val[:dash]
            bowtie_version = [int(x) for x in version_val.split('.')]
        if len(bowtie_version) == 3:
            bowtie_version.append(0)
        
        return bowtie_version
    except OSError, o:
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
       exit(1)

def check_bowtie():
    print >> sys.stderr, "[%s] Checking for Bowtie" % right_now()
    bowtie_version = get_bowtie_version()
    #print >> sys.stderr, "found ", bowtie_version
    if bowtie_version == None:
        print >> sys.stderr, "Error: Bowtie not found on this system"
        exit(1)
    elif bowtie_version[0] == 0 and bowtie_version[1] < 10:
        print >> sys.stderr, "Error: Spats requires Bowtie 0.10.0 or later"
        exit(1)
    print >> sys.stderr, "\tBowtie version:\t\t %s" % ".".join([str(x) for x in bowtie_version])
        


def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds) 

def relabel_reads(params, handle_reads_list, nonhandle_reads_list, treated_handle, untreated_handle):    
    #filter_cmd = ["prep_reads"]
    print >> sys.stderr, "[%s] Relabeling reads in %s and %s" % (right_now(), handle_reads_list, nonhandle_reads_list)
    
    #filter_log = open(logging_dir + "relabel_reads.log", "w")
    
    filter_cmd = [bin_dir + "relabel_reads"]
    filter_cmd.extend(params.cmd())
    if params.read_params.reads_format == "fastq":
        filter_cmd += ["--fastq"]
    elif params.read_params.reads_format == "fasta":
        filter_cmd += ["--fasta"]
    filter_cmd.append(handle_reads_list)
    filter_cmd.append(nonhandle_reads_list) 
    filter_cmd.append(treated_handle)
    filter_cmd.append(untreated_handle) 
    
    #print "\t executing: `%s'" % " ".join(filter_cmd)    
    # files = reads_list.split(',')
    # for reads_file in files:
    try:       
        print >> run_log, " ".join(filter_cmd)
        ret = subprocess.call(filter_cmd)
                              # Bowtie reported an error
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute relabel_reads"
            exit(1)
    # prep_reads not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: relabel_reads not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
    
def match_read_pairs(params, left_in_reads, right_in_reads, left_out_reads, right_out_reads):
    #filter_cmd = ["prep_reads"]
    print >> sys.stderr, "[%s] Rematching read pairs" % (right_now())
    reads_suffix = ".fq"
    left_out_reads_filename = tmp_dir + left_out_reads + reads_suffix
    right_out_reads_filename = tmp_dir + right_out_reads + reads_suffix

    
    #filter_log = open(logging_dir + "relabel_reads.log", "w")
    
    cmd = [bin_dir + "match_read_pairs"]
    cmd.extend([left_in_reads, 
                right_in_reads, 
                left_out_reads_filename, 
                right_out_reads_filename])
       
    # print "\t executing: `%s'" % " ".join(cmd)    
    # files = reads_list.split(',')
    # for reads_file in files:
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
                              # Bowtie reported an error
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute match_read_pairs"
            exit(1)
    # prep_reads not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: match_read_pairs not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
        
    return [left_out_reads_filename, right_out_reads_filename]
    
def trim_read_adapters(params, 
                       five_prime_adapter, 
                       three_prime_adapter, 
                       reads_file, 
                       output_name):    
    #filter_cmd = ["prep_reads"]
      
    print >> sys.stderr, "[%s] Trimming adapters from reads in %s" % (right_now(), reads_file)
    reads_suffix = ".fq"
    trimmed_reads_filename = tmp_dir + output_name + reads_suffix
    
    if os.path.exists(trimmed_reads_filename):
        os.remove(trimmed_reads_filename)
    trimmed_reads = open(trimmed_reads_filename, "a")
    
    trim_log = open(logging_dir + "trim_reads.log", "a")
    try:  
        if five_prime_adapter != None and three_prime_adapter == None:
            cmd = ["fastx_clipper"]
            cmd.extend(["-a", five_prime_adapter])
            cmd.extend(["-i", reads_file])
            cmd.extend(["-o",trimmed_reads_filename])
            print >> run_log, " ".join(cmd)
            ret = subprocess.call(cmd, 
                                  stdout=trimmed_reads,
                                  stderr=trim_log)
        elif five_prime_adapter == None and three_prime_adapter != None:
            cmd = ["fastx_clipper"]
            cmd.extend(["-a", three_prime_adapter])
            cmd.extend(["-i", reads_file])
            cmd.extend(["-o", trimmed_reads_filename])
            print >> run_log, " ".join(cmd)
            ret = subprocess.call(cmd, 
                                  stdout=trimmed_reads,
                                  stderr=trim_log)
        elif five_prime_adapter != None and three_prime_adapter != None:
            five_prime_cmd = ["fastx_clipper"]
            five_prime_cmd.extend(["-a", five_prime_adapter])
            five_prime_cmd.extend(["-i", reads_file])
            
            three_prime_cmd = ["fastx_clipper"]
            three_prime_cmd.extend(["-a", three_prime_adapter])
            three_prime_cmd.extend(["-o", trimmed_reads_filename])
            print >> run_log, " ".join(five_prime_cmd) + " | " + " ".join(three_prime_cmd)
            five_prime_cmd_proc = subprocess.Popen(five_prime_cmd, 
                                                   stdout=subprocess.PIPE, 
                                                   stderr=trim_log)
            three_prime_cmd_proc = subprocess.Popen(three_prime_cmd, 
                                                    stdin=five_prime_cmd_proc.stdout, 
                                                    stderr=trim_log)
            # wait for the whole pipe to finish
            three_prime_cmd_proc.communicate()

        else:
            print >> sys.stderr, fail_str, "Error: trim_read_adapters() called with empty adapter strings"
            exit(1)
            
    # fastx_clipper not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fastx_clipper not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
        
    return trimmed_reads_filename

def bowtie(params,
           bwt_idx_prefix,
           left_reads,
           right_reads,
           reads_format,
           mapped_reads,
           unmapped_reads,
           phred_thresh=70):
    start_time = datetime.now()
    bwt_idx_name = bwt_idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Mapping reads against %s with Bowtie" % (start_time.strftime("%c"), bwt_idx_name)
    
    # Setup Bowtie output redirects
    #bwt_map = output_dir + mapped_reads
    bwt_map = tmp_name()
    tmp_fname = bwt_map.split('/')[-1]
    bwt_log = open(logging_dir + tmp_fname + ".log", "w")
    
    # Launch Bowtie
    try:    
        bowtie_cmd = ["bowtie"]
        
        if reads_format == "fastq":
            bowtie_cmd += ["-q"]
        elif reads_format == "fasta":
            bowtie_cmd += ["-f"]
            
        if unmapped_reads != None:
            unmapped_reads_fasta_name = unmapped_reads
            bowtie_cmd += ["--un", unmapped_reads_fasta_name,
                           "--max", "/dev/null"]
        else:
            unmapped_reads_fasta_name = None
        
        bowtie_cmd += ["--sam",
                       "-m 1",
                       "-y",
                       #"-k 1",
                       "-v 0",
                       #"-l 1000",
                       "--best",
                       "--strata",
                       #"-n", str(params.segment_mismatches),
                       #"-e", str(phred_thresh),
                       "-p", str(params.system_params.bowtie_threads),
                       bwt_idx_prefix, 
                       "-1", left_reads,
                       "-2", right_reads,
                       mapped_reads]
        
        #bowtie_proc = subprocess.Popen(bowtie_cmd, stderr=bwt_log)
             
        print >> run_log, " ".join(bowtie_cmd)
        
        
        print >> run_log, " ".join(bowtie_cmd)
        ret = subprocess.call(bowtie_cmd,
                              stderr=bwt_log)
                              # Bowtie reported an error
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute bowtie"
            exit(1)
            
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
            
    # Success    
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return (bwt_map, unmapped_reads_fasta_name)

def tmp_name():
    tmp_root = output_dir + "tmp/"
    if os.path.exists(tmp_root):
        pass
    else:        
        os.mkdir(tmp_root)
    return tmp_root + os.tmpnam().split('/')[-1] 

def build_target_bwt_index(target_prefix):
    print >> sys.stderr, "[%s] Indexing target sequences" % (right_now())
    bowtie_build_log = open(logging_dir + "index_target.log", "w")
    
    #user_splices_out_prefix  = output_dir + "user_splices_idx"
    
#    bowtie_build_cmd = ["bowtie-build", 
#                        "--version"] 
    
    bowtie_build_cmd = ["bowtie-build", 
                        target_prefix + ".fa",
                        target_prefix]            
    try:    
        print >> run_log, " ".join(bowtie_build_cmd)
        retcode = subprocess.call(bowtie_build_cmd, 
                                 stdout=bowtie_build_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Indexing failed with err =", retcode
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-build not found on this system"
        exit(1)
    return target_prefix

def write_sam_header(read_params, sam_file):
    print >> sam_file, "@HD\tVN:1.0\tSO:sorted"
    
    if read_params.read_group_id and read_params.sample_id:
        rg_str = "@RG\tID:%s\tSM:%s" % (read_params.read_group_id,
                                        read_params.sample_id)
        if read_params.library_id:
            rg_str += "\tLB:%s" % read_params.library_id
        if read_params.description:
            rg_str += "\tDS:%s" % read_params.description
        if read_params.seq_platform_unit:
            rg_str += "\tPU:%s" % read_params.seq_platform_unit
        if read_params.seq_center:
            rg_str += "\tCN:%s" % read_params.seq_center
        if read_params.mate_inner_dist:
            rg_str += "\tPI:%s" % read_params.mate_inner_dist
        if read_params.seq_run_date:
            rg_str += "\tDT:%s" % read_params.seq_run_date
        if read_params.seq_platform:
            rg_str += "\tPL:%s" % read_params.seq_platform
        
        print >> sam_file, rg_str
    print >> sam_file, "@PG\tID:Spats\tVN:%s\tCL:%s" % (get_version(), run_cmd)

def get_version():
   return "0.0.1"

# From http://www.dalkescientific.com/writings/NBN/parsing.html
class FastaRecord(object):
    def __init__(self, title, sequence):
        self.title = title
        self.sequence = sequence

def read_fasta_record(infile):
    # The first line in a FASTA record is the title line.
    # Examples:
    # >third sequence record
    # >gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene
    line = infile.readline()

    if not line:
        # End of file
        return None

    # Double-check that it's a FASTA file by looking for the '>'.
    if not line.startswith(">"):
        raise TypeError("Not a FASTA file: %r" % line)

    # The title starts after the '>' and doesn't include any
    # trailing whitespace.
    title = line[1:].rstrip()

    # Read the sequence lines up to the blank line.
    sequence_lines = []
    while 1:
        # I know I'm at the end of the sequence lines when I reach the
        # blank line or the end of the file.  I can simplify the test
        # by calling rstrip and checking for the empty string.  If it's
        # a sequence line I'll still call rstrip to remove the final
        # newline and any trailing whitespace so I'll go ahead and
        # always rstring the line.
        line = infile.readline().rstrip()
        if line == "":
            # Reached the end of the record or end of the file
            break
        sequence_lines.append(line)

    # Merge the lines together to make a single string containing the sequence
    # (This is faster than doing "sequence = sequence + line" if there are
    # more than a few lines)
    sequence = "".join(sequence_lines)
    
    return FastaRecord(title, sequence)

def read_fasta_records(input_file):
    records = []
    while 1:
        record = read_fasta_record(input_file)
        if record is None:
            break
        records.append(record)
    return records

###########################################
# End of pilfered code
   
def print_fasta(fout, name, seq):
    print >> fout, ">"+name
    i = 0
    while i < len(seq):
        print >> fout,seq[i:i+60]
        i += 60
    print >> fout, seq[i:]
    print >> fout

def index_targets(target_file):
    targets = []
    target_fasta_name = tmp_dir + "targets.fa"
    target_fasta = open(target_fasta_name, "w")
    
    for record in read_fasta_records(open(target_file,"r")):
        tr_name = record.title
        tr_seq = record.sequence
        print_fasta(target_fasta, tr_name, tr_seq)

    target_fasta.close()
    return build_target_bwt_index(tmp_dir + "targets")
  
def compute_profiles(params, target_fasta, treated_alignments, untreated_alignments):
    #filter_cmd = ["prep_reads"]
    print >> sys.stderr, "[%s] Building reactivity profiles" % (right_now())

    #filter_log = open(logging_dir + "relabel_reads.log", "w")
    
    cmd = [bin_dir + "compute_profiles"]
    cmd.extend(params.cmd())
    cmd.extend([target_fasta, treated_alignments, untreated_alignments])
       
    # print "\t executing: `%s'" % " ".join(cmd)    
    # files = reads_list.split(',')
    # for reads_file in files:
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
                              # Bowtie reported an error
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute compute_profiles"
            exit(1)
    # prep_reads not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: compute_profiles not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def reverse_complement(s):
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
        

def main(argv=None):
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
    
    # Initialize default parameter values
    params = SpatsParams()
    
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()
        
        rna_targets_filename = args[0]
        treated_handle_seq = args[1]
        untreated_handle_seq = args[2]
        left_reads_list = args[3]
        right_reads_list = args[4]
            
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning Spats run (v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------" 
        
        start_time = datetime.now()
        prepare_output_dir()
        
        global run_log
        run_log = open(logging_dir + "run.log", "w", 0)
        global run_cmd
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd
        
        # Validate all the input files, check all prereqs before committing 
        # to the run
        #(ref_fasta, ref_seq_dict) = check_index(bwt_idx_prefix)
        
        check_bowtie()
        
        if len(left_reads_list) != len(right_reads_list):
            print >> sys.stderr, "Error: reads must be supplied as matched pairs of files"
            sys.exit(1)
        
        left_kept_reads_list = []
        right_kept_reads_list = []
        
        left_reads_list_filenames = left_reads_list.split(',')
        right_reads_list_filenames = right_reads_list.split(',')
        
        for i in range(0, len(left_reads_list_filenames)):
            left_trimmed_reads = left_reads_list_filenames[i] + ".trimmed"
            right_trimmed_reads = right_reads_list_filenames[i] + ".trimmed"
                    
            left_kept_reads = left_reads_list_filenames[i] + ".kept"
            right_kept_reads = right_reads_list_filenames[i] + ".kept"
            
            left_trimmed_reads = trim_read_adapters(params, 
                                                    params.read_params.left_adapter,
                                                    params.read_params.right_adapter,
                                                    left_reads_list_filenames[i],
                                                    left_trimmed_reads)
            right_trimmed_reads = trim_read_adapters(params, 
                                                     reverse_complement(params.read_params.left_adapter),
                                                     reverse_complement(params.read_params.right_adapter),
                                                     right_reads_list_filenames[i],
                                                     right_trimmed_reads)
            [left_kept_reads, right_kept_reads] = match_read_pairs(params,
                                                                   left_trimmed_reads, 
                                                                   right_trimmed_reads, 
                                                                   left_kept_reads, 
                                                                   right_kept_reads)
            left_kept_reads_list.append(left_kept_reads)
            right_kept_reads_list.append(right_kept_reads)
                                
        # Now start the time consuming stuff
        relabel_reads(params,
                      left_reads_list, 
                      right_reads_list,
                      treated_handle_seq,
                      untreated_handle_seq)
                      
        index_prefix = index_targets(rna_targets_filename)
        
        maps = []
        
        for handle_seq in [treated_handle_seq, untreated_handle_seq]:
            left_labeled_reads = output_dir + "/" + handle_seq + "_1.fq"
            right_labeled_reads = output_dir + "/" + handle_seq + "_2.fq"

            mapped_reads = output_dir + handle_seq + ".sam"                            
            bowtie(params,
                   index_prefix,
                   left_labeled_reads,
                   right_labeled_reads,
                   "fastq",
                   mapped_reads,
                   None,
                   500)
            maps.append(mapped_reads)
                   
        treated_map = maps[0]
        untreated_map = maps[1]
        
        compute_profiles(params, index_prefix + ".fa", treated_map, untreated_map)
               
        if params.system_params.keep_tmp == False:
            tmp_files = os.listdir(tmp_dir)
            for t in tmp_files:
                os.remove(tmp_dir+t)
            os.rmdir(tmp_dir)
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Run complete [%s elapsed]" %  formatTD(duration)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "    for detailed help see http://spats.cbcb.umd.edu/manual.html"
        return 2


if __name__ == "__main__":
    sys.exit(main())

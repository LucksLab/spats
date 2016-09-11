#!/usr/bin/env python
# encoding: utf-8

#python script to take SHAPE-Seq reactivities file and split into individual files

"""
reactivities_split.py

Created by Julius Lucks on 2010-10-20.
Copyright (c) 2012 Julius Lucks. All rights reserved.
"""

import sys
import getopt
import fileinput
import re
import os
from datetime import datetime, date, time

use_message = '''
 reactivities_split takes a spats reactivities output file and divides it into individual reactivity files for each RNA, containing all of the 
 RT start positions (unsplit).  These reactivity files are then split into smaller files that contain only one RT start position each.
 
 
 Usage:
     reactivities_split [options] <reactivities_file.out>
 
 Options: 
     
 -h               Prints this help menus
 -v               Prints the version number
 -q,--quiet       Suppresses output of reactivity filenames being generated (useful for a lot of RNA targets)
 -g,--group       Groups all split reactivities files into one directory, regardless of RT start position (useful for random priming)
 -l,--long-only   Only generates reactivities files for the longest RT start position of any target (most 3'), takes precedence over -g or --group

'''

def get_version():
   return "0.0.3"

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class Params:

    def __init__(self):
        pass
        
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvqgl", 
                                        ["version",
                                         "help",
                                         "quiet",
                                         "group",
                                         "long-only",
                                         "option1="])
        except getopt.error, msg:
            raise Usage(msg)
        
        silence = False
        group = False
        long_only = False
        for option, value in opts:
            if option in ("-v", "--version"):
                print "reactivities_split v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option in ("-q", "--quiet"):
                silence = True
            if option in ("-g", "--group"):
                group = True
            if option in ("-l", "--long-only"):
                long_only = True
            if option == "--option1":
                self.option1 = value
        
        if len(args) < 1:
            raise Usage(use_message)
        return args, silence, group,long_only

    def check(self):
        pass


def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def print_file(name,print_buffer):
    output = file(name+'_reactivities.txt','w')
    output.write(print_buffer)
    output.close()

def main(argv=None):

    # Initialize default parameter values
    params = Params()

    try:
        if argv is None:
            argv = sys.argv
            args,silence,group,long_only = params.parse_options(argv)
            params.check()

        reactivities_file = args[0]

        print >> sys.stderr
        print >> sys.stderr, "[%s] Splitting %s into Individual RNA Files:" % (right_now(), reactivities_file)
        print >> sys.stderr, "---------------------------------------------------------------------------------"
        
        #Block to make directories for the split reactivities files
        unsplit_reactivities_dir="Reactivities_all_RT_starts_unsplit/"
        if not os.path.exists(unsplit_reactivities_dir):
            os.mkdir(unsplit_reactivities_dir[:-1])
        if not group or long_only:
            long_RT_dir="Reactivities_Longest_RT_starts_only_split/"
            if not os.path.exists(long_RT_dir[:-1]):
                os.mkdir(long_RT_dir)
            truncated_RT_dir="Reactivities_Internal_RT_starts_split/"
            if not os.path.exists(truncated_RT_dir) and not long_only:
                os.mkdir(truncated_RT_dir[:-1])
        else:
            long_RT_dir="Reactivities_all_RT_starts_split/"
            truncated_RT_dir=long_RT_dir
            if not os.path.exists(long_RT_dir[:-1]):
                os.mkdir(long_RT_dir)
            if not os.path.exists(truncated_RT_dir):
                os.mkdir(truncated_RT_dir[:-1])
            
        files_gen = []
        first_line = None
        print_buff = ''
        name = ''
        new_name = ''
        for line in fileinput.input(reactivities_file):
            if first_line is None:
                first_line = line
            else:
                parse = re.split('\s+',line)
                new_name = parse[0]
                if not name == new_name:
                    files_gen.append(unsplit_reactivities_dir+new_name+'_reactivities.txt')
                    if not silence:
                        print new_name
                    #new RNA - print out buffer and reset everything
                    if print_buff is not '':
                        #Only print if have something
                        print_file(unsplit_reactivities_dir+name,first_line+print_buff)
                        print_buff = line #reset print_buff to be this line
                        name = new_name
                    else:
                        print_buff = line #don't lose this line
                else:
                    #store this line in the buffer
                    print_buff += line
                name = new_name
        
        #print out the very last one
        print_file(unsplit_reactivities_dir+name,first_line+print_buff)
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Splitting RT Start Positions in RNA Files " % (right_now())
        print >> sys.stderr, "---------------------------------------------------------------------"
        
        for filename in files_gen:
            RT_start = None
            first_line = None
            name = ''
            print_buff = ''
            for line in fileinput.input(filename):
                if first_line is None:
                    first_line = line
                else:
                    parse = re.split('\s+',line)
                    if name == '':
                        name = parse[0]
                    new_RT_start = parse[1]    
                    if not RT_start == new_RT_start:
                        if not silence:
                            print name + "_" + str(new_RT_start)
                        #new RNA - print out buffer and reset everything
                        if print_buff is not '':
                            #Only print if have something
                            if not long_only:
                                print_file(truncated_RT_dir+name+"_"+str(RT_start),first_line+print_buff)
                            print_buff = line #reset print_buff to be this line
                            RT_start = new_RT_start
                        else:
                            print_buff = line #don't lose this line
                    else:
                        #store this line in the buffer
                        print_buff += line
                    RT_start = new_RT_start
                      
            #print out the last RT length of the last
            print_file(long_RT_dir+ name+"_"+str(RT_start),first_line+print_buff)
        
        print >> sys.stderr       
         
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "    for detailed help see http://spats.cbcb.umd.edu/manual.html"
        return 2


if __name__ == "__main__":
    sys.exit(main())
        
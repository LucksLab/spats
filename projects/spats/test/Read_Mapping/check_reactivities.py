from itertools import product, cycle
import getopt
import sys

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

opts = getopts("", ["input=", "permutation=", "adapter=", "sequence="])
spats_reactivities_out = opts["--input"]
permutation = int(opts["--permutation"])
adapter_len = len(opts["--adapter"])
sequence_len = len(opts["--sequence"])
print "Checking permutation %s..."%(permutation)

all_permutations = [seq for seq in product([0,1], repeat=4)]
permutation_cycle = cycle(all_permutations[permutation])
print all_permutations[permutation]

reads = []
with open(spats_reactivities_out, "r") as f:
    header = f.readline()
    for line in f:
        fields = line.split("\t")
        reads += [int(fields[4]), int(fields[5])]

correct = sum([1 if a[0] == a[1] else 0 for a in zip(permutation_cycle, reads[:2*(sequence_len+1)])])
correct += sum([1 if a == 0 else 0 for a in reads[2*(sequence_len+1):]])
incorrect = len(reads) - correct
exp_read_lines = 2 * (sequence_len + adapter_len + 1)

if correct == len(reads) and len(reads) == exp_read_lines:
    print "Permutation %s: OK - %s reads out of %s expected, %s correct"%(permutation, len(reads), exp_read_lines, correct)
else:
    print "Permutation %s: FAILED - %s reads out of %s expected, %s incorrect, %s correct"%(permutation, len(reads), exp_read_lines, incorrect, correct)

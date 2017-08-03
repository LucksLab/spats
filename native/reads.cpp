
#include "db.hpp"
#include "target.hpp"


int
main(int argc, char ** argv)
{
    if (5 != argc  &&  6 != argc) {
        fprintf(stderr, "Usage: reads [target_path] [r1_path] [r2_path] [spats_out_path] [optional: sample_size = 100000]\n");
        return -1;
    }
    PairDB * pdb = new PairDB(argv[4]);
    Targets * t = new Targets();
    t->parse(argv[1]);
    pdb->store_targets(t);
    int sample_size = (6 == argc ? atoi(argv[5]) : 100000);
    pdb->parse_and_sample(argv[2], argv[3], sample_size);
    pdb->close();
    return 0;
}


#include <stdio.h>
#include "spats.hpp"

int
main(int argc, char ** argv)
{
    if (5 != argc) {
        fprintf(stderr, "Usage: cotrans [target_path] [r1_path] [r2_path] [spats_out_path]\n");
        return -1;
    }
    Spats s(true);
    s.addTargets(argv[1]);
    s.process_pair_data(argv[2], argv[3]);
    s.store(argv[4]);
    return 0;
}

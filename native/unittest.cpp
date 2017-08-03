
#import "spats.hpp"

int
main(int argc, char ** argv)
{
    if (4 != argc) {
        fprintf(stderr, "Usage: unittest [target_path] [r1_seq] [r2_seq]\n");
        return -1;
    }
    Spats s(true);
    s.addTargets(argv[1]);
    Case c("", argv[2], argv[3]);
    s.run_case(&c);
    if (MASK_NO_MATCH == c.mask  ||  c.L <= 0  ||  c.site < 0) {
        printf("[ None, None, None ]\n");
        return 1;
    }
    else {
        printf("[ \"%s\", %d, %d ]\n", (MASK_TREATED == c.mask ? "RRRY" : "YYYR"), c.L, c.site);
        return 0;
    }
}

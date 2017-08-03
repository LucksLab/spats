
#import "spats.hpp"

int
main(int argc, char ** argv)
{
    Spats s(true);
    s.addTargets("/Users/jbrink/mos/tasks/1RwIBa/tmp/datasets/cotrans/cotrans_single.fa");
    Case c("", argv[1], argv[2]);
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

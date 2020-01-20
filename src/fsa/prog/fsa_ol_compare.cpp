#include "../overlap_compare.hpp"


int main(int argc, char *argv[]) {
    OverlapCompare oc;

    if (oc.ParseArgument(argc, argv)) {
        oc.Run();
    }
    else {
        oc.Usage();
    }
    return 0;

}

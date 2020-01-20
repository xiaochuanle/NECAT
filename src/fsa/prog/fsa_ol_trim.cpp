#include "../overlap_trim.hpp"

int main(int argc, char *argv[]) {
    fsa::OverlapTrim ot;

    if (ot.ParseArgument(argc, argv)) {
        ot.Run();
    }
    else {
        ot.Usage();
    }
    return 0;

}

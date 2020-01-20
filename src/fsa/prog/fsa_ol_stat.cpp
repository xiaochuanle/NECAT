#include "../overlap_stat.hpp"

int main(int argc, char *argv[]) {
    fsa::OverlapStat os;

    if (os.ParseArgument(argc, argv)) {
        os.Run();
    }
    else {
        os.Usage();
    }
    return 0;

}

#include "../contig_bridge.hpp"


int main(int argc, const char* argv[]) {
    fsa::ContigBridge cb;

    if (cb.ParseArgument(argc, argv)) {
        cb.Run();
    }
    else {
        cb.Usage();
    }
    return 0;

}

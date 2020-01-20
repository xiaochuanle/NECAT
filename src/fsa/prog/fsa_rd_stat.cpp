#include "../read_stat.hpp"

int main(int argc, char *argv[]) {
    fsa::ReadStat rs;

    if (rs.ParseArgument(argc, argv)) {
        rs.Run();
    }
    else {
        rs.Usage();
    }
    return 0;

}

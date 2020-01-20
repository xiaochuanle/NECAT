#include "../read_extract.hpp"

int main(int argc, char *argv[]) {
    fsa::ReadExtract re;

    if (re.ParseArgument(argc, argv)) {
        re.Run();
    }
    else {
        re.Usage();
    }
    return 0;

}

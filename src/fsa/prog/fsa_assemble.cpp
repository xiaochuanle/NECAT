#include "../assembly.hpp"

int main(int argc, char *argv[])
{
    fsa::Assembly ass;
    if (ass.ParseArgument(argc, argv)) {
        ass.Run();
    }
    else {
        ass.Usage();
    }

	return 0;

}

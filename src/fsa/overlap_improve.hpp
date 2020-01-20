#ifndef FSA_OVERLAP_IMPROVER_HPP
#define FSA_OVERLAP_IMPROVER_HPP



#include "argument_parser.hpp"

class OverlapStore;
class Overlap;
class ReadStore;

class OverlapImprover {
public:
    bool ParseArgument(int argc, const char *const argv[]);
    void Usage();
    void Run();

    void Improve(OverlapStore &ol_store, ReadStore &read_store);

protected:

    bool Filter(Overlap &o);
    ArgumentParser GetArgumentParser();
protected:

    int max_overhang_{ 10 };
    int thread_size_{ 1 };
    std::string read_file_;
    std::string ifname_;
    std::string ofname_;

};



#endif // OVERLAP_IMPROVER_HPP
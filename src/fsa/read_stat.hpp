#ifndef FSA_READ_STAT_HPP
#define FSA_READ_STAT_HPP


#include <string>
#include <unordered_map>

#include "sequence.hpp"
#include "overlap_store.hpp"
#include "argument_parser.hpp"


namespace fsa {
class ReadStore;
class ReadStat {
public:
    bool ParseArgument(int argc, const char* const argv[]);
    void Run();
    void Usage();
protected:
    ArgumentParser GetArgumentParser();

    void StatReads();
    void StatN50(const ReadStore &rs);
protected:
    std::string action_ { "N50" };
    std::string ifname_;
    int thread_size_ { 4 };
};


} // namespace fsa {
#endif // FSA_READ_STAT_HPP

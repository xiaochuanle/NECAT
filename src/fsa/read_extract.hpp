#ifndef FSA_READ_EXTRACT_HPP
#define FSA_READ_EXTRACT_HPP

#include "argument_parser.hpp"
#include "read_store.hpp"
#include "logger.hpp"
#include <algorithm>

namespace fsa {
class ReadExtract {
public:
    bool ParseArgument(int argc, const char* const argv[]);
    void Run();
    void Usage();
protected:
    ArgumentParser GetArgumentParser();
    void PrintArguments();

    void StatReads();
    void StatN50(const ReadStore &rs);
    void ExtractLongest(const std::string &ifname, const std::string &ofname, long long size);
    void Split(const std::string &ifname, const std::string &opattern, long long size);
    void Copy(const std::string &ifname, const std::string &ofname);
    size_t FindLongestXHeap(std::vector<int> &lengths, long long goal);
protected:
    std::string action_ { "longest" };
    long long base_size_ {0 };
    std::string ifname_;
    std::string ofname_;
    int thread_size_ { 4 };
};

} // namespace fsa {

#endif // FSA_READ_EXTRACT_HPP
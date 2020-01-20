#ifndef FSA_OVERLAP_STAT_HPP
#define FSA_OVERLAP_STAT_HPP


#include <string>
#include <unordered_map>

#include "sequence.hpp"
#include "overlap_store.hpp"
#include "argument_parser.hpp"

namespace fsa {

class OverlapStat {
public:
    bool ParseArgument(int argc, const char* const argv[]);
    void Run();
    void Usage();
protected:
    ArgumentParser GetArgumentParser();
    
public:
    struct Read {
        int len;
    };

    void Load(const std::string &fname, const std::string &type="");
    void StatReads();
    void FindCommon(const std::string &fname, const std::string &read0, const std::string &read1);
protected:
    std::string ifname_;
    std::string read0_;
    std::string read1_;
    int thread_size_;
 
    OverlapStore ol_store_;
    std::unordered_map<Seq::Id, Read> reads_;
};


} // namespace fsa {

#endif // FSA_OVERLAP_STAT_HPP

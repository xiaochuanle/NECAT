#ifndef FSA_READ_TOOLS_HPP
#define FSA_READ_TOOLS_HPP

#include <string>
#include <unordered_set>

#include "program.hpp"
#include "utility.hpp"

namespace fsa {

class ReadTools : public Program {
public:
    
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();

    void StatN50(const std::string& ifname);
    void Split(const std::string &ifname, const std::string &opattern, long long block_size, long long base_size);
    void SplitName(const std::string &ifname, const std::string &opattern, long long block_size, long long base_size);
    void Longest(const std::string &ifname, const std::string &ofname, long long base_size, int min_length);
    void Check(const std::string &ifname);
protected:
    std::unordered_set<std::array<std::string, 2>, ArrayHash<std::string, 2>, ArrayCompare<std::string,2>> LoadIgnored(const std::string &fname);
    std::unordered_map<std::string, std::string> LoadNamemap(const std::string &fname);
protected:
    std::string cmd_;
    std::string namemap_;
    std::string ifname_;
    std::string ofname_;
    std::string names_;
    std::string id2name_;
    bool discard_illegal_read_ { false };
    int thread_size_ { 1 };
    int size_ { 0 };
    long long block_size_ { 0 };
    long long base_size_ { 0 };
    long long genome_size_ { 0 };
    int min_length_ { 0 }; 
};



} // namespace fsa

#endif // FSA_READ_TOOLS_HPP

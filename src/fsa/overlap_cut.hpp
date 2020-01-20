#ifndef FSA_OVERLAP_CUT_HPP
#define FSA_OVERLAP_CUT_HPP


#include <string>
#include <unordered_map>

#include "sequence.hpp"
#include "overlap_store.hpp"
#include "argument_parser.hpp"


namespace fsa {

class OverlapCut {
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

    struct IterCigar {
        IterCigar(const std::string &c) : cigar(c) {}
        bool Next(size_t &n, char &tpye);
        const std::string &cigar;
        size_t index {0};
    };

    void Load(const std::string &fname, const std::string &type="");
    void StatReads();
    void FindCommon(const std::string &fname, const std::string &read0, const std::string &read1);
    void ComputeIdentity(const std::string& qseq, size_t qs, size_t qe, const std::string &tseq, size_t ts, size_t te, const std::string &cigar);
    std::array<size_t, 2> ComputeIdentity(size_t qs, size_t qe, size_t ts, size_t te, const std::string &cigar);
    std::string CutPafLine(const std::string &line, std::ofstream &of);
protected:
    std::string ifname_;
    std::string ofname_;
    std::string read_file_;
    int window_size_ { 1000 };
    int thread_size_ { 1 };
    int min_aligned_length_ { 2000 };

    std::string local_identity_fname_ ;
 
    OverlapStore ol_store_;
    ReadStore rd_store_;
    std::unordered_map<Seq::Id, Read> reads_;
};

} // namespace fsa {

#endif // FSA_OVERLAP_CUT_HPP


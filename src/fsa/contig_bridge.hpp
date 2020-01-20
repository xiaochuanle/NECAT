#ifndef FSA_CONTIG_BRIDGE_HPP
#define FSA_CONTIG_BRIDGE_HPP

#include "argument_parser.hpp"
#include "contig_link_store.hpp"
#include "contig_graph.hpp"
#include "read_store.hpp"

namespace fsa {

class ContigBridge {

public:

    struct ReadStatInfo {
        double identity{0.0};
        int overhang {0}; 
        int score{0};
        int len {-1};
        int aligned { -1 };
        int count {0};  
        int oh_count {0};
    };

    ArgumentParser GetArgumentParser();
    void Usage() ;
    bool ParseArgument(int argc, const char *const argv[]);
    void Run();

    std::string OutputPath(const std::string &name) {
        return output_directory_ + "//" + name;
    }

    void SaveBridgedContigs(const std::string &fname);
    void PrintArguments() const;

    void AutoSelectCtg2ctgParams();
    void AutoSelectCtg2ctgMinIdentity(const std::unordered_map<Seq::Id, ReadStatInfo> &info);
    void AutoSelectCtg2ctgMaxOverhang(const std::unordered_map<Seq::Id, ReadStatInfo> &info);
    void AutoSelectRead2ctgParams();
    void AutoSelectRead2ctgMinIdentity(const std::unordered_map<Seq::Id, ReadStatInfo> &info);
    void AutoSelectRead2ctgMaxOverhang(const std::unordered_map<Seq::Id, ReadStatInfo> &info);

    void SelectCtgParamtersByReadInfos(const std::string& fname);
    
    std::unordered_map<Seq::Id, ReadStatInfo> StatReadInfo(const std::string &fname, int th_identity, int th_overhang);
    void Dump() ;
protected:
    std::string read_file_{ "" };
    std::string contig_file_{ "" };
    std::string read2ctg_file_{ "" };
    std::string bridged_contig_file_{ "" };
    std::string ctg2ctg_file_{ "" };

    std::string select_branch_{ "one" };
    std::string output_directory_{ "." };

    bool dump_{ true };
    double read2ctg_min_identity_{ 80 };
    double ctg2ctg_min_identity_{ 95 };
    int read_min_length_{ 5000 };
    int ctg_min_length_{ 500 };
    int read2ctg_max_overhang_{ 500 };
    int ctg2ctg_max_overhang_ { 100 };
    int read2ctg_min_aligned_length_{ 5000 };
    int ctg2ctg_min_aligned_length_{ 2000 };
    int read2ctg_min_coverage_{ 3 };
    int min_contig_length_ { 500 };
    int thread_size_ { 4 };
    int ctg2ctg_max_overhang_raw_ { 300 };
    double ctg2ctg_min_identity_raw_ { 90 };
    int read2ctg_max_overhang_raw_ { 1000 };
    double read2ctg_min_identity_raw_ { 70.0 };
    int window_size_ {1000};
    std::string readinfo_fname_ { "" };

    ReadStore read_store_;
    ContigLinkStore contig_links_{ read_store_ };
    ContigGraph contig_graph_{ contig_links_ };
    std::unordered_set<Seq::Id> contigs_;
    
    ArgumentParser ap_ {GetArgumentParser()};
};

} // namespace fsa {
#endif // FSA_CONTIG_BRIDGE_HPP  

#ifndef FSA_ASSEMBLE_HPP
#define FSA_ASSEMBLE_HPP

#include <vector>
#include <map>

#include "argument_parser.hpp"
#include "string_graph.hpp"
#include "path_graph.hpp"
#include "overlap_store.hpp"
#include "read_store.hpp"

namespace fsa {

struct Options {
    // The overlaps checked by fsa_ol_filter are considered to be of good-quality.
    // The default values of min_length, min_aligned_length and min_identity are set to 0.
    int min_length{ 0 };                 
    int min_aligned_length{ 0 }; 
    double min_identity{ 0 };       
    int min_contig_length { 500 };
    int max_bubble_identity{ 96 };
    int max_bubble_coverage{ 97 };
    bool lfc{ false };
    bool remove_chimer{ false };
    bool help{ false };
    std::string overlap_file{ "" };       
    std::string read_file{ "" };         
    std::string output_directory{ "." };
    std::string select_branch{ "no" };
    int run_mode{ 4 };                    
    int dump{ 0 };
    std::string overlap_file_type{ "" };
    int thread_size {1};
    int max_spur_length { 50000 };
};

class Assembly {
public:
    Assembly();

    bool ParseArgument(int argc, char *const argv[]);
    void Run();
    void Usage();

    ArgumentParser GetArgumentParser();

    void LoadOverlaps(const std::string &fname);
    void LoadReads(const std::string &fname);
    std::unordered_set<Seq::Id> CollectReadIdsInContigs();
    void CreateStringGraph();
    void CreatePathGraph();
    void SaveGraph();
    void SaveContigs();

    void SavePContig(FILE* file, int id, const std::list<StringEdge*> &pcontig);
    void SavePContig1(FILE* file, int id, const std::list<StringEdge*> &pcontig);

    void SavePContigTiles(FILE* file, int id, const std::list<StringEdge*> &pcontig);

    void SaveContigs(FILE *fseq, FILE *ftile, int id, const std::list<StringEdge*> &contigs);
    void SaveBubbles(FILE *fseq, FILE* ftile, int ctgid, const std::list<std::pair<CompoundPathEdge*, std::list<std::list<StringEdge*>>>> &bubbles);

protected:
    std::string ConstructContigStraight(const std::list<StringEdge*> &contig);
    std::string ConstructContig(const std::list<StringEdge*> &contig);
    std::string ConstructContigMain(const std::list<StringEdge*> &contig);
    std::vector<std::string> ConstructContig1(const std::list<StringEdge*> &contig);
    std::vector<std::string> ConstructContigAll(const std::list<StringEdge*> &contig);
    std::array<double,2> ComputeSequenceSimilarity(const std::string &qseq, const std::string &tseq);
    std::string EdgeToSeq(const StringEdge *e);
    std::string OutputPath(const std::string &fname) { return options_.output_directory + "/" + fname; }
    void PrintArguments();
protected:
    Options options_;
    ReadStore read_store_;
    OverlapStore ol_store_;
    StringGraph string_graph_;
    PathGraph path_graph_;
};
} // namespace fsa {

#endif // FSA_ASSEMBLE_HPP

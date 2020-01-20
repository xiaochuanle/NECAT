#ifndef FSA_OVERLAP_COMPARE_HPP
#define FSA_OVERLAP_COMPARE_HPP


#include <string>
#include "overlap_store.hpp"
#include "argument_parser.hpp"

namespace fsa {

class OverlapCompare {
public:
    OverlapCompare();
    ~OverlapCompare();
    bool ParseArgument(int argc, const char *const argv[]);
    void Usage();
    void Run();
    void PrintArguments();
protected:
    struct Stats {
        Stats(const OverlapStore& ols, int ee);
        uint64_t OverlapId(const Overlap &o) {
            return ((uint64_t)std::min(o.a_.id, o.b_.id) << 32) + std::max(o.a_.id, o.b_.id);
        }

        bool IsContain(const Overlap &o);
        bool IsNone(const Overlap &o);
        bool IsBetter(const Overlap &a, const Overlap &b);

        const OverlapStore &overlaps;

        std::unordered_map<uint64_t, std::vector<const Overlap*>> map;
        std::unordered_set<const Overlap*> contains;
        std::unordered_set<const Overlap*> nones;
        std::unordered_set<const Overlap*> dups;
        std::unordered_map<uint64_t, const Overlap*> goods;

        int max_overhang;
    };

    ArgumentParser GetArgumentParser();
    void Compare(const Stats &a, const Stats &b, const std::string& name);
    bool IsSimilar(const Overlap &a, const Overlap &b);

    void Load();
protected:
    int thread_size_ { 4 };
    int max_overhang_{ 250 };
    int min_length_ { 2500 };
    int min_aligned_length_ { 2000 };
    std::array<std::string, 2> fnames_;
    std::string output_directory_;
    
    ReadStore read_store_;
    std::array<OverlapStore*, 2> ol_stores_;


};

} // namespace fsa {

#endif // FSA_OVERLAP_COMPARE_HPP

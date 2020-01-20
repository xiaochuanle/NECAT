#ifndef FSA_OVERLAP_TRIM_HPP
#define FSA_OVERLAP_TRIM_HPP

#include <string>

#include "argument_parser.hpp"
#include "overlap_store.hpp"

namespace fsa {

class OverlapTrim {
public:

    bool ParseArgument(int argc, const char *const argv[]);

    void Run();
    void Usage();

    
protected:
    ArgumentParser GetArgumentParser();

    void LoadOverlaps(const std::string &fname);
    void Trim();
    void GroupAndFilterDuplicate();
    
    void FilterCoverage();
    void FilterCoverageMt();
    void FilterCoverage(int min_coverage, int max_coverage, int max_diff_coverge);

    /** Calc coverage param: min_coverage, max_coverage and max_diff_coverage */
    std::array<int, 3> CoverageParam() const;
    std::array<int, 3> CoverageParam1() const;
    std::pair<int, int> CalcMinMaxCoverage(int id, const std::unordered_map<int, const Overlap*>& group);
    int FirstTrough(const std::vector<int> &data, size_t last, size_t k) const;
    void DumpCoverage(const std::string &fname) const;
protected:
    std::string overlap_fname_;
    int thread_size_ { 4 };
    int min_coverage_ { 3 };
    double min_identity_ { 80 };
    OverlapStore ol_store_;
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> groups_;
    std::unordered_map<Seq::Id, std::array<int, 2>> coverages_;                         //!< record min and max base coverages of the reads

    
};

} // namespace fsa {

#endif // FSA_OVERLAP_TRIM_HPP

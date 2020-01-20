#ifndef FSA_OVERLAP_FILTER_HPP
#define FSA_OVERLAP_FILTER_HPP


#include<climits>
#include <sstream>
#include <array>

#include "overlap_store.hpp"
#include "argument_parser.hpp"

namespace fsa {
class OverlapFilter {
public:
    OverlapFilter();
    virtual ~OverlapFilter();

    bool ParseArgument(int argc, char *const argv[]);

    void Run();
    void Usage();

protected:
    ArgumentParser GetArgumentParser();

protected:
    struct RdReason {
        enum Type {
            RS_OK = 0,
            RS_CONTAINED,
            RS_COVERAGE,
            RS_NO_LONGEST,
            RS_UNKNOWN
        };
        static RdReason Contained(int id) { return RdReason(Type::RS_CONTAINED, id); }
        static RdReason Coverage(int low, int high) { return RdReason(Type::RS_COVERAGE, low, high); }
        static RdReason Coverage(const std::array<int, 2>& a) { return RdReason(Type::RS_COVERAGE, a[0], a[1]); }
        static RdReason NoLongest() { return RdReason(Type::RS_NO_LONGEST); }
        Type type;
        std::array<int, 2> sub;

protected:
        RdReason(Type t=RS_UNKNOWN, int s0=0, int s1=0) { type = t; sub[0] = s0;  sub[1] = s1;}
    };

    struct OlReason {
        enum Type {
            RS_OK = 0,
            RS_SIMPLE,
            RS_DUPLICATE,
            RS_BESTN,
            RS_FILTERED_READ,
            RS_LACK_OF_SUPPORT,
            RS_LOCAL,
            RS_UNKNOWN
        };
        static OlReason Ok() { return OlReason(RS_OK); }
        static OlReason Simple() { return OlReason(RS_SIMPLE); }
        static OlReason Duplicate() { return OlReason(RS_DUPLICATE); }
        static OlReason BestN() { return OlReason(RS_BESTN); }
        static OlReason FilteredRead(int id) { return OlReason(RS_FILTERED_READ, id); }
        static OlReason LackOfSupport() { return OlReason(RS_LACK_OF_SUPPORT); }
        static OlReason Local(int t, int v) { return OlReason(RS_LOCAL, t*10000+v); }
        
        OlReason(Type t=RS_UNKNOWN, int p0=0, int p1=0) : type(t), sub{p0, p1} { }

        Type type;
        std::array<int, 2> sub;
    };

    struct ReadStatInfo {
        double identity{0.0};
        double identity_threshold {0.0};
        double overhang_r_threshold {0.0};
        double overhang_l_threshold {0.0};
        int overhang {0}; 
        int len {-1};
        size_t count {0};
        double Identity() { return identity; }
        int CompareIdentity(const ReadStatInfo& a) const {
            if (identity < a.identity) {
                return -1;
            } else if (identity == a.identity) {
                return 0;
            } else {
                return 1;
            }
        }
        int CompareOverhang(const ReadStatInfo& a) const {
            if (overhang < a.overhang) {
                return -1;
            } else if (overhang == a.overhang) {
                return 0;
            } else {
                return 1;
            }
        }
    };

    void LoadOverlaps(const std::string &fname);
    void LoadOverlapsWithoutLowQuality(const std::string &fname);
    void SaveOverlaps(const std::string &fname);

    void StatLowQuality();
    void FilterLowQuality();
    void FilterLowQuality(const Overlap &o);
    void GroupAndFilterDuplicate();
    void GroupAndFilterDuplicateMt();
    void FilterContained();
    void FilterContainedMt();
    void FilterCoverage();
    void FilterCoverageMt();
    void FilterLackOfSupport();
    void FilterLackOfSupportMt();
    void CalcLocalThreshold(Seq::Id id, const std::unordered_map<Seq::Id, const Overlap*>& g);
    void FilterBestN();
    void FilterBestNMt();

    void AutoSelectParams();
    void AutoSelectMinLength(std::vector<int> &lengths);
    double CalcMinIdentity(const std::unordered_map<Seq::Id, ReadStatInfo>& readinfos, double dev1, double dev2, double minvalue);
    int CalcMaxOverhang(const std::unordered_map<Seq::Id, ReadStatInfo>& readinfos, double dev1, double dev2, int maxvalue);

    std::array<double,2> StatIdentity(const std::unordered_map<Seq::Id, ReadStatInfo> &info, int n);
    std::array<double,2> StatOverhang(const std::unordered_map<Seq::Id, ReadStatInfo> &info, int n);

    double CalcLocalIdentityThreshold(std::vector<std::array<double,2>> &idents, int base_lower_limit, int base_upper_limit);
    double CalcLocalOverhangThreshold(std::vector<std::array<double,2>> &overhang, int lower_limit, int upper_limit);
    ReadStatInfo CalcReadInfo(Seq::Id id, const std::unordered_map<Seq::Id, const Overlap*>& group);

    void ModifyEnd(const Overlap &o, int maxoh);

    std::pair<int, int> CalcMinMaxCoverage(int id, const std::unordered_map<int, const Overlap*>& group);

    std::string OutputPath(const std::string &fname) const { return output_directory_+"/"+fname; }
    void PrintArguments();

    /** Calc coverage param: min_coverage, max_coverage and max_diff_coverage */
    std::array<int, 3> CoverageParam() const;
    std::array<int, 3> CoverageParam1() const;

    size_t FindLongestXHeap(std::vector<std::array<int,2>> &lengths, long long goal);
    size_t FindLongestXSort(std::vector<std::array<int,2>> &lengths, long long goal);   
    size_t FindLongestXHeap(std::vector<int> &lengths, long long goal); 

    void FilterCoverage(int min_coverage, int max_coverage, int max_diff_coverge);

    bool HasSupport(const Overlap &o, int count) const;
    bool HasAlignment(int a, int b, int end, int count, bool exceeding) const;
    std::unordered_set<const Overlap*> FindBestN(const std::pair<int, std::unordered_map<int, const Overlap*>> &groud) const;
  
    bool IsContained(const Overlap& o, std::array<int, 2> &rel);
    bool IsSupportContained(const Overlap& o, const std::array<int, 2> &rel);
    bool IsReserved(const Overlap &o) const {
        return GetOlReason(o).type == OlReason::RS_OK;
    }
    void UpdateFilteredRead(const std::unordered_map<Seq::Id, RdReason> &ignored);

    std::unordered_set<Seq::Id> ReservedReads();
    void ExtendReservedReads();

    /** Record internal state and variables */
    void Dump() const;      
    void DumpCoverage(const std::string &fname) const;
    void DumpFilteredReads(const std::string &fname) const;
    void DumpFilteredOverlaps(const std::string &fname) const;
    void DumpReadInfos(const std::string &fname, const std::unordered_map<int, ReadStatInfo> &readInfos) const;

    static bool BetterAlignedLength(const Overlap &o0, const Overlap &o1) { return o0.AlignedLength() > o1.AlignedLength(); }
    static void SetOlReason(const Overlap &o, OlReason rs);
    static OlReason GetOlReason(const Overlap &o);

public:
    static bool ParamToGenomeSize(const std::string& str, long long *v);
    static int Percentile(const std::vector<int>& data, double percent);
    static int FirstTrough(const std::vector<int>& data, size_t last, size_t k);
protected:
    double min_identity_raw_ { 70 };
    int max_overhang_raw_ { 1000 };
    double min_identity_{ -1 };         //!< 
    int min_length_{ 2500 };            //!< 
    int max_length_{ INT_MAX };         //!< 
    int min_aligned_length_{ 2500 };    //!< 
    int max_overhang_{ -1 };               
    int min_coverage_{ -1 };             //!< 
    int max_coverage_{ -1 };           //!< 
    int max_diff_coverage_{ -1 };      //!< 
    double coverage_discard_ { 0.01 };
    int bestn_{ 10 };                   //!< 
    std::string overlap_file_type_{ "" };
    int thread_size_{ 4 };           //!< 
    std::string output_directory_ {"."};
    std::string coverage_fname_ { "coverage.txt.gz" };  //!< variable this->coverages_

    double identity_global_deviation1_ { 98.0 };
    double identity_global_deviation2_ { 6 };
    double overhang_global_deviation1_ { 30 };
    double overhang_global_deviation2_ { 6 };
    double identity_local_deviation1_ { 99.0 };
    double identity_local_deviation2_ { 6 };
    double overhang_local_deviation1_ { 10 };
    double overhang_local_deviation2_ { 6 };
    int identity_local_condition_ {0};

    int local_low_coverage_ { 25 };

    int overhang_limit {0};

    std::string filtered_read_fname { "filtered_reads.txt.gz"};
    std::string filtered_overlap_fname {"filtered_overlaps.txt.gz"};
    std::string ifname_;                       
    std::string ofname_;                        
    long long genome_size_ {0};
    int coverage_ {40};
    std::array<int, 3> coverage_params_;

    OverlapStore ol_store_;
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> groups_;

    std::unordered_map<Seq::Id, std::array<int, 2>> coverages_;                         //!< record min and max base coverages of the reads
    std::unordered_map<Seq::Id, RdReason> filtered_reads_;                              //!< record filtered reads and reason for filtering
    std::unordered_map<Seq::Id, ReadStatInfo> read_infos_;
};

} // namespace fsa {
    
#endif // FSA_OVERLAP_FILTER_HPP

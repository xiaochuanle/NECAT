#ifndef FSA_CONTIG_LINK_HPP
#define FSA_CONTIG_LINK_HPP

#include <array>
#include <cassert>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "overlap_store.hpp"
#include "sequence.hpp"
#include "read_store.hpp"
#include "overlap.hpp"

namespace fsa {

class ContigLink {
public:
    /** Location of one contig in the other */
    enum class Loc{
        Left, Right, Containing, Contained, Equal, Abnormal
    };
public:
    struct Link {
        Link(const Overlap::Read &s, const Overlap::Read &t) : source(s), target(t) {}

        std::array<int, 2> Strand() { return std::array<int, 2>{strand_s2t, target.strand}; }

        bool SameStrand() const { return strand_s2t == target.strand; }
        static Loc Reverse(Loc t, bool direct);
        Loc Location(Seq::Id id, int err);
        Loc Location(const Overlap::Read &r, int err) { return Location(r.id, err);}
        Loc Location(int err) ;
        
        void CalcRealExpectOverlap();

        virtual size_t TotalLength() const {
            if (strand_s2t == target.strand) {
                return (size_t)(std::max(target.len, pos_s2t[3]) - std::min(0, pos_s2t[0]));
            }
            else {
                return (size_t)(std::max(target.len, pos_s2t[0]) - std::min(0, pos_s2t[3]));
            }
        }

        size_t TargetLength() const { return (size_t) target.len; }
        size_t SourceLength() const { return (size_t) source.len; }

        virtual std::vector<Seq::Area> GetSeqArea(Seq::Id id, char end) const = 0  ;

        std::array<int, 2> ol_expect;   // position at target 
        std::array<int, 2> ol_real;     // position at target 
        std::array<int, 4> pos_s2t;     // 
        int strand_s2t {-1};                 // 

        Overlap::Read source;
        Overlap::Read target;
        std::vector<const Overlap*> ols;
    };

    struct C2r2cLink : public Link{
        C2r2cLink(const Overlap& s, const Overlap &t);
        virtual std::vector<Seq::Area> GetSeqArea(Seq::Id id, char end) const;
        static bool Valid(const Overlap &s, const Overlap &t, int err);
        int AlignmentScore() const;
        int GapLength() const;
        Seq::Id ReadId() const { return ols[0]->a_.id; }
    };

    struct C2cLink : public Link {
        C2cLink(const Overlap& o, const Overlap::Read& s, const Overlap::Read& t);
        virtual std::vector<Seq::Area> GetSeqArea(Seq::Id id, char end) const;
        static bool Valid(const Overlap& o, int err);
    };

    size_t TotalLength() const {

        if (best_group != nullptr) {
            int sum = 0;
            for (auto l : (*best_group)) {
                sum += l->TotalLength();
            }
            return sum / (*best_group).size();
        } else {
            return 0;
        }
    }
    size_t TargetLength() const { return best_group != nullptr ?  (*best_group)[0]->TargetLength() : 0; }
    size_t SourceLength() const { return best_group != nullptr ?  (*best_group)[0]->SourceLength() : 0; }

    double Score() const;
    std::unordered_set<Seq::Id> GetRawreads() const;

    static bool SimpleValid(const Overlap& o, int end_err);
    static bool SimpleValid(const Overlap& o0, const Overlap &o1, int err);

    void Add(const Overlap& o, const Overlap::Read &source, const Overlap::Read& target, int err);
    void Add(const Overlap &source, Overlap &target);



    void AnalyzeLinks(int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size);
    void AnalyzeC2cLinks();
    void AnalyzeC2r2cLinks(int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size);
    

    std::vector<std::array<int,2>> GroupC2r2cLinks(std::vector<C2r2cLink*> &links, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size);
    std::vector<std::array<int,2>> GroupC2r2cLinks2(std::vector<C2r2cLink*> &links, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size);
    bool CheckGroupDeviation(std::vector<C2r2cLink*> &links, std::array<int, 2> &group, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length);
    bool GetBestGroupInWindow(std::vector<C2r2cLink*> &links, std::array<int, 2> &group, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size);
    bool CheckGroupCoverage(std::vector<C2r2cLink*>& links,  const std::array<int,2>& g, int read2ctg_max_overhang, int read2ctg_min_coverage);

    Link* Best() const { return best_c2c != nullptr ? (Link*)best_c2c : BestC2r2c(); }
    C2r2cLink* BestC2r2c() const {
        return best_group != nullptr ? (*best_group)[0] : nullptr;
    }

    //bool Valid() const { return best_c2c != nullptr || best_group != nullptr; }
    bool Valid() const { return best_group != nullptr; }
    Seq::Id Target() const;
    Seq::Id Source() const;
    
    void Removed(bool torf) { removed_ = torf; }
    bool Removed() const { return removed_; }
    
    std::vector<C2cLink> c2c_links;
    std::vector<C2r2cLink> c2r2c_links;

    std::vector<std::vector<C2r2cLink*>> groups;
    C2cLink *best_c2c {nullptr};
    std::vector<C2r2cLink*>* best_group {nullptr};

    bool removed_ { false };

    
};

} // namespace fsa {

#endif // FSA_CONTIG_GRAPH_HPP  

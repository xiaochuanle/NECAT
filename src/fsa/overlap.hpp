#ifndef FSA_OVERLAP_HPP
#define FSA_OVERLAP_HPP

#include <string>
#include <array>

#include "sequence.hpp"

namespace fsa {

class Overlap {
public:
	struct Read {
		int id;
		int strand;
		int start;
		int end;
		int len;
	};
public:


    enum class Loc {
        Left, Right, Contained, Equal, Containing, Abnormal
    };
public:
    Overlap() = default;

    std::string ToM4Line() const;
    bool FromM4Line(const std::string &line);
    
    bool CheckEnd(int error) const;

    bool SameDirect() const { return a_.strand == b_.strand; }
    size_t AlignedLength() const { return (a_.end - a_.start + b_.end - b_.start) / 2;  }

    bool IsContaining(int err) const;
    bool IsContained(int err) const;
    bool IsAbnormal(int err) const;
    bool IsProper(int err) const;
    
    Loc Location(int err) const;
    Loc Location(Seq::Id id, int err) const ;
    Loc Location(const Read& r, int err) const { return Location(r.id, err);}

    static std::array<int, 2> Overhang(const Overlap &o, Loc loc);

    std::array<int, 2> Overhang() const;
    std::array<int, 2> Overhang2() const;

    template<int N>
    std::array<int, N>  MappingToSource(const std::array<int, N>& pos) const {
        return Mapping<N>(a_, b_, pos);
    }

    template<int N>
    std::array<int, N>  MappingToTarget(const std::array<int, N>& pos) const {
        return Mapping<N>(b_, a_, pos);
    }

    template<int N>
    static std::array<int, N> Mapping(const Read& a, const Read& b, const std::array<int, N> &bpos);

    static int ReverseStrand(int s) { return s == 0 ? 1 : 0;  }
    static Loc ReverseLocation(Loc loc, bool direct=true);

	Read a_;
	Read b_;
	int score_;
	double identity_;
    mutable long long attached {0};   // used by other, and avoid building map (Overlap -> info)
};


template<int N>
std::array<int, N> Overlap::Mapping(const Read& a, const Read& b, const std::array<int, N> &bpos) {
    std::array<int, N> apos;
    if (a.strand == b.strand) {
        for (size_t i=0; i<N; ++i) {
            if (bpos[i] < b.start) {
                apos[i] = a.start + (bpos[i] - b.start);
            } else if (bpos[i] >= b.start && bpos[i] < b.end) {
                apos[i] = a.start + (long long)(bpos[i] - b.start) * (a.end - a.start) / (b.end - b.start);
            } else {
                assert(bpos[i] >= b.end);
                apos[i] = a.end + (bpos[i] - b.end);
            }
        }
    }
    else {
        for (size_t i=0; i<N; ++i) {
            if (bpos[i] > b.end) {
                apos[i] = a.start - (bpos[i] - b.end);
            } else if (bpos[i] <= b.end && bpos[i] > b.start) {
                apos[i] = a.start - (long long)(bpos[i]-b.end) * (a.end - a.start) / (b.end - b.start);
            } else { 
                assert(bpos[i] <= b.start);
                apos[i] = a.end - (bpos[i] - b.start);
            }
        }
    }
    return apos;
}

} // namespace fsa {

#endif  // FSA_OVERLAP_HPP

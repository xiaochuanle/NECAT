#ifndef FSA_ALIGN_SIMPLE_ALIGN_HPP
#define FSA_ALIGN_SIMPLE_ALIGN_HPP

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace fsa {

class SimpleAlign {
public:
    struct KmerList {
        KmerList(const std::string& target, uint8_t k);

        struct Node {
            Node(size_t h=-1, size_t t=-1, size_t c=0) : head(h), tail(t), count(c) {}
            size_t head;
            size_t tail;
            size_t count;
        };

        size_t FirstPosition(uint32_t bkmer) const {
            auto it = list.find(bkmer);
            return it != list.end() ? it->second.head : next_same.size();
        }
        size_t NextPosition(size_t pos) const { return next_same[pos]; }
        bool ValidPosition(size_t pos) const { return pos < next_same.size(); }

        size_t Size(uint32_t bkmer) const {
            auto it = list.find(bkmer);
            return it != list.end() ? it->second.count : 0;
        }

        

        std::vector<size_t> next_same;

        std::unordered_map<uint32_t, Node> list;
    };

    struct KmerMatch {
        void Add(size_t qpos, size_t tpos) { matchs.push_back({qpos, tpos}); }
        size_t Size() const { return matchs.size(); }
        const std::array<size_t,2>& Get(size_t i) const { return matchs[i]; }
        int Diff(size_t i) const { return Diff(matchs[i]);}
        static int Diff(const std::array<size_t, 2> &m) { return int(m[0]) - int(m[1]);}
        std::array<int, 2> DiffMinMax() const;

        std::vector<std::array<size_t, 2>> matchs;  // (query target)
    };

    struct Range {
        //Range(size_t qs, size_t qe, size_t ts, size_t te) : query({qs, qe}), target({ts, te}) {}
        size_t QueryLength() const { return query[1] - query[0];}
        size_t QueryStart() const { return query[0];}
        size_t QueryEnd() const { return query[1];}
        size_t TargetLength() const { return target[1] - target[0];}
        size_t TargetStart() const { return target[0];}
        size_t TargetEnd() const { return target[1];}
        std::array<size_t, 2> query ;
        std::array<size_t, 2> target;
    };

    struct Result {
        int distance {0};
        int query_start {0};
        int query_end {0};
        int target_start {0};
        int target_end {0};
        std::string aligned_query;
        std::string aligned_target;
    };

public:
    SimpleAlign(const std::string &target, uint8_t k);

    Result Align(const std::string &query, int band, bool detail=false);
    Result Align(const std::string &query, const Range &range, int band, bool detail=false);
    Result Align(const std::string &query, const std::string &target, int band, bool detail=false);
    
public:    
    static uint32_t KmerMask(const uint8_t k);
    static uint32_t CalcKmer(const char* seq, size_t k);
    static uint32_t MoveKmer(uint32_t bkmer, size_t k, char e);

    KmerMatch FindKmerMatch(const std::string &query, size_t stride);
    Range FindCandidateRange(const KmerMatch &match, size_t binsize, size_t threshold);

    const std::string &target_;
    uint8_t k_;
    KmerList target_bkmers_;
    uint8_t query_stride_;
};

} // namespace fsa {

#endif // FSA_ALIGN_SIMPLE_ALIGN_HPP

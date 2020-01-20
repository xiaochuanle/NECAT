#include "simple_align.hpp"

#include <algorithm>
#include <cassert>
#include <list>

#include "utility.hpp"

namespace fsa {

inline uint8_t Acgt2Num(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: assert(!"never come here"); return 0xFF;
    }
}

SimpleAlign::KmerList::KmerList(const std::string &target, uint8_t k) {

    assert(k <= 16 && target.size() >= k);

    next_same.assign(target.size() - k + 1, target.size() - k + 1);

    assert(next_same.size() == target.size() - k + 1);

    uint32_t bkmer = CalcKmer(target.c_str(), k);
    list[bkmer] = Node(0, 0, 1);

    for ( size_t i=1; i< target.size()-k+1; ++i) {
        bkmer = MoveKmer(bkmer, k, target[i+k-1]);
        auto it = list.find(bkmer);
        if (it == list.end()) {
            list[bkmer] = Node(i, i, 1);
        } else {
            it->second.count++;
            next_same[it->second.tail] = i;
            it->second.tail = i;
        }
    }
}



std::array<int, 2> SimpleAlign::KmerMatch::DiffMinMax() const {
    assert(Size() > 0);
    auto minmax = std::minmax_element(matchs.begin(), matchs.end(), [](const std::array<size_t, 2> &a, const std::array<size_t, 2> &b) {
        return Diff(a) < Diff(b); 
    });

    return {Diff(*minmax.first), Diff(*minmax.second)};
}

SimpleAlign::SimpleAlign(const std::string &target, uint8_t k) 
    : target_(target), k_(k), target_bkmers_(target, k) {

    assert(k <= 16 && target_.size() >= k);
    query_stride_ = k_ / 2;
    
}

SimpleAlign::KmerMatch SimpleAlign::FindKmerMatch(const std::string &query, size_t stride) {
    KmerMatch match;
    for (size_t i=0; i<query.size() - k_ + 1; i += stride) {
        auto bkmer = CalcKmer(query.c_str()+i, k_);

        for (size_t pos = target_bkmers_.FirstPosition(bkmer); 
            target_bkmers_.ValidPosition(pos); 
            pos = target_bkmers_.NextPosition(pos)) {

            match.Add(i, pos);
        }
    }
    return match;
}

SimpleAlign::Result SimpleAlign::Align(const std::string &query, int band_tol, bool detail) {
    KmerMatch match = FindKmerMatch(query, query_stride_);
    auto range = FindCandidateRange(match, k_*5, 12);
    return Align(query, range, band_tol, detail);

}

SimpleAlign::Result SimpleAlign::Align(const std::string &query, const Range &range, int band_tol, bool detail) {
    Result r = Align(query.substr(range.QueryStart(), range.QueryLength()), 
                 target_.substr(range.TargetStart(), range.TargetLength()), band_tol, detail);

    r.query_start += range.QueryStart();
    r.query_end += range.QueryStart();
    r.target_start += range.TargetStart();
    r.target_end += range.TargetStart();

    return r;
}

SimpleAlign::Result SimpleAlign::Align(const std::string &query, const std::string &target_, int band_tolerance, bool detail) {

    Result result = Result();

    struct PathNode {
        int pdxy;   //    prev diff bwtreen x and y
        int x;
        int y;
        int len { 0 };
    };

    const int MAX_DIST = (int) (0.3*(query.size() + target_.size()));
    const int LOCAL_BAND_WIDTH = band_tolerance * 2;

    std::vector<int> xAtBandBuff(MAX_DIST*2+1, 0);
    int* xAtBand = &xAtBandBuff[0] + MAX_DIST;  // support  negative index [-MAX_DIST, MAX_DIST] 


    // We should probably use hashmap to store the backtracing information to save memory allocation time
    // This O(MN) block allocation scheme is convient for now but it is slower for very long sequences
    std::unordered_map<std::array<int,2>, PathNode, ArrayHash<int,2>, ArrayCompare<int,2>> path_nodes;
    std::array<int,2> last_node {0, -1}; // (diff of x & y, distance), distance == -1 mean don't find aligned strings 

    int best_xy = 0;
    std::array<int,2> dxy_range {0, 0};
    
    //printf("new len %zd, %zd\n", query.size(), target_.size());
    for (int d=0; d < MAX_DIST; d ++ ) {
        for (int dxy = dxy_range[0]; dxy <= dxy_range[1];  dxy += 2) {
            PathNode node;

            if ( dxy == dxy_range[0] || (dxy != dxy_range[1]  && xAtBand[dxy+1] > xAtBand[ dxy-1])) {
                node.pdxy = dxy + 1;
                node.x = xAtBand[ dxy + 1];
            }  else {
                assert(dxy == dxy_range[1] || (dxy != dxy_range[0] && xAtBand[dxy+1] <= xAtBand[ dxy-1]));
                node.pdxy = dxy - 1;
                node.x = xAtBand[ dxy - 1 ] + 1;
            }

            node.y = node.x - dxy;

            while ( size_t(node.x + node.len) < query.size() && 
                    size_t(node.y + node.len) < target_.size() && 
                    query[node.x+node.len] == target_[node.y+node.len]) {
                node.len ++;
            }

            assert(path_nodes.find({dxy,d}) ==path_nodes.end());          
            path_nodes[{dxy,d}] = node;

            xAtBand[dxy] = node.x + node.len;
            if ( node.x + node.len + node.y + node.len> best_xy) {
                best_xy = node.x + node.len + node.y + node.len;
            }

            if ( node.x + node.len  >= query.size() || node.y + node.len  >= target_.size()) {
                last_node[0] = dxy;
                last_node[1] = d;
                break;
            }
        }

        if (!(last_node[1] >= 0)) { 
            // 
            for (int dxy = dxy_range[0]; dxy <= dxy_range[1]; dxy+=2) {
                int xy = xAtBand[dxy] + xAtBand[dxy] - dxy; // x+y
                if (xy >= best_xy - band_tolerance) {
                    dxy_range[0] = dxy-1;
                    break;
                }
            }

            for (int dxy = dxy_range[1]; dxy>= dxy_range[0]; dxy-=2) {
                int xy = xAtBand[dxy] + xAtBand[dxy] - dxy; // x+y
                if (xy >= best_xy - band_tolerance) {
                    dxy_range[1] = dxy+1;
                    break;
                }
            }

            
            if (dxy_range[1]-dxy_range[0] > LOCAL_BAND_WIDTH ) {
                break;
            }
            //printf("new k range %d, %d\n", dxy_range[0], dxy_range[1]);
            
        } else {
            break;
        }
    }

    if (last_node[1] >= 0) {
        const auto& node = path_nodes[last_node];
        result.query_end = node.x + node.len;
        result.target_end = node.y + node.len;
        result.distance = last_node[1];
        result.query_start = 0;
        result.target_start = 0;

        if (detail > 0) {
            // compute aligned strings, like
            //  query:  ATCG-GCAT
            // target:  A-CGTGC-T

            std::list<std::array<int,2>> path; // list of (x,y)
            std::array<int,2> curr_node = last_node;
            while (curr_node[1] >= 0) {

                assert(path_nodes.find(curr_node) !=path_nodes.end());    
                auto &node = path_nodes[curr_node];         
   
                path.push_front({node.x+node.len, node.y + node.len});
                path.push_front({node.x, node.y});
                curr_node = {node.pdxy, curr_node[1]-1};
            }
            
            auto cxy = path.front();
            assert(result.query_start == cxy[0] && result.target_start == cxy[1]);

            for (const auto& nxy : path) {
                assert(nxy[0] >= cxy[0] && nxy[1] >= cxy[1]);
                    //printf("new nx,ny %d, %d\n", nxy[0], nxy[1]);
 
                if (nxy[0] == cxy[0] && nxy[1] != cxy[1]){ //advance in y
                    result.aligned_query.insert(result.aligned_query.size(), nxy[1]-cxy[1], '-');
                    result.aligned_target.insert(result.aligned_target.size(), target_, cxy[1], nxy[1]-cxy[1]);
                } else if (nxy[0] != cxy[0] && nxy[1] == cxy[1]){ //advance in x
                    result.aligned_query.insert(result.aligned_query.size(), query, cxy[0], nxy[0]-cxy[0]);
                    result.aligned_target.insert(result.aligned_target.size(), nxy[0]-cxy[0], '-');
                } else if (nxy[0] != cxy[0] && nxy[1] != cxy[1]) {
                    result.aligned_query.insert(result.aligned_query.size(), query, cxy[0], nxy[0]-cxy[0]);
                    result.aligned_target.insert(result.aligned_target.size(), target_, cxy[1], nxy[1]-cxy[1]);
                } else {
                    assert(cxy[0] == nxy[0] && cxy[1] == nxy[1]);
                    // No need to deal with it
                }
                cxy = nxy;
            }
            assert(result.query_end == cxy[0] && result.target_end == cxy[1]);
        }
        
    }


    return result;


}

SimpleAlign::Range SimpleAlign::FindCandidateRange(const KmerMatch &match, size_t binsize, size_t threshold) {

    Range range{{0,0},{0,0}};

    auto diff_minmax = match.DiffMinMax();
    //printf("new minmax %d, %d, %zd\n", diff_minmax[0], diff_minmax[1], binsize);
    std::vector<size_t> coor;
    std::vector<size_t> bins((diff_minmax[1] - diff_minmax[0]) / binsize + 1);

    for (size_t i=0; i< match.Size(); i++) {
        bins[(match.Diff(i) - diff_minmax[0] ) / binsize] += 1;
    }


    auto bin_minmax = std::minmax_element(bins.begin(), bins.end());
    size_t bin_max_index = bin_minmax.second - bins.begin();

    if (bins[bin_max_index] > threshold) {
        for (size_t i=0; i < match.Size(); ++i) {
            auto bin_index = (match.Diff(i) - diff_minmax[0]) / binsize;
            if (bin_index + 5 >= bin_max_index && bin_index <= bin_max_index + 5 && bins[bin_index] >= threshold) {
                coor.push_back(i);
            }
        }

    }

    //printf("new ss %zd\n", coor.size());;
    if (coor.size() > 0) {
        range = {{match.Get(coor[0])[0], match.Get(coor[0])[0]}, {match.Get(coor[0])[1], match.Get(coor[0])[1]} };

        int max_score = 0;
        int score = 0;
        size_t start = 0;
        for (size_t i=1; i<coor.size(); ++i) {
            score += 32 - (match.Get(coor[i])[0] - match.Get(coor[i-1])[0]);
            if (score < 0) {
                score = 0;
                start = i;
            } else if (score > max_score) {
                range = {{match.Get(coor[start])[0], match.Get(coor[i])[0]}, {match.Get(coor[start])[1], match.Get(coor[i])[1]} };
                max_score = score;
            }
        }
    }

    return range;
}


uint32_t SimpleAlign::KmerMask(uint8_t k) { 
    assert(k <= 16);
    static uint32_t s_mark[] = {
        0x00, 0x03, 0x0f, 0x3f, 0xff, 
        0x03ff, 0x0fff, 0x3fff, 0xffff,
        0x03ffff, 0x0fffff, 0x3fffff, 0xffffff,
        0x03ffffff, 0x0fffffff, 0x3fffffff, 0xffffffff};
    return s_mark[k];
    //return (1 << (2*k)) - 1;  
}

uint32_t SimpleAlign::CalcKmer(const char* seq, size_t k) {
    assert(k <= 16);
    uint32_t kmer = 0;
    for (size_t i=0; i<k; ++i) {
        kmer <<= 2;
        kmer |= Acgt2Num(seq[i]) & 0x03;
    }
    return kmer;
}




uint32_t SimpleAlign::MoveKmer(uint32_t bkmer, size_t k, char e) {
    return ((bkmer << 2) | Acgt2Num(e)) & KmerMask(k);
}

} // namespace fsa {

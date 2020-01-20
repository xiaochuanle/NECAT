#include "contig_link.hpp"

#include <array>
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <iostream>
#include <tuple>
#include <cmath>

#include "overlap.hpp"
#include "utility.hpp"

namespace fsa {

ContigLink::Loc ContigLink::Link::Reverse(Loc t, bool direct) {
    switch (t) {
        case Loc::Left: return direct ? Loc::Right : Loc::Left;
        case Loc::Right: return direct ? Loc::Left : Loc::Right;
        case Loc::Equal: return Loc::Equal;
        case Loc::Containing: return Loc::Contained;
        case Loc::Contained: return Loc::Containing;
        case Loc::Abnormal: return Loc::Abnormal;
        default:   { assert(!"never come here");    return Loc::Abnormal; }
    }
}

ContigLink::Loc ContigLink::Link::Location(Seq::Id id, int err) {
    assert(id == source.id || id == target.id);
    Loc loc = Location(err);
    return id == source.id ? loc : Reverse(loc, SameStrand());
}

ContigLink::Loc ContigLink::Link::Location(int err) {
    auto compare = [err](int a, int b) { 
        if (a < b - err) return -1;
        else if (a > b + err) return 1;
        else return 0;
    };


    if (strand_s2t == target.strand) {
        auto s_cmp_t_start = compare(pos_s2t[0], 0);
        auto s_cmp_t_end = compare(pos_s2t[3], target.len);

        if (s_cmp_t_start < 0 && s_cmp_t_end < 0) {
            return Loc::Left;
        }

        if (s_cmp_t_start > 0 && s_cmp_t_end > 0) {
            return Loc::Right;
        }

        if (s_cmp_t_start == 0 && s_cmp_t_end == 0) {
            return Loc::Equal;
        }

        if (s_cmp_t_start <=0 && s_cmp_t_end >= 0 && (s_cmp_t_start != 0 || s_cmp_t_end != 0)) {
            return Loc::Containing;
        }

        if (s_cmp_t_start >=0 && s_cmp_t_end <= 0 && (s_cmp_t_start != 0 || s_cmp_t_end != 0)) {
            return Loc::Contained;
        }

    } else {
        
        auto s_cmp_t_start = compare(pos_s2t[3], 0);
        auto s_cmp_t_end = compare(pos_s2t[0], target.len);


        if (s_cmp_t_start < 0 && s_cmp_t_end < 0) {
            return Loc::Left;
        }

        if (s_cmp_t_start > 0 && s_cmp_t_end > 0) {
            return Loc::Right;
        }

        if (s_cmp_t_start == 0 && s_cmp_t_end == 0) {
            return Loc::Equal;
        }

        if (s_cmp_t_start <=0 && s_cmp_t_end >= 0 && (s_cmp_t_start != 0 || s_cmp_t_end != 0)) {
            return Loc::Containing;
        }

        if (s_cmp_t_start >=0 && s_cmp_t_end <= 0 && (s_cmp_t_start != 0 || s_cmp_t_end != 0)) {
            return Loc::Contained;
        }

    }

    return Loc::Abnormal;
}

void ContigLink::Link::CalcRealExpectOverlap() {
    assert(strand_s2t == 0 || strand_s2t == 1);

    if (SameStrand())  {
        ol_real[0] = pos_s2t[1] >= target.start ? pos_s2t[1] : target.start;
        ol_real[1] = pos_s2t[2] <= target.end ? pos_s2t[2] : target.end;

        ol_expect[0] = pos_s2t[0] >= 0 ? pos_s2t[0] : 0;
        ol_expect[1] = pos_s2t[3] <= target.len ? pos_s2t[3] : target.len;
    
    } else {
        ol_real[0] = pos_s2t[2] >= target.start ? pos_s2t[2] : target.start;
        ol_real[1] = pos_s2t[1] <= target.end ? pos_s2t[1] : target.end;

        ol_expect[0] = pos_s2t[3] >= 0 ? pos_s2t[3] : 0;
        ol_expect[1] = pos_s2t[0] <= target.len ? pos_s2t[0] : target.len;
    }
}

ContigLink::C2cLink::C2cLink(const Overlap &o, const Overlap::Read &s, const Overlap::Read &t) 
    : Link(s, t) {

    ols.push_back(&o);

    pos_s2t = o.a_.id == s.id ? o.MappingToTarget<4>(std::array<int, 4>{0, s.start, s.end, s.len})
        : o.MappingToSource<4>(std::array<int, 4>{0, s.start, s.end, s.len});
    strand_s2t = s.strand;

    CalcRealExpectOverlap();
}

ContigLink::C2r2cLink::C2r2cLink(const Overlap &s, const Overlap &t) 
    : Link(s.b_, t.b_) {

    ols.push_back(&s);
    ols.push_back(&t);

    pos_s2t = t.MappingToTarget<4>(s.MappingToSource<4>(std::array<int, 4>{0, s.b_.start, s.b_.end, s.b_.len}));
    strand_s2t = s.a_.strand == t.a_.strand ? s.b_.strand : Overlap::ReverseStrand(s.b_.strand);

    CalcRealExpectOverlap();
}

std::vector<Seq::Area> ContigLink::C2r2cLink::GetSeqArea(Seq::Id id, char end) const {
    if (target.id == id) {

        auto pos_s2t = ols[1]->MappingToTarget<4>(ols[0]->MappingToSource<4>(std::array<int, 4>{0, source.start, source.end, source.len}));
        auto strand_s2t = ols[0]->a_.strand == ols[1]->a_.strand ? source.strand : Overlap::ReverseStrand(source.strand);

        if (end == 'B') {
            int s = strand_s2t == target.strand ? pos_s2t[0] : pos_s2t[3];
            //int e = strand_s2t == target.strand ? pos_s2t[3] : pos_s2t[0];
            assert(pos_s2t[0] >= 0);

            if (s >= target.len) {
                auto t2r = ols[1]->MappingToSource<2>(std::array<int, 2>{s, target.len });
                if (t2r[0] < t2r[1]) {
                    return { Seq::Area{ ols[1]->a_.id, 0, t2r[0] , t2r[1] },  Seq::Area{ id, 1, 0, target.len } };
                }
                else {
                    return { Seq::Area{ ols[1]->a_.id, 1, t2r[1], t2r[0] },  Seq::Area{ id, 1, 0, target.len } };
                }
            }
            else {

                return { Seq::Area{ id, 1, 0, s } };
            }
        }
        else {
            assert(end == 'E');
            //int s = strand_s2t == target.strand ? pos_s2t[0] : pos_s2t[3];
            int e = strand_s2t == target.strand ? pos_s2t[3] : pos_s2t[0];
            assert(e <= target.len);

            if (e < 0) {
                auto t2r = ols[1]->MappingToSource<2>(std::array<int, 2>{e, 0 });
                if (t2r[0] < t2r[1]) {
                    return { Seq::Area{ ols[1]->a_.id, 0, t2r[0] , t2r[1] },  Seq::Area{ id, 0, 0, target.len } };
                }
                else {
                    return { Seq::Area{ ols[1]->a_.id, 1, t2r[1], t2r[0] },  Seq::Area{ id, 0, 0, target.len } };
                }
            }
            else {

                return { Seq::Area{ id, 0, e, target.len } };
            }
        }
    }
    else {
        assert(source.id == id);
        auto pos_t2s = ols[0]->MappingToTarget<4>(ols[1]->MappingToSource<4>(std::array<int, 4>{0, target.start, target.end, target.len}));
        auto strand_t2s = ols[0]->a_.strand == ols[1]->a_.strand ? target.strand : Overlap::ReverseStrand(target.strand);

        if (end == 'B') {
            int s = strand_t2s == source.strand ? pos_t2s[0] : pos_t2s[3];
            //int e = strand_t2s == source.strand ? pos_t2s[3] : pos_t2s[0];

            assert(s >= 0);
            if (s >= source.len) {
                auto t2r = ols[0]->MappingToSource<2>(std::array<int, 2>{s, source.len });
                if (t2r[0] < t2r[1]) {
                    return { Seq::Area{ ols[0]->a_.id, 0, t2r[0] , t2r[1] },  Seq::Area{ id, 1, 0, source.len } };
                }
                else {
                    return { Seq::Area{ ols[0]->a_.id, 1, t2r[1], t2r[0] },  Seq::Area{ id, 1, 0, source.len } };
                }
            }
            else {

                return { Seq::Area{ id, 1, 0, s } };
            }
        }
        else {
            assert(end == 'E');
            //int s = strand_t2s == source.strand ? pos_t2s[0] : pos_t2s[3];
            int e = strand_t2s == source.strand ? pos_t2s[3] : pos_t2s[0];

            assert(e <= source.len);

            if (e < 0) {
                auto t2r = ols[0]->MappingToSource<2>(std::array<int, 2>{e, 0});
                if (t2r[0] < t2r[1]) {
                    return { Seq::Area{ ols[0]->a_.id, 0, t2r[0] , t2r[1] },  Seq::Area{ id, 0, 0, source.len } };
                }
                else {
                    return { Seq::Area{ ols[0]->a_.id, 1, t2r[1], t2r[0] },  Seq::Area{ id, 0, 0, source.len } };
                }
            }
            else {

                return { Seq::Area{ id, 0, e, source.len } };
            }
        }
    }
}


int ContigLink::C2r2cLink::AlignmentScore() const {
    return (ols[0]->AlignedLength()*ols[0]->identity_ + ols[1]->AlignedLength()*ols[1]->identity_)/100;
}

int ContigLink::C2r2cLink::GapLength() const {
    if (SameStrand()) {
        return pos_s2t[0] < 0 ? 0 - pos_s2t[3] : pos_s2t[0] - target.len; 
    } else {
        return pos_s2t[3] < 0 ? 0 - pos_s2t[0] : pos_s2t[3] - target.len; 
    }
}

std::vector<Seq::Area> ContigLink::C2cLink::GetSeqArea(Seq::Id id, char end) const {
    if (target.id == id) {
        if (end == 'B') {
            return { Seq::Area{id, 1, 0, target.start} };
        }
        else {
            assert(end == 'E');

            return { Seq::Area{ id, 0, target.end, target.len } };
        }
    }
    else {
        if (end == 'B') {
            return { Seq::Area{ id, 1, 0, source.start } };
        }
        else {
            assert(end == 'E');

            return { Seq::Area{ id, 0, source.end, source.len } };
        }
    }
}

double ContigLink::Score() const {
    return best_group != nullptr ? 
        std::accumulate(best_group->begin(), best_group->end(), 0.0, [](double a, const C2r2cLink* b) { return a + b->AlignmentScore();}) :
        0.0;
}
std::unordered_set<Seq::Id> ContigLink::GetRawreads() const {
    std::unordered_set<Seq::Id> rr;

    if (best_group != nullptr) {
        for (auto i : *best_group) {
            rr.insert(i->ReadId());
        }
    }
    return rr;
}

bool ContigLink::SimpleValid(const Overlap& o, int err) {
    Overlap::Loc loc = o.Location(err);
    return loc == Overlap::Loc::Left || loc ==  Overlap::Loc::Right ||
        loc == Overlap::Loc::Contained || loc ==  Overlap::Loc::Containing || loc ==  Overlap::Loc::Equal;
}


bool ContigLink::SimpleValid(const Overlap &s, const Overlap &t, int err) {
    Overlap::Loc loc_s = s.Location(s.b_, err);
    Overlap::Loc loc_t = t.Location(t.b_, err);

    return loc_s != Overlap::Loc::Abnormal && loc_t != Overlap::Loc::Abnormal && (
           (loc_s == Overlap::Loc::Left && loc_t == Overlap::Loc::Right) || 
           (loc_s == Overlap::Loc::Right && loc_t == Overlap::Loc::Left) || 
           (loc_s == Overlap::Loc::Equal || loc_s == Overlap::Loc::Contained) ||
           (loc_t == Overlap::Loc::Equal || loc_t == Overlap::Loc::Contained));
}


void ContigLink::Add(const Overlap& o, const Overlap::Read &source, const Overlap::Read& target, int err) {
    C2cLink link(o, source, target);

    for (const auto &i : c2c_links) {
        if (i.SameStrand() == link.SameStrand() && i.pos_s2t[0] - link.pos_s2t[0] < err) {
            // remove duplicated overlaps
            return;
        }
    }
    c2c_links.push_back(link);
}


void ContigLink::Add(const Overlap &source, Overlap &target) {
    // duplicated overlaps were removed by grouping operation
    C2r2cLink link(source, target);

    c2r2c_links.push_back(link);
}



void ContigLink::AnalyzeLinks(int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size) {
    AnalyzeC2cLinks();
    AnalyzeC2r2cLinks(read2ctg_max_overhang, ctg2ctg_max_overhang, read2ctg_min_coverage, read2ctg_min_aligned_length, window_size);
  
  
    if (c2c_links.size() > 0) {
        for (auto &g : groups) {
            for (auto &c : c2c_links) {
                if (g[0]->SameStrand() == c.SameStrand() && g[0]->pos_s2t[0] - c.pos_s2t[0] < read2ctg_max_overhang) {
                    best_group = &g;
                    best_c2c = &c;
                    break;
                }
            }
            if (best_group != nullptr) {
                break;
            }
        }
    } else {
        if (groups.size() > 0) {
            best_group = &groups[0];
        }
    }

    if (best_group != nullptr) {
        std::sort(best_group->begin(), best_group->end(), [](const C2r2cLink* a, const C2r2cLink* b) {
            return a->AlignmentScore() > b->AlignmentScore();
        });
    }
}

void ContigLink::AnalyzeC2cLinks() {
    std::sort(c2c_links.begin(), c2c_links.end(), [](const C2cLink&a, const C2cLink&b) { 
        return a.ols[0]->AlignedLength() > b.ols[0]->AlignedLength();
    });
}

void ContigLink::AnalyzeC2r2cLinks(int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size) {

    std::array<std::vector<C2r2cLink*>, 2> offset;
    for (auto &c : c2r2c_links) {
        // TODO assert(!c.ols[0]->IsContained(read2ctg_max_overhang) && !c.ols[1]->IsContained(read2ctg_max_overhang));
        
        if (c.Location(read2ctg_max_overhang) != ContigLink::Loc::Abnormal) {
            if (c.SameStrand()) {
                offset[0].push_back(&c);
            } else {
                offset[1].push_back(&c);
            }
        }
    }

    auto group0 = GroupC2r2cLinks(offset[0], read2ctg_max_overhang, ctg2ctg_max_overhang, read2ctg_min_coverage, read2ctg_min_aligned_length, window_size);
    auto group1 = GroupC2r2cLinks(offset[1], read2ctg_max_overhang, ctg2ctg_max_overhang, read2ctg_min_coverage, read2ctg_min_aligned_length, window_size);

    for (auto r : group0) {
        std::vector<C2r2cLink*> g;
        for (int i=r[0]; i<r[1]; ++i) {
            g.push_back(offset[0][i]);
        }
        groups.push_back(g);
    }

    for (auto r : group1) {
        std::vector<C2r2cLink*> g;
        for (int i=r[0]; i<r[1]; ++i) {
            g.push_back(offset[1][i]);
        }
        groups.push_back(g);
    }

    std::sort(groups.begin(), groups.end(), [](const std::vector<C2r2cLink*> &a, const std::vector<C2r2cLink*> &b) {
        return std::accumulate(a.begin(), a.end(), 0.0, [](double a, const C2r2cLink* b) { return a + b->AlignmentScore();}) > 
               std::accumulate(b.begin(), b.end(), 0.0, [](double a, const C2r2cLink* b) { return a + b->AlignmentScore();});
    });
    
}


std::vector<std::array<int,2>>  ContigLink::GroupC2r2cLinks(std::vector<C2r2cLink*> &links, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size) {

    int split_threshold = 1000;//read2ctg_max_overhang * 2;
        
    std::sort(links.begin(), links.end(), [](const C2r2cLink* a, const C2r2cLink* b) { return a->pos_s2t[0] < b->pos_s2t[0]; });

    std::vector<int> split_pos(1, 0);
    for (size_t i = 1; i<links.size(); ++i) {
        if (links[i]->pos_s2t[0] - links[i-1]->pos_s2t[0] > split_threshold) {
            split_pos.push_back(i);
        }
    }
    split_pos.push_back(links.size());

    std::vector<std::array<int,2>> groups;
    for (size_t i=1; i<split_pos.size(); ++i) {
        std::array<int, 2> g = {split_pos[i-1], split_pos[i]};
        if (GetBestGroupInWindow(links, g, read2ctg_max_overhang, ctg2ctg_max_overhang, read2ctg_min_coverage, read2ctg_min_aligned_length, window_size) &&
            CheckGroupCoverage(links, g, read2ctg_max_overhang, read2ctg_min_coverage)) {
            groups.push_back(g);
        }
    }

    return groups;
}


std::vector<std::array<int,2>>  ContigLink::GroupC2r2cLinks2(std::vector<C2r2cLink*> &links, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size) {

    std::sort(links.begin(), links.end(), [](const C2r2cLink* a, const C2r2cLink* b) { return a->pos_s2t[0] < b->pos_s2t[0]; });

    std::array<int, 2> range{0, 0};
    double score = -1;
    double curr_score = 0.0;
    size_t curr_start = 0;
    for (size_t i=0; i < links.size(); ++i) {
        curr_score += links[i]->AlignmentScore();
        for (; curr_start <= i; ++curr_start) {
            if (links[i]->pos_s2t[0] - links[curr_start]->pos_s2t[0] > window_size) {
                curr_score -= links[curr_start]->AlignmentScore();
            } else {
                break;
            }
        }

        if (curr_score > score) {
            range[0] = curr_start;
            range[1] = i + 1;
            score = curr_score;
        }

    }

    std::vector<std::array<int,2>> groups;
    if (range[1] - range[0] >= read2ctg_min_coverage) {
        groups.push_back(range);
    }

    return groups;
}

bool ContigLink::GetBestGroupInWindow(std::vector<C2r2cLink*> &links, std::array<int, 2> &group, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length, int window_size) {
     auto get_aligned_length = [](const std::vector<C2r2cLink*>& offset, const std::array<int, 2>& range) {
        std::array<int, 2> aligned = {0, 0};    
        for (auto i=range[0]; i< range[1]; ++i) {
            if (aligned[0] < offset[i]->source.end -  offset[i]->source.start)
                aligned[0] = offset[i]->source.end -  offset[i]->source.start;
            if (aligned[1] < offset[i]->target.end -  offset[i]->target.start)
                aligned[1] = offset[i]->target.end -  offset[i]->target.start;
        }
        return aligned;
    };



    auto best_in_window = [&links, window_size](std::array<int, 2> &g) {

        std::array<int, 2> range{0, 0};
        double best_score = -1;
        double score = 0.0;
        int start = g[0];
        for (int i=g[0]; i < g[1]; ++i) {
            score += links[i]->AlignmentScore();
            for (; start <= i; ++start) {
                if (links[i]->pos_s2t[0] - links[start]->pos_s2t[0] > window_size) {
                    score -= links[start]->AlignmentScore();
                } else {
                    break;
                }
            }

            if (score > best_score) {
                range[0] = start;
                range[1] = i + 1;
                best_score = score;
            }

        }

        g = range;
        return g[1] - g[0];
    };
    
    if (group[1]-group[0] >= read2ctg_min_coverage && best_in_window(group) >= read2ctg_min_coverage ) {
        
        std::array<int, 2> aligned = get_aligned_length(links, group);
        if (((aligned[0] >= read2ctg_min_aligned_length) || (links[0]->source.len <= read2ctg_min_aligned_length && aligned[0] > 0.98*links[0]->source.len)) && 
            ((aligned[1] >= read2ctg_min_aligned_length) || (links[0]->target.len <= read2ctg_min_aligned_length && aligned[1] > 0.98*links[0]->target.len))) {
            return true;
        }

        return false;

    }  else {
        return false;
    }
}

bool ContigLink::CheckGroupDeviation(std::vector<C2r2cLink*> &links, std::array<int, 2> &group, int read2ctg_max_overhang, int ctg2ctg_max_overhang, int read2ctg_min_coverage, int read2ctg_min_aligned_length) {
    auto get_aligned_length = [](const std::vector<C2r2cLink*>& offset, const std::array<int, 2>& range) {
        std::array<int, 2> aligned = {0, 0};    
        for (auto i=range[0]; i< range[1]; ++i) {
            if (aligned[0] < offset[i]->source.end -  offset[i]->source.start)
                aligned[0] = offset[i]->source.end -  offset[i]->source.start;
            if (aligned[1] < offset[i]->target.end -  offset[i]->target.start)
                aligned[1] = offset[i]->target.end -  offset[i]->target.start;
        }
        return aligned;
    };

    auto mean = [&links](int s, int e) {
        int sum = 0;
        for (int i=s; i<e; ++i) sum += links[i]->pos_s2t[0];
        return sum/(e-s);
    };

    auto sd = [&links](int s, int e, int m) {
        int sum = 0;
        for (int i=s; i<e; ++i) sum += (links[i]->pos_s2t[0] - m)*(links[i]->pos_s2t[0] - m);
        return (int)sqrt(sum / (e-s));
    };

    auto confidence = [&links, mean, sd, read2ctg_min_coverage](std::array<int, 2> &g) {
        int m = mean(g[0], g[1]);
        int d = sd(g[0], g[1], m);

        for (; g[0]<g[1]; ++g[0]) {
            if (links[g[0]]->pos_s2t[0] >= m-3*d) break;
        }

        for (; g[1] > g[0]; --g[1]) {
            if (links[g[1]-1]->pos_s2t[0] <= m+3*d) break;
        }

        return g[1]-g[0] >= read2ctg_min_coverage;
    };
    
    if (group[1]-group[0] >= read2ctg_min_coverage && confidence(group)) {
        
        std::array<int, 2> aligned = get_aligned_length(links, group);
        if (((aligned[0] >= read2ctg_min_aligned_length) || (links[0]->source.len <= read2ctg_min_aligned_length && aligned[0] > 0.98*links[0]->source.len)) && 
            ((aligned[1] >= read2ctg_min_aligned_length) || (links[0]->target.len <= read2ctg_min_aligned_length && aligned[1] > 0.98*links[0]->target.len))) {
            return true;
        }

        return false;

    }  else {
        return false;
    }
}

bool ContigLink::CheckGroupCoverage(std::vector<C2r2cLink*>& links, const std::array<int,2>& g, int read2ctg_max_overhang, int read2ctg_min_coverage) {
    std::array<int, 2> area{0,0};
    for (int i=g[0]; i<g[1]; ++i) {
        auto link = links[i];
        area[0] += link->ol_expect[0];
        area[1] += link->ol_expect[1];
    }

    assert(g[0] < g[1]);
    
    area[0] /= g[1] - g[0];
    area[1] /= g[1] - g[0];

    if (area[1] - area[0] > 2*read2ctg_max_overhang) {
        std::vector<int> cov(area[1]-area[0]+1);

        for (int i=g[0]; i<g[1]; ++i) {
            auto link = links[i];
            int start = link->ol_real[0] - area[0];
            if (start < 0) start = 0;
            int end = link->ol_real[1] - area[0];
            if (end > area[1]-area[0]) end = area[1]-area[0];
            if (end > start) {
                cov[start] ++;
                cov[end] --;
            }
        }

        for (size_t i=1; i<cov.size(); ++i) {
            cov[i] += cov[i-1];
        }
        assert(cov.back() == 0);

        return *std::min_element(cov.begin()+read2ctg_max_overhang, cov.end()-read2ctg_max_overhang-1) > read2ctg_min_coverage;  


    } else {
        return true;
    }
}

Seq::Id ContigLink::Target() const {
    assert(c2c_links.size() > 0 || c2r2c_links.size() > 0);
    return c2c_links.size() > 0 ? c2c_links[0].target.id : c2r2c_links[0].target.id;
}

Seq::Id ContigLink::Source() const {
    assert(c2c_links.size() > 0 || c2r2c_links.size() > 0);
    return c2c_links.size() > 0 ? c2c_links[0].source.id : c2r2c_links[0].source.id;
}

} // namespace fsa {

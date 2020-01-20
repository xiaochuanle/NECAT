#include "overlap_compare.hpp"


#include <algorithm>
#include <cassert>

#include "logger.hpp"
#include "read_store.hpp"

namespace fsa {

OverlapCompare::OverlapCompare() {
    ol_stores_[0] = new OverlapStore(read_store_);
    ol_stores_[1] = new OverlapStore(read_store_);
}

OverlapCompare::~OverlapCompare() {
    delete ol_stores_[0];
    delete ol_stores_[1];
}

bool OverlapCompare::ParseArgument(int argc, const char *const argv[]) {
    ArgumentParser &&ap = GetArgumentParser();

    if (ap.ParseArgument(argc, argv)) {


    }
    else {
        return false;
    }
}

void OverlapCompare::Usage() {
    GetArgumentParser().Usage();
}

void OverlapCompare::Run() {

    LOG(INFO)("Start");
    PrintArguments();
    Load();


    std::vector<Stats> stats_list;
    for (const auto overlaps : ol_stores_) { 
        stats_list.push_back(Stats(*overlaps, max_overhang_));
    }

    Compare(stats_list[0], stats_list[1], "a_b");
    Compare(stats_list[1], stats_list[0], "b_a");

    LOG(INFO)("End");
}

void OverlapCompare::Load() {
    auto filter_simple = [&](Overlap& o) {
        return o.a_.len >= min_length_ && o.b_.len >= min_length_  &&
            o.AlignedLength() >= (size_t)min_aligned_length_ &&
            o.Location(max_overhang_) != Overlap::Loc::Abnormal;     
    };

    LOG(INFO)("Load file %s", fnames_[0].c_str()) ;
    ol_stores_[0]->Load(fnames_[0], "", (size_t)thread_size_, filter_simple);
    LOG(INFO)("Load file %s", fnames_[1].c_str());
    ol_stores_[1]->Load(fnames_[1], "", (size_t)thread_size_, filter_simple);
}
ArgumentParser OverlapCompare::GetArgumentParser() {
    ArgumentParser ap;
    ap.AddNamedOption(max_overhang_, "max_overhang", "max_overhang");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    ap.AddNamedOption(output_directory_, "output_directory", "Directory for storing output files");

    ap.AddPositionOption(fnames_[0], "fname0", "the first file name");
    ap.AddPositionOption(fnames_[1], "fname1", "the second file name");
    return ap;
}



void OverlapCompare::Compare(const Stats &a, const Stats &b, const std::string& name) {
    std::unordered_set<const Overlap*> a_not_in_b;
    std::unordered_set<const Overlap*> a_in_b_similar;
    std::unordered_set<const Overlap*> a_in_b_not_similar;

    for (const auto &o : a.goods) {
        const auto &it = b.goods.find(o.first);
        if (it != b.goods.end()) {
            if (IsSimilar(*o.second, *it->second)) {
                a_in_b_similar.insert(o.second);
            }
            else {
                a_in_b_not_similar.insert(o.second);
            }
        }
        else {
            a_not_in_b.insert(o.second);
        }
    }

    printf("A not in B: %zd\n", a_not_in_b.size());
    printf("A in B similar: %zd\n", a_in_b_similar.size());
    printf("A in B not similar: %zd\n", a_in_b_not_similar.size());

    a.overlaps.SaveM4aFile(name + "_not_in.m4", [&a_not_in_b](const Overlap&o) {
        return a_not_in_b.find(&o) != a_not_in_b.end(); });
    a.overlaps.SaveM4aFile(name + "_similar.m4", [&a_in_b_similar](const Overlap&o) {
        return a_in_b_similar.find(&o) != a_in_b_similar.end(); });
    a.overlaps.SaveM4aFile(name + "_no_similar.m4", [&a_in_b_not_similar](const Overlap&o) {
        return a_in_b_not_similar.find(&o) != a_in_b_not_similar.end(); });
}

bool OverlapCompare::IsSimilar(const Overlap &a, const Overlap &b) {


    if (a.a_.id == b.a_.id) {
        return (a.SameDirect() == b.SameDirect()) &&
            std::abs(a.a_.start - b.a_.start) <= max_overhang_ &&
            std::abs(a.a_.end - b.a_.end) <= max_overhang_ &&
            std::abs(a.b_.start - b.b_.start) <= max_overhang_ &&
            std::abs(a.b_.end - b.b_.end) <= max_overhang_;

    } else {
        return (a.SameDirect() == b.SameDirect()) &&
            std::abs(a.a_.start - b.b_.start) <= max_overhang_ &&
            std::abs(a.a_.end - b.b_.end) <= max_overhang_ &&
            std::abs(a.b_.start - b.a_.start) <= max_overhang_ &&
            std::abs(a.b_.end - b.a_.end) <= max_overhang_;
    }
}

bool OverlapCompare::Stats::IsNone(const Overlap &o) {
    // TODO
    return false;
    int a_start = o.a_.strand == 0 ? o.a_.start : o.a_.len - o.a_.end;
    int a_end = o.a_.strand == 0 ? o.a_.end : o.a_.len - o.a_.start;

    int b_start = o.b_.strand == 0 ? o.b_.start : o.b_.len - o.b_.end;
    int b_end = o.b_.strand == 0 ? o.b_.end : o.b_.len - o.b_.start;

    return ((a_start > max_overhang && b_start > max_overhang) ||
        (a_end < o.a_.len - max_overhang && b_end < o.b_.len - max_overhang));
}

bool OverlapCompare::Stats::IsContain(const Overlap &o) {
    // TODO
    return false;
    int a_start = o.a_.strand == 0 ? o.a_.start : o.a_.len - o.a_.end;
    int a_end = o.a_.strand == 0 ? o.a_.end : o.a_.len - o.a_.start;

    int b_start = o.b_.strand == 0 ? o.b_.start : o.b_.len - o.b_.end;
    int b_end = o.b_.strand == 0 ? o.b_.end : o.b_.len - o.b_.start;

    return (a_start < max_overhang && a_end > o.a_.len - max_overhang) || 
        (b_start < max_overhang && b_end > o.b_.len - max_overhang);
}

bool OverlapCompare::Stats::IsBetter(const Overlap &a, const Overlap &b) {
    return a.a_.end - a.a_.start > b.a_.end - b.a_.start;
}

OverlapCompare::Stats::Stats(const OverlapStore &ols, int ee): overlaps(ols), max_overhang(ee){
    printf("Size: %zd\n", overlaps.Size());

    for (const auto &o : overlaps.Get()) {
        uint64_t id = OverlapId(o);

        auto it = map.find(id);
        if (it == map.end()) {
            map[id] = std::vector<const Overlap*>({ &o });
        }
        else {
            map[id].push_back(&o);
        }
    }

    printf("Dup Size: %zd\n", overlaps.Size() - map.size());

    for (const auto &m : map) {
        const Overlap* best = nullptr;

        for (const auto &o : m.second) {
            if (IsNone(*o)) {
                nones.insert(o);
            }
            else if (IsContain(*o)) {
                contains.insert(o);
            }
            else {
                if (best == nullptr){
                    best = o;
                }
                else if (IsBetter(*o, *best)) {
                    dups.insert(best);
                    best = o;
                }
                else {
                    dups.insert(o);
                }
            }
        }
        if (best != nullptr) {
            goods[m.first] = best;
        }
    }

    printf("Contain Size: %zd\n", contains.size());
    printf("None Size: %zd\n", nones.size());
    printf("Dups Size: %zd\n", dups.size());
    printf("Good Size: %zd\n", goods.size());
}

void OverlapCompare::PrintArguments() {
    LOG(INFO)("Arguments: \n%s", GetArgumentParser().PrintOptions().c_str());
}

} // namespace fsa {
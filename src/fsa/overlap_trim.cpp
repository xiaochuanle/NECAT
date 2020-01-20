#include "overlap_trim.hpp"

#include <iostream>
#include "file_io.hpp"

namespace fsa {

bool OverlapTrim::ParseArgument(int argc, const char *const argv[]) {
    ArgumentParser &&ap = GetArgumentParser();

    if (ap.ParseArgument(argc, argv)) {

        return true;
    }
    else {
        return false;
    }
}

ArgumentParser OverlapTrim::GetArgumentParser() {
    ArgumentParser ap;

    ap.AddPositionOption(overlap_fname_, "overlap_fname", "overlap file path");
    //ap.AddPositionOption(read_name_, "read_name", "read name");
    return ap;
}

void OverlapTrim::Usage() {
    std::cout << GetArgumentParser().Usage();
}

void OverlapTrim::Run() {
    LoadOverlaps(overlap_fname_);

    GroupAndFilterDuplicate();
    FilterCoverageMt();
    DumpCoverage("./coverage_trim.txt.gz");
}

void OverlapTrim::LoadOverlaps(const std::string &fname) {
    auto filter_simple = [&](Overlap& o) {
        if (o.identity_ >= min_identity_ && o.a_.id != o.b_.id) {

            return true;
        } else {
            return false;
        }
    };

    ol_store_.Load(fname, "", (size_t)thread_size_, filter_simple);

    if (ol_store_.Size() > 0) {
        LOG(INFO)("Overlap size: %zd", ol_store_.Size());
    } else {
        LOG(FATAL)("No overlap was loaded");
    }
}


void OverlapTrim::GroupAndFilterDuplicate() {
    std::mutex mutex;

    auto add_overlap = [](int low, int a, int b, const Overlap& o, std::vector<std::unordered_map<Seq::Id, const Overlap*>>& group) {
        
        auto it = group[a-low].find(b);
        if (it == group[a-low].end()) {
            group[a-low][b] = &o;
        } else {
            if (o.AlignedLength() > it->second->AlignedLength()) {
                it->second = &o;
            }
        }

    };

    auto split_func = [this]() {
        auto r = ol_store_.GetReadIdRange();
        return SplitRange(thread_size_, r[0], r[1]);
    };
    auto comb_func = [this, &mutex](int low, std::vector<std::unordered_map<Seq::Id, const Overlap*>>&& group) {
        std::lock_guard<std::mutex> lock(mutex);
        for (size_t i=0; i<group.size(); ++i) {
            if (group[i].size() > 0) {
                groups_[low+(int)i] = std::move(group[i]);
            }

        }
    };

    auto work_func = [this, add_overlap, comb_func](std::array<Seq::Id, 2> r) {
        std::vector<std::unordered_map<Seq::Id, const Overlap*>> group(r[1] - r[0]);

        for (const auto &o : ol_store_.Get()) {
            if (o.a_.id >= r[0] && o.a_.id < r[1]) {
                add_overlap(r[0], o.a_.id, o.b_.id, o, group);
            }
            if (o.b_.id >= r[0] && o.b_.id < r[1]) {
                add_overlap(r[0], o.b_.id, o.a_.id, o, group);
            }
        }

        comb_func(r[0], std::move(group));
    };

    MultiThreadRun((int)thread_size_, split_func, work_func);

}


void OverlapTrim::FilterCoverageMt() {
    auto work_func = [&](const std::vector<int>& input) -> std::unordered_map<int, std::array<int, 2>> {
        std::unordered_map<int, std::array<int, 2>> output;
        
        for (auto i : input) {
            auto minmax = CalcMinMaxCoverage(i, groups_[i]);
            printf("minmax:, %d, %d\n", minmax.first, minmax.second);
            output.insert(std::make_pair(i, std::array<int,2>{minmax.first, minmax.second}));
        }
        return output;
    };

    coverages_ = MultiThreadRun(thread_size_, groups_, 
        SplitMapKeys<decltype(groups_)>, 
        work_func, 
        MoveCombineMapOrSet<std::unordered_map<int, std::array<int, 2>>>);   
    std::array<int, 3> coverage_params_ = CoverageParam1();
    printf("coverage_params_ %d, %d, %d\n", coverage_params_[0], coverage_params_[1], coverage_params_[2]);
}


std::pair<int, int> OverlapTrim::CalcMinMaxCoverage(int id, const std::unordered_map<int, const Overlap*>& group) {
    if (group.size() > 0) {
        std::vector<int> cov((group.begin()->second->a_.id == id ? group.begin()->second->a_.len :
            group.begin()->second->b_.len) + 1, 0);

        for (const auto &ig : group) {
            const Overlap& o = *ig.second;
            if (o.a_.id == id) {
                cov[o.a_.start] ++;
                cov[o.a_.end] --;
            } else {
                assert(o.b_.id == id);
                cov[o.b_.start] ++;
                cov[o.b_.end] --;
            }
        }
        for (size_t i = 1; i < cov.size(); ++i) {
            cov[i] += cov[i - 1];
        }
        assert(cov.back() == 0);

        auto c_minmax = std::minmax_element(cov.begin(), cov.end()-1);
        return std::make_pair(*c_minmax.first, *c_minmax.second);

    }
    else {
        return std::make_pair(0, 0);
    }


}

std::array<int, 3> OverlapTrim::CoverageParam1() const {
    std::array<int, 3> param { 0, 0, 0};

    std::vector<int> cov(coverages_.size());
    if (param[0] < 0 && coverages_.size() > 0) {
        std::transform(coverages_.begin(), coverages_.end(), cov.begin(), [](const decltype(coverages_)::value_type& a) {
            return a.second[0];
        });

        param[0] = FirstTrough(cov, 100, 9);
        //param[0] = Percentile(cov, coverage_percent_);
    }


    if (param[1] < 0 && coverages_.size() > 0) {
        std::transform(coverages_.begin(), coverages_.end(), cov.begin(), [](const decltype(coverages_)::value_type& a) {
            return a.second[1];
        });

        //param[1] = Percentile(cov, 100-coverage_discard_);
    }

    

    if (param[2] < 0 && coverages_.size() > 0) {
        std::transform(coverages_.begin(), coverages_.end(), cov.begin(), [](const decltype(coverages_)::value_type& a) {
            return a.second[1]-a.second[0];
        });

        //param[2] = Percentile(cov, 100-coverage_discard_);
    }

    return param;
}

int OverlapTrim::FirstTrough(const std::vector<int> &data, size_t last, size_t k) const {
    assert(k % 2 == 1 && data.size() >= k);

    auto minmaxv = std::minmax_element(data.begin(), data.end());

    auto minv = (*minmaxv.first);
    auto maxv = (*minmaxv.second);

    std::vector<int> counts(maxv-minv+1, 0);
    for (auto c : data) {
        counts[c-minv]++;
    }
    
    // calc the starting poistion
    size_t s = 0;
    for (size_t i=1; i<last/10; ++i) {
        if (counts[i-1] > counts[i]) {
            s = i - 1;
            break;
        }
    }

    int value = std::accumulate(counts.begin()+s, counts.begin()+s+k, 0);
    std::pair<size_t, int> best(s, value);
    for (size_t i=s+1; i<counts.size()+k-1 && i<last; ++i) {
        value += -counts[i-1] + counts[i+k-1];
        if (value < best.second*1.00) {
            best.first = i;
            best.second = value;
        } else {
            //best.first = i;
            //best.second = value;
            break;
        }
    }

    assert(best.first >= s);
    size_t bestbest = best.first;
    for (auto i=bestbest+1; i<best.first+k; ++i) {
        if (counts[i] < counts[bestbest]) bestbest = i;
    }

    double best_mean = double(best.second)  / k;
    for (auto i=best.first; i<best.first+k; ++i) {
        // if (counts[i] < counts[bestbest]) bestbest = i;
        if (counts[i] <= best_mean && counts[i] - counts[bestbest] <= counts[bestbest]*0.1) {
            bestbest = i; 
            break;
        }
    }
    return bestbest;
}

void OverlapTrim::DumpCoverage(const std::string &fname) const {
    gzFile of = gzopen(fname.c_str(), "w");
    if (of != nullptr) {
        for (auto c : coverages_) {
            //of << c.first << " " << c.second[0] << " " << c.second[1] << " " <<  c.second[1] -  c.second[0] << "\n";
            gzprintf(of, "%d %d %d %d\n", c.first, c.second[0], c.second[1], c.second[1] -  c.second[0]);
        }
        gzclose(of);
    } else {
        LOG(ERROR)("Fail to open coverage file %s", fname.c_str());
    }
}

} // namespace fsa {
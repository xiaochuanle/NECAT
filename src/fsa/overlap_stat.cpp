#include "overlap_stat.hpp"
#include "logger.hpp"

#include <numeric> 
#include <algorithm>
#include <list>
#include <iostream>
#include <cmath>

namespace fsa {

void test_load(const std::string &fname) {
    std::ifstream in(fname);
    std::string line;
    size_t N = 5000;
    bool r = true;
    while (r) {
        std::list<std::string> data;
        for (size_t i=0; i<N && r; ++i) {
            r = (bool)std::getline(in, line);
            if (r)
                data.push_back(line);
        }
    }
}
bool OverlapStat::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}

void Stat(std::vector<double>& identities) {
    auto mean = [&identities]() {
        double sum = std::accumulate(identities.begin(), identities.end(), 0.0);
        return sum / identities.size();
    };

    auto sd = [&identities](int m) {
        double sum = 0;
        for (auto i : identities) {
            sum += (i - m)*(i-m);
        }
        return sqrt(sum/identities.size());
    };

    auto MAD = [&identities]() {
        std::nth_element(identities.begin(), identities.begin() + identities.size()/2, identities.end());
        double m = identities[identities.size()/2];

        for (size_t i=0; i<identities.size(); ++i) {
            identities[i] = std::abs(identities[i]-m);
        }
        std::nth_element(identities.begin(), identities.begin() + identities.size()/2, identities.end());
        return identities[identities.size()/2];
    };

    double m = mean();
    double s = sd(m);
    double mm = MAD();
    printf("median: %f, mean: %f, sd: %f\n", mm, m, s);

}

void OverlapStat::Run() {
    std::mutex mutex;
    size_t block_size = 5000;
    std::vector<double> identites;
    std::vector<int> overhangs;
    std::vector<double> overhangs2;
    std::unordered_map<int, int> readlens;
    OverlapStore ol;

    std::vector<std::vector<double>> works(thread_size_);
    auto alloc_work = [&]() -> std::vector<double>& {
        std::lock_guard<std::mutex> lock(mutex);
        works.push_back(std::vector<double>());
        works.reserve(block_size);
        return works.back();

    };

    auto combine = [&](std::vector<double>& idens) {
        std::lock_guard<std::mutex> lock(mutex);
        identites.insert(identites.end(), idens.begin(), idens.end());
        idens.clear();
    };
    auto scan_overlap = [&](Overlap& o) {
        std::vector<double> thread_local &idens = alloc_work();
        if (o.identity_ > 0) idens.push_back(o.identity_);
        
        if (idens.size() >= block_size) {
            combine(idens);
        }
        return false;
    };


    ol.Load(ifname_, "", thread_size_, scan_overlap);
    for (auto & idens : works) {
        combine(idens);
    }

    Stat(identites);

    //for (auto i : identites) {
    //    printf("%.02f\n", i);
    //}
    //for (const auto& i : readlens) {
    //    printf("%i\n", i.second);
    //}
    //for (const auto i : overhangs) {
    //    printf("%d\n", i);
    //}
    //for (const auto i : overhangs2) {
    //    printf("%g\n", i);
    //}
}

void OverlapStat::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser OverlapStat::GetArgumentParser() {
    ArgumentParser ap;
    ap.AddPositionOption(ifname_, "ifname", "overlap file name");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    return ap;
}


void OverlapStat::Load(const std::string &fname, const std::string &type) {

    ol_store_.Load(fname, type, 1, [](Overlap&o) {return false;});
}

void OverlapStat::StatReads() {
    for (auto &o : ol_store_.Get()) {
        reads_[o.a_.id].len = o.a_.len;
        reads_[o.b_.id].len = o.b_.len;
    }

    printf("Read count: %zd", reads_.size());

    std::vector<int> read_lens(reads_.size());
    auto i = reads_.begin();
    std::generate(read_lens.begin(), read_lens.end(), [&i]() {
        return (i++)->second.len;
    });

    int total_length = std::accumulate(read_lens.begin(), read_lens.end(), 0);
    printf("Read total length: %d\n", total_length);
    printf("Read avarage length: %d\n", total_length / (int)read_lens.size());

    
}

void OverlapStat::FindCommon(const std::string &fname, const std::string &read0, const std::string &read1) {
    OverlapStore ol_store;

    Seq::Id id0 =  ol_store.GetReadId(read0);
    Seq::Id id1 = ol_store.GetReadId(read1);

    auto filter = [id0, id1](Overlap& o) ->bool { 
        return o.a_.id == id0 || o.b_.id == id0 || o.a_.id == id1 || o.b_.id == id1;
    };

    ol_store.Load(fname, "m4a", 48, filter);
    

    std::unordered_map<int, std::vector<const Overlap*>> set0;
    std::unordered_map<int, std::vector<const Overlap*>> set1;

    for (auto& o : ol_store.Get()) {
        if (o.a_.id == id0) {
           set0[o.b_.id].push_back(&o);
        }    
        if (o.b_.id == id0) {
           set0[o.a_.id].push_back(&o);
        }    
        if (o.a_.id == id1) {
           set1[o.b_.id].push_back(&o);
        }    
        if (o.b_.id == id1) {
           set1[o.a_.id].push_back(&o);
        }    
    }

    if (set0.find(id1) != set0.end()) {
        for (auto o : set0[id1]) {
            printf("%s", ol_store.ToM4aLine(*o).c_str());
        }
        printf("\n");
    }

    for (auto &i : set0) {
        if (set1.find(i.first) != set1.end()) {
            for (auto o : i.second) {
                printf("%s", ol_store.ToM4aLine(*o).c_str());
            }
            for (auto o : set1[i.first]) {
                printf("%s", ol_store.ToM4aLine(*o).c_str());
            }
        }
    }
}

} // namespace fsa {

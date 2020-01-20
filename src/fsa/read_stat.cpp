#include "read_stat.hpp"

#include <numeric> 
#include <algorithm>
#include <list>
#include <iostream>

#include "logger.hpp"
#include "read_store.hpp"

namespace fsa {

bool ReadStat::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}


void ReadStat::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser ReadStat::GetArgumentParser() {
    ArgumentParser ap;
    //ap.AddPositionOption(action_, "action", "");
    ap.AddPositionOption(ifname_, "ifname", "read file name");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    return ap;
}

void ReadStat::Run() {
    LOG(INFO)("Start");

    LOG(INFO)("Load read file %s", ifname_.c_str());
    ReadStore rs;
    rs.Load(ifname_, "", 4);

    if (action_ == "N50") {
        StatN50(rs);
    } else {
        LOG(WARNING)("Unrecognize action %s", action_.c_str());
    }
    LOG(INFO)("END");
}

void ReadStat::StatN50(const ReadStore &rs) {
    std::vector<int> read_lens(rs.GetIdRange()[1]);
    for (size_t i=0; i<read_lens.size(); ++i) {
        read_lens[i] = rs.GetSeqLength(i);
    }

    std::sort(read_lens.begin(), read_lens.end(), [](int a, int b) { return a > b; });

    long long total_length = std::accumulate(read_lens.begin(), read_lens.end(), (long long)0);

    std::cout << "Count: " << read_lens.size() << "\n";
    std::cout << "Total: " << total_length << "\n";
    std::cout << "Max: " << read_lens.front() << "\n";
    std::cout << "Min: " << read_lens.back() << "\n";
    
    long long accu = 0;
    int ns[] = { 25, 50, 75};
    size_t ins = 0;
    for (size_t i=0; i<read_lens.size(); ++i) {
        accu += read_lens[i];
        for (; ins < sizeof(ns) / sizeof(ns[0]) && accu > (long long)ns[ins]*1.0/100 * total_length; ++ins) {
            std::cout << "N" << ns[ins] << ": " << read_lens[i] << "\n";
            std::cout << "L" << ns[ins] << ": " << i+1 << "\n";
        }
        if (ins >= sizeof(ns) / sizeof(ns[0])) {
            break;
        }
    }
}



} // namespace fsa {
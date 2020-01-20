#include "overlap_improve.hpp"

#include <cassert>
#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_set>
#include <thread>

#include "read_store.hpp"
#include "overlap_store.hpp"
#include "logger.hpp"
#include "utility.hpp"


bool OverlapImprover::ParseArgument(int argc, const char *const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}

void OverlapImprover::Usage() {
    std::cout << GetArgumentParser().Usage();
}

void OverlapImprover::Run() {

    LOG(FATAL)("The procedure is being developed");
    LOG(INFO)("Start");


    ReadStore read_store;

    LOG(INFO)("Load fasta file, %s", read_file_.c_str());
    read_store.Load(read_file_);

    OverlapStore ol_store(read_store);

    LOG(INFO)("Load rough overlap file, %s",ifname_.c_str());
    ol_store.LoadM4File(ifname_);

    Improve(ol_store, read_store);



    LOG(INFO)("Save improved overlaps file, %s", ofname_.c_str());
    ol_store.SaveM4File(ofname_, [](const Overlap&){return true;});
    LOG(INFO)("End");
}

void OverlapImprover::Improve(OverlapStore &ol_store, ReadStore &read_store) {

    auto work_func = [&](const std::array<size_t, 2>& input) {


        LOG(INFO)("Start Worker Task: %d - %d", (int)input[0], (int)input[1]);
        for (size_t i = input[0]; i < input[1]; ++i) {
            Overlap &o = ol_store.Get(i);
    
            assert(!"TODO code modified");
            if (!Filter(o)) {
            }
            else {
                assert(!"TODO");
            }
        }
        LOG(INFO)("End Worker Task: %d - %d", (int)input[0], (int)input[1]);
        return nullptr;
    };

    MultiThreadRun(thread_size_, ol_store.Get(), SplitVectorKeys<decltype(ol_store.Get())>, work_func); 
}

bool OverlapImprover::Filter(Overlap &o) {

    int a_start = o.a_.strand == 0 ? o.a_.start : o.a_.len - o.a_.end;
    int a_end = o.a_.strand == 0 ? o.a_.end : o.a_.len - o.a_.start;

    int b_start = o.b_.strand == 0 ? o.b_.start : o.b_.len - o.b_.end;
    int b_end = o.b_.strand == 0 ? o.b_.end : o.b_.len - o.b_.start;

    return (a_start > max_overhang_ && b_start > max_overhang_) ||
        (a_end < o.a_.len - max_overhang_ && b_end < o.b_.len - max_overhang_);
}

ArgumentParser OverlapImprover::GetArgumentParser() {
    ArgumentParser ap;
    ap.AddNamedOption(max_overhang_, "max_overhang", "end error", "INT");
    ap.AddNamedOption(thread_size_, "thread_size", "thread size", "INT");
    ap.AddNamedOption(read_file_, "read_file", "read file path");

    ap.AddPositionOption(ifname_, "ifname", "Input file name");
    ap.AddPositionOption(ifname_, "ofname", "Output file name");

    return ap;
}


int main(int argc, const char *const argv[])
{
    OverlapImprover overlap_improver;

    if (overlap_improver.ParseArgument(argc, argv)) {
        overlap_improver.Run();
    }
    else {
        overlap_improver.Usage();
    }
    return 0;

}
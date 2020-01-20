#include "overlap_show.hpp"

#include <cassert>
#include <algorithm>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include <sstream>
#include <iostream>

#include "utility.hpp"

OverlapShow::OverlapShow() {
   
}


bool OverlapShow::ParseArgument(int argc, const char *const argv[]) {
    ArgumentParser &&ap = GetArgumentParser();

    if (ap.ParseArgument(argc, argv)) {


//       type_ = ol_store_.DetectFileType(ifname_);
        return true;
    }
    else {
        return false;
    }
}

ArgumentParser OverlapShow::GetArgumentParser() {
    ArgumentParser ap;

    ap.AddPositionOption(ifname_, "ifname", "overlap file path");
    ap.AddPositionOption(read_name_, "read_name", "read name");
    return ap;
}

void OverlapShow::Usage() {
    std::cout << GetArgumentParser().Usage();
}


void OverlapShow::Run() {
//    LOG(INFO)("Load overlap file: %s", ifname_.c_str());
//    LoadOverlaps(ifname_);
 
    std::unordered_set<std::string> names;
    GzFileReader reader(read_name_);
    std::string n = reader.GetStrippedLine();
    while (!n.empty()) {
        names.insert(n);
        n = reader.GetStrippedLine();
    }
    LOG(INFO)("dup name size: %zd", names.size());

    GzFileReader reader2(ifname_);
    std::string line = reader2.GetLine();
    while (!line.empty()) {
        auto its = SplitStringBySpace(line);
        if (names.find(its[0]) == names.end() && names.find(its[1]) == names.end()) {
            printf("%s", line.c_str());
        }
        line = reader2.GetLine();
    }

    LOG(INFO)("Load overlap file: %s", ifname_.c_str());
}

void OverlapShow::LoadOverlaps(const std::string &fname) {

    ol_store_.Load(ifname_);
    

}

std::string OverlapShow::ToString(const Overlap &o) const {
    auto t = ol_store_.DetectFileType(ifname_);
    if (t == "m4") return ol_store_.ToM4Line(o);
    else if (t == "m4a") return ol_store_.ToM4aLine(o);
    else if (t == "paf") return ol_store_.ToPafLine(o);
    else return ol_store_.ToM4Line(o);
}

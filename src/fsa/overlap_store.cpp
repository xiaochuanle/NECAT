#include "overlap_store.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cassert>

#include "read_store.hpp"
#include "logger.hpp"


namespace fsa {


std::vector<std::string> SplitString(const std::string &str, const std::string sep) {
    std::vector<std::string> substrs;

    std::string::size_type begin = str.find_first_not_of(sep);
    while (begin != std::string::npos) {
        std::string::size_type end = str.find(sep, begin);
        if (end != std::string::npos) {
            substrs.push_back(str.substr(begin, end - begin));
            begin = end + sep.length();
        }
        else {
            substrs.push_back(str.substr(begin));
            begin = end;
        }
    }
    return substrs;
}

std::string OverlapStore::DetectFileType(const std::string &fname) {
    if (fname.size() >= 3 && fname.substr(fname.size()-3) == ".m4") {
        return "m4";
    } else if (fname.size() >= 6 && fname.substr(fname.size()-6) == ".m4.gz") {
        return "m4.gz";
    } else if (fname.size() >= 4 && fname.substr(fname.size()-4) == ".m4a") {
        return "m4a";
    } else if (fname.size() >= 7 && fname.substr(fname.size()-7) == ".m4a.gz") {
        return "m4a.gz";
    } else if (fname.size() >= 4 && fname.substr(fname.size()-4) == ".paf") {
        return "paf";
    } else if (fname.size() >= 7 && fname.substr(fname.size()-7) == ".paf.gz") {
        return "paf.gz";
    } else if (fname.size() >= 4 && fname.substr(fname.size()-4) == ".ovl") {
        return "ovl";
    } else if (fname.size() >= 7 && fname.substr(fname.size()-7) == ".ovl.gz") {
        return "ovl.gz";
    } else {
        auto i = fname.find_last_of('.');
        return fname.substr(i == fname.npos ? 0 : i+1);
    }
}

bool OverlapStore::FromM4Line(const std::string &line, Overlap& o) {

    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 12) {

        // M4文件的Id就是read在fasta文件的序号。
        o.a_.id = atoi(items[0].c_str()) - 1;
        o.b_.id = atoi(items[1].c_str()) - 1;

        o.identity_ = atof(items[2].c_str());
        o.score_ = atoi(items[3].c_str());

        o.a_.strand = atoi(items[4].c_str());
        o.a_.start = atoi(items[5].c_str());
        o.a_.end = atoi(items[6].c_str());
        o.a_.len = atoi(items[7].c_str());

        o.b_.strand = atoi(items[8].c_str());
        o.b_.start = atoi(items[9].c_str());
        o.b_.end = atoi(items[10].c_str());
        o.b_.len = atoi(items[11].c_str());

        // 调整strand，保证a.strand = 0, 先设置b，再设置a
        o.b_.strand = o.a_.strand == o.b_.strand ? 0 : 1;
        o.a_.strand = 0;

        o.score_ = -((o.a_.end - o.a_.start) + (o.b_.end - o.b_.start)) / 2;

        return true;
    }
    else {
        return false;
    }
}


bool OverlapStore::FromM4aLine(const std::string &line, Overlap& o) {

    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 12) {

        o.a_.id = read_store_.GetIdByNameSafe(items[0]);
        o.b_.id = read_store_.GetIdByNameSafe(items[1]);

        o.identity_ = atof(items[2].c_str());
        o.score_ = atoi(items[3].c_str());

        o.a_.strand = atoi(items[4].c_str());
        o.a_.start = atoi(items[5].c_str());
        o.a_.end = atoi(items[6].c_str());
        o.a_.len = atoi(items[7].c_str());

        o.b_.strand = atoi(items[8].c_str());
        o.b_.start = atoi(items[9].c_str());
        o.b_.end = atoi(items[10].c_str());
        o.b_.len = atoi(items[11].c_str());

        o.score_ = -((o.a_.end - o.a_.start) + (o.b_.end - o.b_.start)) / 2;

        return true;
    }
    else {
        return false;
    }
}

bool OverlapStore::FromOvlLine(const std::string &line, Overlap& o) {
    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 13) {

        o.a_.id = atoi(items[0].c_str());
        o.b_.id = atoi(items[1].c_str());
        //o.a_.id = read_store_.NameToId(items[0]);
        //o.b_.id = read_store_.NameToId(items[1]);

        o.score_ = atoi(items[2].c_str());
        o.identity_ = atof(items[3].c_str());

        o.a_.strand = atoi(items[4].c_str());
        o.a_.start = atoi(items[5].c_str());
        o.a_.end = atoi(items[6].c_str());
        o.a_.len = atoi(items[7].c_str());

        o.b_.strand = atoi(items[8].c_str());
        o.b_.start = atoi(items[9].c_str());
        o.b_.end = atoi(items[10].c_str());
        o.b_.len = atoi(items[11].c_str());

        return true;
    } else {
        return false;
    }
}


bool OverlapStore::FromPafLine(const std::string &line, Overlap& o) {
    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 12) {
        // query_name, query_length, query_start, query_end, 
        // relative_strand, 
        // target_name, target_lenght, target_start, target_end, 
        // number_residue_matches, alignment_block_length, mapping_quality


        // query_name, query_length, query_start, query_end, 
        o.a_.id = read_store_.GetIdByNameSafe(items[0]);
        o.a_.len = atoi(items[1].c_str());
        o.a_.start = atoi(items[2].c_str());
        o.a_.end = atoi(items[3].c_str());
        
        // relative_strand, 
        o.a_.strand = items[4] == "+" ? 0 : 1;
        o.b_.strand = 0;
        
        // target_name, target_lenght, target_start, target_end, 
        o.b_.id = read_store_.GetIdByNameSafe(items[5]);
        o.b_.len = atoi(items[6].c_str());
        o.b_.start = atoi(items[7].c_str());
        o.b_.end = atoi(items[8].c_str());

        // number_residue_matches, alignment_block_length, mapping_quality
        
        o.identity_ = atof(items[9].c_str())*100 / atoi(items[10].c_str());
        o.score_ = -atoi(items[10].c_str());
        // items[11]

        return true;
    } else {
        return false;
    }
}

 std::array<Seq::Id, 2> OverlapStore::GetReadIdRange() const { 
    auto range = read_store_.GetIdRange(); 
    if (range[1] == 0) {
        auto cmp_range = [&range](Seq::Id i) {
            if (i < range[0]) range[0] = i;
            else if (i >= range[1]) range[1] = i+1;
        };
        for (auto &o : overlaps_) {
            cmp_range(o.a_.id);
            cmp_range(o.b_.id);
            
        }
    }

    return range;
}

std::unordered_map<int, std::unordered_map<int, Overlap*>> OverlapStore::Group() {
    std::unordered_map<int, std::unordered_map<int, Overlap*>> groups;

    for (auto &o : overlaps_) {

        auto add_overlap = [&](int a, int b, Overlap& o) {
            auto iter = groups.find(a);
            if (iter != groups.end()) {
                assert(iter->second.find(b) == iter->second.end()); // 已经删除了重复的overlap
                iter->second[b] = &o;
            }
            else {
                auto iter = groups.insert(std::make_pair(a, std::unordered_map<int, Overlap*>()));
                iter.first->second[b] = &o;
            }
        };

        add_overlap(o.a_.id, o.b_.id, o);
        add_overlap(o.b_.id, o.a_.id, o);
        
    }

    return groups;
}

std::unordered_map<int, std::unordered_map<int, const Overlap*>> OverlapStore::Group(bool(*better)(const Overlap& a, const Overlap &b)) const {
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> groups;

    for (auto &o : overlaps_) {

        auto add_overlap = [&](int a, int b, const Overlap& o) {
            auto iter = groups.find(a);
            if (iter != groups.end()) {
                auto itit = iter->second.find(b);
                if (itit == iter->second.end()) {
                    iter->second.insert(std::make_pair(b, &o));
                } else {
                    if (better(o, *itit->second)) {
                        itit->second = &o;
                    }
                }
            }
            else {
                auto iter = groups.insert(std::make_pair(a, std::unordered_map<int, const Overlap*>()));
                iter.first->second[b] = &o;
            }
        };

        add_overlap(o.a_.id, o.b_.id, o);
        add_overlap(o.b_.id, o.a_.id, o);
        
    }

    return groups;
}

std::unordered_map<int, std::unordered_map<int, Overlap*>> OverlapStore::GroupTarget(bool(*better)(const Overlap* a, const Overlap *b)) {
    std::unordered_map<int, std::unordered_map<int, Overlap*>> groups;

    for (auto &o : overlaps_) {
        auto add_overlap = [&](int a, int b, Overlap& o) {
            auto iter = groups.find(a);
            if (iter != groups.end()) {
                auto biter = iter->second.find(b);
                if (biter == iter->second.end()) {
                    iter->second[b] = &o;
                }
                else {
                    if (!better(biter->second, &o)) {
                        iter->second[b] = &o;
                    }
                }
            }
            else {
                auto iter = groups.insert(std::make_pair(a, std::unordered_map<int, Overlap*>()));
                iter.first->second[b] = &o;
            }
        };

        add_overlap(o.b_.id, o.a_.id, o);
    
    }

    return groups;
}

std::string OverlapStore::ToM4aLine(const Overlap& o) const{

    std::ostringstream oss;

    oss << read_store_.IdToName(o.a_.id) << " " << read_store_.IdToName(o.b_.id) << " " 
        << o.identity_ << " " << o.score_ << " "
        << o.a_.strand << " " << o.a_.start << " " << o.a_.end << " " << o.a_.len << " "
        << o.b_.strand << " " << o.b_.start << " " << o.b_.end << " " << o.b_.len << "\n";

    return oss.str();
}

std::string OverlapStore::ToM4Line(const Overlap &o) const {
    return o.ToM4Line();
}

std::string OverlapStore::ToPafLine(const Overlap &o) const {
    return ToM4aLine(o);
}

} // namespace fsa {


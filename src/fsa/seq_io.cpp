#include "seq_io.hpp"

#include <cassert>
#include <algorithm>
#include <cassert>

namespace fsa {

bool FastaReader::Next(Item &item) {
    item.id = Tell();
    item.quality = "";
    return GetHead(item.head, item.sub_head) && GetSeq(item.seq);
}

bool FastaReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = in_.QueryStrippedLine();
    in_.ConsumeStrippedLine();
    //printf("xxx：%s\n", line.c_str()); fflush(stdout);
    if (!line.empty() && line[0] == '>') {
        std::string::size_type s = 1;
        while (s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s+1, line.size());
        while (e < line.size() && !::isspace(line[e])) e++;

        if (e > s) {
            head = line.substr(s, e-s);
            while (e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool FastaReader::GetSeq(std::string &seq) {
    seq = "";
    std::string line = in_.QueryStrippedLine();

    while (!line.empty() && line[0] != '>') {
        seq += line;
        in_.ConsumeStrippedLine();
        line = in_.QueryStrippedLine();
    }

    return true;// !seq.empty(); allow empty sequence

}



bool FastqReader::Next(Item &item) {
    assert(IsValid());

    item.id = Tell();
    return GetHead(item.head, item.sub_head) && GetSeq(item.seq) &&
           GetHead1() && GetQuality(item.quality);
}

bool FastqReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = in_.GetStrippedLine();

    if (!line.empty() && line[0] == '@') {
        std::string::size_type s = 1;
        while (s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s+1, line.size());
        while (e < line.size() && !::isspace(line[e])) e++;

        if (e > s) {
            head = line.substr(s, e-s);
            while (e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

} // namespace fsa {
    
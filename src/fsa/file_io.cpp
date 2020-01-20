#include "file_io.hpp"

#include <cassert>
namespace fsa {
/*
std::string SeqReader::NextNonEmptyLine(std::ifstream &in) {
    std::string line;
    bool r = (bool)std::getline(in, line);
    while (r) {
        std::string::size_type s = 0;
        while (s< line.size() && ::isspace(line[s])) s++;
        
        std::string::size_type e = line.size();
        while (e > s && ::isspace(line[e-1])) e--;
        
        if (e > s) return line.substr(s, e-s);
        
        r = (bool)std::getline(in, line);
    }
    return std::string();
}
*/

std::string StripString(const std::string &line) {
    
    std::string::size_type s = 0;
    while (s< line.size() && ::isspace(line[s])) s++;
        
    std::string::size_type e = line.size();
    while (e > s && ::isspace(line[e-1])) e--;
        
    return line.substr(s, e-s);
        
}


GzFileReader::GzFileReader(const std::string &fname) {
    in_ = gzopen(fname.c_str(), "rb");
}

std::string GzFileReader::GetLine() {
    std::string str;

    const int BUF_SIZE = 1024*2;
    char buffer[BUF_SIZE];

    char* line = gzgets(in_, buffer, BUF_SIZE);
    while (line != nullptr) {
        str += line;
        if (str.back() == '\n') {
            break;
        }
        line = gzgets(in_, buffer, BUF_SIZE);
    }

    return str;
}

bool GzFileReader::GetLine(std::string &str) {
    str.clear();

    const int BUF_SIZE = 1024*2;
    char buffer[BUF_SIZE];

    char* line = gzgets(in_, buffer, BUF_SIZE);
    while (line != nullptr) {
        str += line;
        if (str.back() == '\n') {
            break;
        }
        line = gzgets(in_, buffer, BUF_SIZE);
    }

    return !str.empty();
}

size_t GzFileReader::GetLines(std::vector<std::string> &lines) {

    const int BUF_SIZE = 1024*2;
    char buffer[BUF_SIZE];

    size_t index = 0; 
    while (index < lines.size()) {
        std::string &str = lines[index];
        str.clear();

        char* line = gzgets(in_, buffer, BUF_SIZE);
        while (line != nullptr) {
            str += line;
            if (str.back() == '\n') {
                break;
            }
            line = gzgets(in_, buffer, BUF_SIZE);
        }

        if (str.empty()) break;

        index ++;
    }

    return index;
}


std::string GzFileReader::GetStrippedLine() {
    if (!HasStrippedLine()) {
        std::string str = StripString(GetLine());
        while (str.empty() && !IsEnd()) {
            str = StripString(GetLine());
        } 
        return str;
    } else {
        next_line_pos_ = -1;
        return next_line_;
    }
}


} // namespace fsa {
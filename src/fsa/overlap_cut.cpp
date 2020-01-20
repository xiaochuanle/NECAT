#include "overlap_cut.hpp"
#include "logger.hpp"

#include <numeric> 
#include <algorithm>
#include <list>
#include <iostream>
#include <fstream>
#include <cmath>

namespace fsa {
bool OverlapCut::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}

void OverlapCut::Run() {

    std::ofstream of(local_identity_fname_);

    std::string line;
    while(getline(std::cin, line)) {
        std::string nline = CutPafLine(line, of);
        //std::string nline = line;
        if (!nline.empty())  std::cout << nline << "\n";
    }
}

void OverlapCut::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser OverlapCut::GetArgumentParser() {
    ArgumentParser ap;
    ap.AddPositionOption(ifname_, "ifname", "overlap file name");
    ap.AddPositionOption(ofname_, "ofname", "overlap file name");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    ap.AddNamedOption(read_file_, "read_file", "read file");
    ap.AddNamedOption(local_identity_fname_, "local_identity_fname", "");
    return ap;
}


void OverlapCut::Load(const std::string &fname, const std::string &type) {

    ol_store_.Load(fname, type, 1, [](Overlap&o) {return false;});
}

bool OverlapCut::IterCigar::Next(size_t &n, char &t) {
    n = 0;
    t = 0;

    if (index < cigar.size()) {
        size_t start = index;
        for (; index<cigar.size(); ++index) {
            if (!::isdigit(cigar[index])) break;
        }

        n = stoul(cigar.substr(start, index-start));
        
        //n = stoul(cigar, &index);
        //LOG(INFO)("INDEX %zd", index);
        if (index < cigar.size()) {
            t = cigar[index++];
        }

    }
    return n != 0 && t != 0;
}

void OverlapCut::ComputeIdentity(const std::string& qseq, size_t qs, size_t qe, const std::string &tseq, size_t ts, size_t te, const std::string &cigar) {
    size_t cigar_num = 0; 
    char cigar_type = 0;


    size_t qi = qs, ti = ts;
    size_t mat = 0, ins = 0, del = 0, mismat = 0;
    std::vector<std::array<size_t, 2>> all_mats {{ti, mat}};

    IterCigar iter_cigar(cigar);
    while (iter_cigar.Next(cigar_num, cigar_type)) {
        switch (cigar_type) {
        case 'M':
            for (size_t i=0; i<cigar_num; ++i) {
                if (tseq[ti+i] == qseq[qi+i]) {
                    mat += 1;
                } else {
                    mismat += 1;
                }
            }
            ti += cigar_num;
            qi += cigar_num;
            break;

        case '=':
            mat += cigar_num;
            qi += cigar_num;
            ti += cigar_num;
            break;

        case 'X':
            mismat += cigar_num;
            qi += cigar_num;
            ti += cigar_num;
            break;

        case 'I':
            ins += cigar_num;
            qi += cigar_num; 
            break;

        case 'D':
            del += cigar_num;
            ti += cigar_num; 
            break;

        case 'S':
        case 'H':
        default:
            assert("never come here");
            break;
        }
        all_mats.push_back({ti, mat});
    }

    
}

std::array<size_t,2> OverlapCut::ComputeIdentity(size_t qs, size_t qe, size_t ts, size_t te, const std::string &cigar) {
    size_t cigar_num = 0; 
    char cigar_type = 0;


    size_t qi = qs, ti = ts;
    size_t mat = 0, ins = 0, del = 0, mismat = 0;

    
    std::vector<size_t> errors(te-ts, 0);

    IterCigar iter_cigar(cigar);
    while (iter_cigar.Next(cigar_num, cigar_type)) {
        switch (cigar_type) {
        case 'M':
            assert("never come here");  // M会转换成X和=
            break;

        case '=':
            mat += cigar_num;
            qi += cigar_num;
            ti += cigar_num;
            
            break;

        case 'X':
            std::fill(errors.begin()+ti-ts, errors.begin()+ti-ts+cigar_num, 1);
            mismat += cigar_num;
            qi += cigar_num;
            ti += cigar_num;
            break;

        case 'I':
            errors[ti-ts] += cigar_num;
            ins += cigar_num;
            qi += cigar_num; 
            break;

        case 'D':
            std::fill(errors.begin()+ti-ts, errors.begin()+ti-ts+cigar_num, 1);
            del += cigar_num;
            ti += cigar_num; 
            break;

        case 'S':
        case 'H':
        default:
            assert("never come here");
            break;
        }
    }

    
    assert(errors.size() >= window_size_);
    size_t cur_local_error = std::accumulate(errors.begin(), errors.begin()+window_size_, 0);

    size_t pos_max_local_error = 0;
    size_t max_local_error = cur_local_error;
    for (size_t i=window_size_; i<errors.size(); ++i) {
        cur_local_error += errors[i];
        cur_local_error -= errors[i-window_size_];

        if (cur_local_error > max_local_error) {
            max_local_error = cur_local_error;
            pos_max_local_error = i + 1 - max_local_error;
        }
    }

    return {max_local_error, pos_max_local_error};
}


std::string OverlapCut::CutPafLine(const std::string &line, std::ofstream &of) {
    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 12) {
        size_t qs = stoll(items[2]);
        size_t qe = stoll(items[3]);
        size_t ts = stoll(items[7]);
        size_t te = stoll(items[8]);

        if (te - ts >= (size_t)min_aligned_length_) {
            std::string cigar;
            for (size_t i=12; i<items.size(); ++i) {
                if (items[i].size() >= 5 && items[i].compare(0, 5, "cg:Z:") == 0) {
                    cigar = items[i].substr(5);
                }
            }

            auto iden = ComputeIdentity(qs, qe, ts, te, cigar);
            of << items[0] << "\t" << items[5] << "\t" << iden[0]*1.0/window_size_ << "\t" <<iden[1] << "\n";

            std::string newline = items[0];
            for (size_t i= 1; i<12; ++i) {
                newline.append("\t");
                newline.append(items[i]);
            }
            return newline;
        } else {
            return "";
        }
    } else {
        return "";
    }
}

} // namespace fsa {

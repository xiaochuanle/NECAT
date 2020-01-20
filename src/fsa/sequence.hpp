#ifndef FSA_SEQUENCE_HPP
#define FSA_SEQUENCE_HPP

#include <cassert>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>   

namespace fsa {

struct Seq {
    typedef int Id;
    typedef int EndId;
    static EndId IdToEndId(Id id, int end) {
        assert(id >= 0 && (end ==0 || end == 1));
        return end == 0 ? id + 1 : -id - 1;
    }


    static Id EndIdToId(EndId id) {
        assert(id != 0);
        return id > 0 ? id - 1 : (-id) - 1;
    }

    static EndId ReverseEndId(EndId id) { return -id; }

    struct Area {
        Seq::Id id;
        int strand;
        int start;
        int end;
    };

    static std::string ReverseComplement(const std::string &seq);
};


class DnaSerialTable {
public:
    DnaSerialTable() { 
        std::fill(table_, table_+256, -1); 
        table_['A'] = table_['a'] = 0;
        table_['C'] = table_['c'] = 1;
        table_['G'] = table_['g'] = 2;
        table_['T'] = table_['t'] = 3;
    }
    uint8_t operator[] (char c) const {
        return table_[(uint8_t)c];
    }

protected:
    uint8_t table_[256];
};

class DnaComplementTable {

public:
    DnaComplementTable() {
        table['A'] = 'T';
        table['C'] = 'G';
        table['G'] = 'C';
        table['T'] = 'A';
        table['a'] = 't';
        table['c'] = 'g';
        table['g'] = 'c';
        table['t'] = 'a';
        table['-'] = '-';
    }

    char operator [](int c) { return table[c]; }
    char table[256] = {0};
};


class DnaSeq {
public:
    DnaSeq() : len_(0) {}
    DnaSeq(const std::string& str) : len_(str.size()), data_((str.size()+3)/4, 0) { Reset(str);  }

    void Reset(const std::string &str) {
        
        len_ = str.size();
        data_.assign((str.size()+3)/4, 0);
        size_t index = 0;
        for (index=0; index+3<str.size(); index+=4) {
            assert(index/4 < data_.size());
            // assert(s_SerialTable[str[index+3]] != (uint8_t) -1);

            data_[index/4] = (s_SerialTable[str[index+3]] << 6) + (s_SerialTable[str[index+2]] << 4) + 
                            (s_SerialTable[str[index+1]] << 2) + (s_SerialTable[str[index]]);

   
        }
        if (index < str.size()) {
        assert(index/4 < data_.size());
            data_[index/4] = 0;
            for (size_t i=index; i<str.size(); ++i) {
                data_[index/4] += s_SerialTable[str[i]] << (i-index)*2;
            }
        }
        
    }

    size_t Size() const { return len_;}
    
    uint8_t operator [](size_t i) const {
        assert(i < len_);
        return (data_[i/4] >> ((i%4)*2)) & 0x3;
    }
    std::shared_ptr<std::string> ToString(bool upper=true) const {
        std::shared_ptr<std::string> str(new std::string(len_, 'N'));
        const char* base = upper ? "ACGT" : "acgt";
        for (size_t i=0; i<len_; ++i) {
            str->operator[](i) =  base[(*this)[i]];
        }
        return str;
    }

    static uint8_t Serial(char c) { return s_SerialTable[c]; }
    static char Complement(char c) { return s_ComplTable[c]; }
    static bool Check(const std::string &seq) {
        return std::find_if(seq.begin(), seq.end(), [](char c)->bool{
                return s_SerialTable[c] == (uint8_t) -1;
            }) == seq.end(); 
    }


//protected:
    static DnaSerialTable s_SerialTable;
    static DnaComplementTable s_ComplTable;
    size_t len_;
    std::vector<uint8_t> data_;
};

} // namespace fsa {
    
#endif // FSA_SEQUENCE_HPP

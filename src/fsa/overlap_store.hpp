#ifndef FSA_OVERLAP_STORE_HPP
#define FSA_OVERLAP_STORE_HPP


#include <array>
#include <vector>
#include <list>
#include <deque>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>

#include "overlap.hpp"
#include "logger.hpp"
#include "read_store.hpp"
#include "sequence.hpp"
#include "utility.hpp"
#include "file_io.hpp"


namespace fsa {
class ReadStore;
class OverlapStore {
public:
    OverlapStore() : OverlapStore(empty_read_store_) {}
    OverlapStore(ReadStore &rs) : read_store_(rs) {}

    template<typename C=bool (*)(Overlap &o)>
    void Load(const std::string &fname, const std::string &type="", size_t thread_size=1, C check=[](Overlap &o) {return true; });
    static std::string DetectFileType(const std::string &fname);

    template<typename S, typename C=bool (*)(const Overlap &o)>
    void Append(const std::string &fname, const std::string &type="", const S& s=S(), C check=[](const Overlap &o) {return true; });

    template<typename C=bool (*)(Overlap &o)>
    void LoadM4aFile(const std::string &fname, size_t thread_size=1, C check= [](Overlap &o) {return true; }) {
        if (thread_size > 1) LoadFileMt(fname, &OverlapStore::FromM4aLine, check, thread_size);
        else                 LoadFile(fname, &OverlapStore::FromM4aLine, check);
    }

    template<typename C = bool(*)(Overlap &o)>
    void LoadM4File(const std::string &fname, size_t thread_size=1,  C check = [](Overlap &o) {return true; }) {
        if (thread_size > 1) LoadFileMt(fname, &OverlapStore::FromM4Line, check, thread_size);
        else                 LoadFile(fname, &OverlapStore::FromM4Line, check);
    }
    
    template<typename C = bool(*)(Overlap &o)>
    void LoadOvlFile(const std::string &fname, size_t thread_size=1, C check = [](Overlap &o) {return true; }) {        
        if (thread_size > 1) LoadFileMt(fname, &OverlapStore::FromOvlLine, check, thread_size);
        else                 LoadFile(fname, &OverlapStore::FromOvlLine, check);
    }
    
    template<typename C = bool(*)(Overlap &o)>
    void LoadPafFile(const std::string &fname, size_t thread_size=1, C check = [](Overlap &o) {return true; }) {
        if (thread_size > 1) LoadFileMt(fname, &OverlapStore::FromPafLine, check, thread_size);
        else                 LoadFile(fname, &OverlapStore::FromPafLine, check);
    }
    
    template<typename C = bool(*)(const Overlap &o)>
    void Save(const std::string &fname, const std::string &type, C check=[](const Overlap &o){return true;}); 

    template<typename C =bool (*)(const Overlap &o)>
    void SaveM4File(const std::string &fname, C check=[](const Overlap &o){return true;}) const {
        SaveFile(fname, &OverlapStore::ToM4Line, check);
    }
    
    template<typename C =bool (*)(const Overlap &o)>
    void SaveM4aFile(const std::string &fname, C check=[](const Overlap &o){return true;}) const {
        SaveFile(fname, &OverlapStore::ToM4aLine, check);
    }

    template<typename C =bool (*)(const Overlap &o)>
    void SavePafFile(const std::string &fname, C check=[](const Overlap &o){return true;}) const {
        SaveFile(fname, &OverlapStore::ToPafLine, check);
    }


    template<typename S, typename C>
    void AppendM4File(const std::string &fname, const S& s, C check) {
        AppendFile(fname, s, &OverlapStore::ToM4Line, check);
    }

    template<typename S, typename C>
    void AppendM4aFile(const std::string &fname, const S& s, C check){
        AppendFile(fname, s, &OverlapStore::ToM4aLine, check);
    }


    template<typename S, typename C>
    void AppendPafFile(const std::string &fname, const S& s, C check){
        AppendFile(fname, s, &OverlapStore::ToPafLine, check);
    }



    size_t Size() const { return overlaps_.size(); }
    Overlap& Get(size_t i) { return  overlaps_[i]; }
    const std::deque<Overlap>& Get() const { return overlaps_; }
    std::deque<Overlap>& Get() { return overlaps_; }

    std::string GetReadName(int id) {
        return read_store_.IdToName(id);
    }

    Seq::Id GetReadId(const std::string &name) {
        return read_store_.GetIdByNameSafe(name);
    }

    std::array<Seq::Id, 2> GetReadIdRange() const;
    const ReadStore& GetReadStore() const { return read_store_; }

    std::unordered_map<int, std::unordered_map<int, Overlap*>> Group() const;
    std::unordered_map<int, std::unordered_map<int, Overlap*>> Group();
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> Group(bool (*better)(const Overlap&, const Overlap&)) const;
    std::unordered_map<int, std::unordered_map<int, Overlap*>> GroupTarget(bool (*better)(const Overlap* a, const Overlap *b));
    std::unordered_map<int, std::unordered_map<int, Overlap*>> GroupQuery();

    template<typename F, typename C>
    void LoadFile(const std::string &fname, F lineToOl, C check);

    template<typename F, typename C>
    void LoadFileMt(const std::string &fname, F lineToOl, C check, size_t thread_size=1);

    template<typename L, typename C>
    void SaveFile(const std::string &fname, L toLine, C check) const;

    template<typename S, typename L, typename C>
    void AppendFile(const std::string &fname, const S& s, L toLine, C check);

    bool FromM4Line(const std::string &line, Overlap &o);
    bool FromM4aLine(const std::string &line, Overlap &o);
    bool FromOvlLine(const std::string &line, Overlap &o);
    bool FromPafLine(const std::string &line, Overlap &o);
    
    std::string ToM4aLine(const Overlap& o) const;
    std::string ToM4Line(const Overlap& o) const;
    std::string ToPafLine(const Overlap &o) const;

protected:
    std::deque<Overlap> overlaps_;

    ReadStore &read_store_;
    ReadStore empty_read_store_;

};


template<typename S, typename C>
void OverlapStore::Append(const std::string &fname, const std::string &type, const S& s, C check) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "m4") {
        AppendM4File(fname, s, check);
    } else if (t == "m4a") {
        AppendM4aFile(fname, s, check);
    } else if (t == "paf") {
        AppendPafFile(fname, s, check);
    } else if (t == "ovl") {
        //SaveOvlFile(fname, check);
        LOG(FATAL)("TODO");
    } else {
        LOG(ERROR)("Failed to recognize overlap files type: %s", t.c_str());
    }
}


template<typename F, typename C>
void OverlapStore::LoadFile(const std::string &fname, F lineToOl, C check) {
    std::string line;
    //std::ifstream in(fname);
    GzFileReader in(fname);
    //if (in.is_open()) {
    if (in.Valid()) {
        line = in.GetStrippedLine();
        while (!line.empty()) {
        //while (std::getline(in, line)) {
            Overlap o;
            if ((this->*lineToOl)(line, o)) {
                if (check(o)) {
                    overlaps_.push_back(o);
                }
            } else {
                LOG(FATAL)("Failed to convert line to overlap \n   %s", line.c_str());
            }
            line = in.GetStrippedLine();
        }

        if (!in.IsEnd()) {
            LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
        }
    }
    else {
        LOG(FATAL)("Failed to load file: %s", fname.c_str());
    }
}

template<typename F, typename C>
void OverlapStore::LoadFileMt(const std::string &fname, F lineToOl, C check, size_t thread_size) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;
    GzFileReader in(fname);
    //std::ifstream in(fname);
    const size_t block_size = 500;

    auto generate_func = [&mutex_gen, &in](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_gen);

        size_t index = 0;
        /*
        for (; index<lines.size(); ++index) {
            //if (!std::getline(in, lines[index])) break;
            lines[index] = in.GetStrippedLine();
            if (lines[index].empty()) break;
        }
        */
        return in.GetLines(lines);
        
        return index;
    };

    auto combine_func = [&mutex_comb, this](const std::vector<Overlap> &ols, size_t sz) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        overlaps_.insert(overlaps_.end(), ols.begin(), ols.begin()+sz);
    };

    auto work_func = [&check, this, lineToOl, block_size, generate_func, combine_func, &fname, &in](size_t id) {
        std::vector<std::string> lines(block_size);
        std::vector<Overlap> ols(block_size);
        while (true) {
            size_t line_size = generate_func(lines);
            
            if (line_size > 0) {
                size_t ol_size = 0;
                for (size_t i=0; i<line_size; ++i) {
                    Overlap o;
                    if ((this->*lineToOl)(lines[i], o)) {
                        if (check(o)) {
                            ols[ol_size++] = o;
                        }
                    }
                    else {
                        LOG(FATAL)("Failed to convert line to overlap \n   %s", lines[i].c_str());
                    }
                }
                combine_func(ols, ol_size);;
            } else {
                break;
            }

        }   
        if (!in.IsEnd()) {
            LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
        }
    };

    //if (in.is_open()) {
    if (in.Valid()) {
        MultiThreadRun(thread_size, work_func);
    } else {
        LOG(FATAL)("Failed to load file: %s", fname.c_str());
    }
}

template<typename L, typename C>
void OverlapStore::SaveFile(const std::string &fname, L toLine, C check) const {
    std::ofstream of(fname);

    if (of.is_open()) {
        for (const auto &o : overlaps_) {
            if (check(o)) {
                of << (this->*toLine)(o);
            }
        }
    }
}


template<typename C>
void OverlapStore::Load(const std::string &fname, const std::string &type, size_t thread_size, C check) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "m4" || t == "m4.gz") {
        LoadM4File(fname, thread_size, check);
    } else if (t == "m4a" || t == "m4a.gz") {
        LoadM4aFile(fname, thread_size, check);
    } else if (t == "paf" || t == "paf.gz") {
        LoadPafFile(fname, thread_size, check);
    } else if (t == "ovl" || t == "ovl.gz") {
        LoadOvlFile(fname, thread_size, check);
    } else {
        LOG(ERROR)("Failed to recognize overlap files type: %s", t.c_str());
    }
}

template<typename C>
void OverlapStore::Save(const std::string &fname, const std::string &type, C check) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "m4") {
        SaveM4File(fname, check);
    } else if (t == "m4a") {
        SaveM4aFile(fname, check);
    } else if (t == "paf") {
        SavePafFile(fname, check);
    } else if (t == "ovl") {
        //SaveOvlFile(fname, check);
        LOG(FATAL)("TODO");
    } else {
        LOG(ERROR)("Failed to recognize overlap files type: %s", t.c_str());
    }
}



template<typename S, typename L, typename C>
void OverlapStore::AppendFile(const std::string &fname, const S& s, L toLine, C check) {
    std::ofstream of(fname, std::fstream::app);

    if (of.is_open()) {
        for (const auto &o : s) {
            if (check(o)) {
                of << (this->*toLine)(o);
            }
        }
    }
}

} // namespace fsa {

#endif // FSA_OVERLAP_STORE_HPP

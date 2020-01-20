#ifndef FSA_CONTIG_LINK_STORE_HPP
#define FSA_CONTIG_LINK_STORE_HPP

#include <array>
#include <cassert>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "overlap_store.hpp"
#include "sequence.hpp"
#include "read_store.hpp"
#include "overlap.hpp"
#include "contig_link.hpp"

namespace fsa {

class ContigLinkStore {
public:


    ContigLinkStore(ReadStore &rs):read_store_(rs), ctg2ctg_(rs), read2ctg_(rs) {}

    void SetParameter(const std::string& name, int v); 
    void SetParameter(const std::string& name, double v); 
    template<typename T> T GetParameter(const std::string& name);

    void LoadR2cFile(const std::string &fname);
    void LoadC2cFile(const std::string &fname);
    void AnalyzeSupport();
    ContigLink::Loc Location(const ContigLink& link) const;
    std::unordered_map<int, std::unordered_map<int, ContigLink>>& Get() { return links_; }
  
    void PurgeLinks();
    void IdentifyPaths();
    void Dump(const std::string &fname);
protected:

    double read2ctg_min_identity_{ 80 };
    double ctg2ctg_min_identity_{ 90 };

    int read2ctg_max_overhang_{ 100 };
    int ctg2ctg_max_overhang_{ 100 };

    int read_min_length_{ 10000 };
    int ctg_min_length_{ 50000 };

    int read2ctg_min_aligned_length_{ 5000 };
    int ctg2ctg_min_aligned_length_{ 5000 };

    int read2ctg_min_coverage_ { 3 };
    int window_size_ { 1000 };

    int thread_size_ {1};

    std::unordered_map<int, std::unordered_map<int, Overlap*>> read2ctg_group_;
    std::unordered_map<int, std::unordered_map<int, ContigLink>> links_;

    ReadStore &read_store_;
    OverlapStore ctg2ctg_;
    OverlapStore read2ctg_;
    
};


template<> int ContigLinkStore::GetParameter<int>(const std::string& name);
template<> double ContigLinkStore::GetParameter<double>(const std::string& name);

} // namespace fsa {

#endif // FSA_CONTIG_LINK_STORE_HPP  


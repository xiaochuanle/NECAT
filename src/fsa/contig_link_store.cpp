#include "contig_link_store.hpp"

#include <array>
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <iostream>
#include <tuple>

#include "file_io.hpp"
#include "overlap.hpp"
#include "utility.hpp"

namespace fsa {

void ContigLinkStore::SetParameter(const std::string& name, int v) {
    if (name == "read_min_length") read_min_length_ = v;
    else if (name == "ctg_min_length") ctg_min_length_ = v;
    else if (name == "read2ctg_max_overhang") read2ctg_max_overhang_ = v;
    else if (name == "ctg2ctg_max_overhang") ctg2ctg_max_overhang_ = v;
    else if (name == "read2ctg_min_aligned_length") read2ctg_min_aligned_length_ = v;
    else if (name == "ctg2ctg_min_aligned_length") ctg2ctg_min_aligned_length_ = v;
    else if (name == "read2ctg_min_coverage") read2ctg_min_coverage_ = v;
    else if (name == "thread_size") thread_size_ = v;
    else if (name == "window_size") window_size_ = v;
    else assert(!"never come here");
}


void ContigLinkStore::SetParameter(const std::string& name, double v) {
    if (name == "read2ctg_min_identity") read2ctg_min_identity_ = v;
    else if (name == "ctg2ctg_min_identity") ctg2ctg_min_identity_ = v;
    else assert(!"never come here");
}

template<>
int ContigLinkStore::GetParameter<int>(const std::string& name) {
    if (name == "read_min_length") return read_min_length_;
    else if (name == "ctg_min_length") return ctg_min_length_;
    else if (name == "read2ctg_max_overhang") return read2ctg_max_overhang_;
    else if (name == "ctg2ctg_max_overhang") return ctg2ctg_max_overhang_;
    else if (name == "read2ctg_min_aligned_length") return read2ctg_min_aligned_length_;
    else if (name == "ctg2ctg_min_aligned_length") return ctg2ctg_min_aligned_length_;
    else if (name == "read2ctg_min_coverage") return read2ctg_min_coverage_;
    else if (name == "thread_size") return thread_size_;
    else if (name == "window_size") return window_size_;
    else assert(!"never come here");

    return 0;
}

template<>
double ContigLinkStore::GetParameter<double>(const std::string& name) {
    if (name == "read2ctg_min_identity")    return read2ctg_min_identity_;
    else if (name == "ctg2ctg_min_identity") return ctg2ctg_min_identity_;
    else assert(!"never come here");

    return 0;
}

void ContigLinkStore::LoadR2cFile(const std::string &fname) {

    auto filter_simple = [&](const Overlap& o)->bool {
        return o.identity_ >= read2ctg_min_identity_ && o.a_.id != o.b_.id &&
               o.a_.len >= read_min_length_ && o.b_.len >= ctg_min_length_ && 
               o.AlignedLength() >= (size_t)std::min(read2ctg_min_aligned_length_/3, ctg_min_length_) && 
               o.Location(read2ctg_max_overhang_) != Overlap::Loc::Abnormal;
    };

    read2ctg_.Load(fname, "", thread_size_, filter_simple);

    auto better = [](const Overlap* a, const Overlap *b) { return a->AlignedLength() > b->AlignedLength(); };
    read2ctg_group_ = read2ctg_.GroupTarget(better);

    for (auto &ctg0 : read2ctg_group_) {
        for (auto &ctg1 : read2ctg_group_) {

            if (ctg0.first < ctg1.first) {

                auto keys = FindIntersectKeys(ctg0.second, ctg1.second);

                for (auto k : keys) {
                    auto &o0 = *ctg0.second[k];
                    auto &o1 = *ctg1.second[k];
                    
                    if (ContigLink::SimpleValid(o0, o1, read2ctg_max_overhang_)) {
                        links_[ctg0.first][ctg1.first].Add(o0, o1);
                    }
                }
            }
        }

    }
}


void ContigLinkStore::LoadC2cFile(const std::string &fname) {

    auto filter_simple = [&](const Overlap& o) -> bool {
        return o.identity_ >= ctg2ctg_min_identity_ && o.a_.id != o.b_.id &&
               o.a_.len >= ctg_min_length_ && o.b_.len >= ctg_min_length_ && 
               o.AlignedLength() >= (size_t)ctg2ctg_min_aligned_length_ && 
               o.Location(ctg2ctg_max_overhang_) != Overlap::Loc::Abnormal;
    };

    ctg2ctg_.Load(fname, "", thread_size_, filter_simple);

    for (auto &o : ctg2ctg_.Get()) {
        if (ContigLink::SimpleValid(o, ctg2ctg_max_overhang_)) {
            assert(o.a_.id != o.b_.id);
            if (o.a_.id < o.b_.id) {
                links_[o.a_.id][o.b_.id].Add(o, o.a_, o.b_, ctg2ctg_max_overhang_);
            }
            else {
                links_[o.b_.id][o.a_.id].Add(o, o.b_, o.a_, ctg2ctg_max_overhang_);
            }
        }
    }
}

void ContigLinkStore::AnalyzeSupport() {
    for (auto &i0 : links_) {
        for (auto &i1 : i0.second) {
            assert (i0.first < i1.first);
            i1.second.AnalyzeLinks(read2ctg_max_overhang_, ctg2ctg_max_overhang_, read2ctg_min_coverage_, read2ctg_min_aligned_length_, window_size_);
        }
    }
}

ContigLink::Loc ContigLinkStore::Location(const ContigLink& bunch) const {
    if (bunch.best_c2c != nullptr) {
        return bunch.best_c2c->Location(ctg2ctg_max_overhang_);
    } else if (bunch.BestC2r2c() != nullptr) {
        return bunch.BestC2r2c()->Location(read2ctg_max_overhang_);
    } else {
        return ContigLink::Loc::Abnormal;
    }
}

void ContigLinkStore::PurgeLinks() {
    // Collect all links
    std::vector<ContigLink*> links;
    for (auto &ctg0 : links_) {
        for (auto &ctg1 : ctg0.second) {
            links.push_back(&ctg1.second);
        }
    }

    std::sort(links.begin(), links.end(), [](ContigLink* a, ContigLink *b) {
        return a->Score() > b->Score();
    });

    
    std::unordered_set<Seq::EndId> done;
    for (auto link : links) {
        Seq::EndId sB = Seq::IdToEndId(link->Best()->source.id, 0);
        Seq::EndId sE = Seq::IdToEndId(link->Best()->source.id, 1);
        Seq::EndId tB = Seq::IdToEndId(link->Best()->target.id, 0);
        Seq::EndId tE = Seq::IdToEndId(link->Best()->target.id, 1);

        ContigLink::Loc loc = Location(*link);
        if (link->Best()->SameStrand()) {
            if (loc ==  ContigLink::Loc::Left) {
                //  s -------->
                //  t              ----------->
                if (done.find(sE) == done.end() && done.find(tB) == done.end()) {
                    done.insert(sE);
                    done.insert(tB);
                } else {
                    link->Removed(true);
                }
            } else if (loc ==  ContigLink::Loc::Right) {
                // s               --------->
                // t -------->
                if (done.find(sB) == done.end() && done.find(tE) == done.end()) {
                    done.insert(sB);
                    done.insert(tE);
                } else {
                    link->Removed(true);
                }
            } else {
                // TODO
                link->Removed(true);
            }
        } else {
            if (loc ==  ContigLink::Loc::Left) {
                //  s <--------
                //  t              ----------->
                if (done.find(sB) == done.end() && done.find(tB) == done.end()) {
                    done.insert(sB);
                    done.insert(tB);
                } else {
                    link->Removed(true);
                }

            } else if (loc == ContigLink::Loc::Right){
                assert(loc == ContigLink::Loc::Right);
                // s               <---------
                // t -------->
                if (done.find(sE) == done.end() && done.find(tE) == done.end()) {
                    done.insert(sE);
                    done.insert(tE);
                } else {
                    link->Removed(true);
                }
            } else {
                    link->Removed(true);
            }
        }
    }
}

void ContigLinkStore::IdentifyPaths() {
    std::unordered_set<ContigLink*> done;
    for (auto &ctg0 : links_) {
        for (auto &ctg1 : ctg0.second) {
            auto link = &ctg1.second;
            if (!link->Removed() && done.find(link) == done.end()) {
                std::list<ContigLink*> path;
                // Left
                printf("link %d, %d\n", link->Best()->source.id, link->Best()->target.id);

            }
        }
    }
}
void ContigLinkStore::Dump(const std::string &fname) {
    gzFile of = gzopen(fname.c_str(), "w");
    if (of != nullptr) {
        for (auto &ctg0 : links_) {
            for (auto &ctg1 : ctg0.second) {
                auto& link = ctg1.second;
                gzprintf(of, "%d %s %d %s\n", ctg0.first, read_store_.IdToName(ctg0.first).c_str(), ctg1.first, read_store_.IdToName(ctg1.first).c_str());
                //of << ctg0.first << " " << read_store_.IdToName(ctg0.first) << " " << ctg1.first << " " << read_store_.IdToName(ctg1.first) << "\n";
                gzprintf(of, "Contig_Contig\n");
                //of << "Contig_Contig\n";
                for (auto &i : link.c2c_links) {
                    for (auto& o : i.ols) {
                        gzprintf(of, "%s", ctg2ctg_.ToM4aLine(*o).c_str());
                        //of << ctg2ctg_.ToM4aLine(*o);
                    }
                    gzprintf(of, "ol_expect:%d %d\n", i.ol_expect[0], i.ol_expect[1]);
                    //of << "ol_expect:" << i.ol_expect[0] << " " << i.ol_expect[1] << "\n";
                }
                gzprintf(of, "Read_Contig\n");
                //of << "Read_Contig\n";
                for (auto &i : link.c2r2c_links) {
                    for (auto& o : i.ols) {
                        gzprintf(of, "%s", read2ctg_.ToM4aLine(*o).c_str());
                        //of << read2ctg_.ToM4aLine(*o);
                    }
                    gzprintf(of, "s2t:%d %d %d %d\n", i.pos_s2t[0], i.pos_s2t[1], i.pos_s2t[2], i.pos_s2t[3]);
                    //of << "s2t:" << i.pos_s2t[0] << " " << i.pos_s2t[1] << " " << i.pos_s2t[2] << " " << i.pos_s2t[3] << "\n";
                    gzprintf(of, "ol_expect:%d %d\n", i.ol_expect[0], i.ol_expect[1]);
                    //of << "ol_expect:" << i.ol_expect[0] << " " << i.ol_expect[1] << "\n";
                }

                if (link.Best() != nullptr) {
                    gzprintf(of, "Best\n");
                    //of << "Best\n";
                    for (auto o : link.Best()->ols) {
                        gzprintf(of, "%s\n", read2ctg_.ToM4aLine(*o).c_str());
                        //of << read2ctg_.ToM4aLine(*o);

                    }
                }
            }
        }
        gzclose(of);
    }
    else {
        LOG(ERROR)("Failed to write file %s", fname.c_str());
    }

}

} // namespace fsa {
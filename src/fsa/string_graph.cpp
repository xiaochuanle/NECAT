#include "string_graph.hpp"

#include <algorithm>
#include <cstdio>
#include <cassert>
#include <list>
#include <unordered_set>
#include <iterator> 

#include "overlap.hpp"
#include "logger.hpp"
#include "file_io.hpp"

namespace fsa {

std::vector<StringNode*> StringNode::GetOutNodes() {
    std::vector<StringNode*> nodes;
    for (auto e : out_edges_) {
        nodes.push_back(e->out_node_);
    }
    return nodes;
}

std::vector<StringNode*> StringNode::GetInNodes() {
    std::vector<StringNode*> nodes;
    for (auto e : in_edges_) {
        nodes.push_back(e->in_node_);
    }
    return nodes;
}


std::vector<StringNode*> StringNode::GetAllOutNodes() {
    std::vector<StringNode*> nodes;
    for (auto e : out_edges_) {
        nodes.push_back(e->out_node_);
    }
    for (auto e : reduced_out_edges_) {
        nodes.push_back(e->out_node_);
    }
    return std::move(nodes);
}
std::vector<StringNode*> StringNode::GetAllInNodes() {
    std::vector<StringNode*> nodes;
    for (auto e : in_edges_) {
        nodes.push_back(e->in_node_);
    }
    for (auto e : reduced_in_edges_) {
        nodes.push_back(e->in_node_);
    }
    return std::move(nodes);

}

void StringNode::ReduceEdge(std::vector<StringEdge*> &src, std::vector<StringEdge*> &dst, const StringEdge *e) {
    auto it = std::find(src.begin(), src.end(), e);
    
    assert(it != src.end());
    dst.push_back(*it);
    src.erase(it);
}



std::tuple<int, bool, int, int> StringEdge::GetSeqArea() const {
    return std::make_tuple(out_node_->ReadId(), out_node_->id_ > 0, start_, end_);
}

StringGraph::~StringGraph() {
    for (auto& i : nodes_) {
        delete i.second;
    }
    for (auto& i : edges_) {
        delete i.second;
    }
}

std::string StringGraph::NodeIdString(int id) {
    char buff[24];
    if (id > 0) {
        sprintf(buff, "%09d:B", id-1);
    }
    else if(id < 0){
        sprintf(buff, "%09d:E", -id-1);
    }
    else {
        sprintf(buff, "NA");
    }

    return buff;   
}

std::string StringGraph::ReadIdString(int id) {
    char buff[24];
    sprintf(buff, "%09d", id);
    return buff;
}


void StringGraph::AddOverlap(const Overlap* overlap) {

	Seq::EndId fB = Seq::IdToEndId(overlap->a_.id, 0);
	Seq::EndId fE = Seq::IdToEndId(overlap->a_.id, 1);
	Seq::EndId gB = Seq::IdToEndId(overlap->b_.id, 0);
	Seq::EndId gE = Seq::IdToEndId(overlap->b_.id, 1);
    
	if (overlap->a_.start > 0) {
        assert(overlap->a_.end == overlap->a_.len);
		if (overlap->SameDirect()) {
			//  f.B         f.E
			//	f----------->
			//	g       ------------->
			//	      g.B           g.E
            assert(overlap->b_.start == 0 && overlap->b_.end < overlap->b_.len);
           
			AddEdge(fE, gE, overlap->b_.len-overlap->b_.end, -overlap->score_, overlap->identity_, overlap->b_.id, overlap->b_.end, overlap->b_.len);
			AddEdge(gB, fB, overlap->a_.start, -overlap->score_, overlap->identity_, overlap->a_.id, overlap->a_.start, 0);
		}
		else {
			// f.B         f.E
			//	f----------->
			//	g         <------------ -
			//	          g.E           g.B
            assert(overlap->b_.start > 0 && overlap->b_.end == overlap->b_.len);

			AddEdge(fE, gB, overlap->b_.start-0, -overlap->score_, overlap->identity_, overlap->b_.id, overlap->b_.start, 0);
			AddEdge(gE, fB, overlap->a_.start-0, -overlap->score_, overlap->identity_, overlap->a_.id, overlap->a_.start, 0);
			
		}
	}
	else {
        assert(overlap->a_.end < overlap->a_.len);
		if (overlap->SameDirect()) {
			
			//       f.B         f.E
			//  f     ----------->
			//	g------------->
			//	g.B           g.E
            assert(overlap->b_.start > 0 && overlap->b_.end == overlap->b_.len);
			AddEdge(fB, gB, overlap->b_.start, -overlap->score_, overlap->identity_, overlap->b_.id, overlap->b_.start, 0);
			AddEdge(gE, fE, overlap->a_.len-overlap->a_.end, -overlap->score_, overlap->identity_, overlap->a_.id, overlap->a_.end, overlap->a_.len);
		}
		else {
			//        f.B         f.E
			// f       ----------->
			//	g <------------ 
			//	g.E           g.B
            assert(overlap->b_.start == 0 && overlap->b_.end < overlap->b_.len);
			AddEdge(fB, gE, overlap->b_.len-overlap->b_.end, -overlap->score_, overlap->identity_, overlap->b_.id, overlap->b_.end, overlap->b_.len);
			AddEdge(gB, fE, overlap->a_.len-overlap->a_.end, -overlap->score_, overlap->identity_, overlap->a_.id, overlap->a_.end, overlap->a_.len);
		}
	}

}

void StringGraph::AddOverlaps(const std::deque<Overlap> &overlaps, int min_length, int min_aligned_lenght, float min_identity) {

    std::unordered_set<StringNode::ID> contained;

    for (auto iter = overlaps.begin(); iter != overlaps.end(); ++iter) {
        if (iter->IsContaining(0)) {
            contained.insert(iter->b_.id);
        } else if (iter->IsContained(0)) {
            contained.insert(iter->a_.id);
        }
    }

    LOG(INFO)("contained = %zd", contained.size());
	std::unordered_set<StringEdge::ID, StringEdge::Hash, StringEdge::Compare> done;
    for (auto &o : overlaps) {
        auto ids = std::minmax(o.a_.id, o.b_.id);
        StringEdge::ID id_pair = {ids.first, ids.second};

        if (done.find(id_pair) == done.end()) {

            if (!FilterOverlap(o, contained, min_length, min_aligned_lenght, min_identity)) {
                AddOverlap(&(o));
                done.insert(id_pair);

            }
        }
    }
    LOG(INFO)("Done = %zd", done.size());
}

bool StringGraph::FilterOverlap(const Overlap &ovlp, const std::unordered_set<StringNode::ID> &contained, int min_length, int min_aligned_length, float min_identity) {
    if (contained.find(ovlp.a_.id) != contained.end() || contained.find(ovlp.b_.id) != contained.end())
        return true;

    if (ovlp.a_.id == ovlp.b_.id) 
        return true;

    if (!ovlp.IsProper(0))
        return true;
    if (ovlp.identity_ < min_identity)
        return true;
    if (ovlp.a_.len < min_length || ovlp.b_.len < min_length)
        return true;

    if (ovlp.AlignedLength() < (size_t)min_aligned_length)
        return true;
    return false;
}

void StringGraph::AddEdge(int in_node, int out_node, int len, int score, double identity, int read, int start, int end) {
	auto in = nodes_.find(in_node);
	if (in == nodes_.end()) {
		auto r = nodes_.insert(std::make_pair(in_node, new StringNode(in_node)));
		assert(r.second);
		in = r.first;
	}

	auto out = nodes_.find(out_node);
	if (out == nodes_.end()) {
		auto r = nodes_.insert(std::make_pair(out_node, new StringNode(out_node)));
		assert(r.second);
		out = r.first;
	}
	
	StringEdge *e = new StringEdge(in->second, out->second);
	edges_[StringEdge::ID{in_node, out_node}] = e;
    e->length_ = len;
    e->score_ = score;
    e->identity_ = identity;
    e->read_ = read;
    e->start_ = start;
    e->end_ = end;

    in->second->out_edges_.push_back(e);
    out->second->in_edges_.push_back(e);
}

void StringGraph::MarkTransitiveEdges() {

	const int FUZZ = 500;

	for (auto it = nodes_.begin(); it != nodes_.end(); ++it) {
		it->second->mark_ = 'V';
	}

	for (auto it : nodes_) {
		if (it.second->out_edges_.size() == 0) continue;


		std::vector<StringEdge*> &out_edges = it.second->out_edges_;



		std::sort(out_edges.begin(), out_edges.end(), [](StringEdge* a, StringEdge *b) { return a->length_ < b->length_; });

		for (auto e : out_edges) {
			e->out_node_->mark_ = 'I';
		}

		int max_len = out_edges.back()->length_  + FUZZ;

		for (auto e : out_edges) {
			StringNode* w = e->out_node_;
			if (w->mark_ == 'I') {
				std::sort(w->out_edges_.begin(), w->out_edges_.end(), [](StringEdge* a, StringEdge *b) { return a->length_ < b->length_; });
				for (auto e2 : w->out_edges_) {
					if (e2->length_ + e->length_ < max_len) {
						if (e2->out_node_->mark_ == 'I') {
							e2->out_node_->mark_ = 'E';
						}
					}
				}
			}
		}
		for (auto e : out_edges) {
			StringNode* w = e->out_node_;
			std::sort(w->out_edges_.begin(), w->out_edges_.end(), [](StringEdge* a, StringEdge *b) { return a->length_ < b->length_; });
			if (w->out_edges_.size() > 0) {
				if (w->out_edges_.front()->out_node_->mark_ == 'I') {
					w->out_edges_.front()->out_node_->mark_ = 'E';
                }
			}				
			for (auto e2 : w->out_edges_) {
				if (e2->length_ < FUZZ) {
					if (e2->out_node_->mark_ == 'I') {
						e2->out_node_->mark_ = 'E';
					}
				}
			}
		}

		for (auto e : out_edges) {
			if (e->out_node_->mark_ == 'E') {
                e->reduce_ = true;
                StringGraph::ReverseEdge(e)->reduce_ = true;
			}
			e->out_node_->mark_ = 'V';
		}
	}

    for (auto& i : edges_) {
        if (i.second->reduce_) i.second->Reduce(StringEdge::TRANSITIVE);
    }

}

/** 
 * 移除毛刺节点：有多个入（出）节点，它入（出）节点没有入（出）节点则认为是毛刺�?
 */
void StringGraph::MarkSpurEdges() {

    std::unordered_set<StringEdge*> removed;
    for (auto &n :  nodes_) {
        
        if (n.second->out_edges_.size() > 1) {
            size_t count = 0;           
            for (auto e : n.second->out_edges_) {
                assert(!e->reduce_);
                if (!e->reduce_) { // TODO 
                    if (e->out_node_->out_edges_.size()+ e->out_node_->reduced_out_edges_.size() == 0) {
                        removed.insert(e);
                        removed.insert(ReverseEdge(e));
                        count++;
                    }
                }
                if (count +1 == n.second->out_edges_.size()) break;
            }
        }
/*
        if (n.second->in_edges_.size() > 1) {
            for (auto e : n.second->in_edges_) {

                assert(!e->reduce_);
                if (!e->reduce_) {  // TODO assert(!e->reduce_);
                    
                    if (e->in_node_->in_edges_.size() + e->in_node_->reduced_in_edges_.size() == 0) {
                        removed.insert(e);
                        removed.insert(ReverseEdge(e));
                    }
                }
             }
        }
*/
    }
    for (auto e : removed) {
        e->Reduce(StringEdge::SPUR);
    }
}

template <typename T, typename U>
bool HasCommon(const T &a, const U &b) 
{
    for (const auto &i : a) {
        if (std::find(b.begin(), b.end(), i) != b.end()) return true;
    }
    return false;
}

template<typename T>
void Intersection(const std::unordered_set<T> &a, std::unordered_set<T> &b, std::unordered_set<T> &c) {

    for (const auto &i : a) {
        if (b.find(i) != b.end())
            c.insert(i);
    }
}

/** 
 * 有多个入边和出边的节点，如果它的入节点的延伸和出节点的延伸没�?
 * 交集，则认为该节点为Chimer节点，标记与该节点相连的入边和出�?
 * 延伸：指一定深度的出节�?
 * TODO：是否需要考虑所有的边，reduce的边是否需要考虑。（初步判断应该考虑�?
 */
void StringGraph::MarkChimerEdges() {
    
    std::unordered_set<StringEdge*> removed;
    std::unordered_set<StringNode*> multi_in_nodes;
    std::unordered_set<StringNode*> multi_out_nodes;

    for (auto &n : nodes_) {
        std::vector<StringNode*>&& out_nodes = n.second->GetOutNodes();
        if (out_nodes.size() >= 2) {
            multi_out_nodes.insert(out_nodes.begin(), out_nodes.end());
        }

        std::vector<StringNode*>&& in_nodes = n.second->GetInNodes();
        if (in_nodes.size() >= 2) {
            multi_in_nodes.insert(in_nodes.begin(), in_nodes.end());
        }
    }
    std::unordered_set<StringNode*> chimer_candidates;

    Intersection(multi_out_nodes, multi_in_nodes, chimer_candidates);

    for (auto n : chimer_candidates) {
        std::vector<StringNode*>&& out_nodes = n->GetAllOutNodes();
        std::vector<StringNode*>&& in_nodes = n->GetAllInNodes();

        std::unordered_set<StringNode*> test;
        for (auto in_node : in_nodes) {
            std::vector<StringNode*>&& nodes = in_node->GetAllOutNodes();
            test.insert(nodes.begin(), nodes.end());
        }

        test.erase(n);
        
        if (!HasCommon(test, out_nodes)) {
            std::unordered_set<StringNode*> flow_node1;
            for (auto v : out_nodes) {
                std::unordered_set<StringNode*>&& nodes = BfsNodes(v, n, 5);
                flow_node1.insert(nodes.begin(), nodes.end());
            }
            std::unordered_set<StringNode*> flow_node2;
            for (auto v : test) {
                std::unordered_set<StringNode*>&& nodes = BfsNodes(v, n, 5);
                flow_node2.insert(nodes.begin(), nodes.end());
            }
            
            if (!HasCommon(flow_node1, flow_node2)) {
                for (auto e : n->out_edges_) {
                    assert(!e->reduce_);
                    removed.insert(e);
                    removed.insert(ReverseEdge(e));
                }
                for (auto e : n->in_edges_) {
                    assert(!e->reduce_);
                    removed.insert(e);
                    removed.insert(ReverseEdge(e));
                }
            }
            
        }
    }

    for (auto e : removed) {
        e->Reduce(StringEdge::CHIMER);
    }
    
}

std::unordered_set<StringNode*> StringGraph::BfsNodes(StringNode* n, StringNode *exclude, int depth) {
	std::unordered_set<StringNode*> result;
	result.insert(n);

	std::list<StringNode*> cand;
	cand.push_back(n);
    StringNode* depth_node = n;
	int dp = 1;
	while (dp < depth && cand.size() > 0) {
		StringNode* v = cand.front();
		cand.pop_front();

		for (auto e : v->out_edges_) {
			if (e->out_node_ != exclude) {
                if (result.find(e->out_node_) == result.end()) {
                    result.insert(e->out_node_);
                    if (e->out_node_->out_edges_.size() > 0) {
                        cand.push_back(e->out_node_);
                    }
                }
			}
		}

        for (auto e : v->reduced_out_edges_) {
            if (e->out_node_ != exclude) {
                if (result.find(e->out_node_) == result.end()) {
                    result.insert(e->out_node_);
                    if (e->out_node_->out_edges_.size() > 0) {
                        cand.push_back(e->out_node_);
                    }
                }
            }
        }
        if (v == depth_node) {

		    dp++;
            depth_node = cand.size() > 0 ? cand.back() : nullptr;
        }

	}
	return result;
}


void StringGraph::MarkBestOverlap() {
    std::unordered_set<StringEdge*> best_edges;

    auto best_cmp_func = [](StringEdge *a, StringEdge *b) { return a->score_ < b->score_; };
    //auto best_cmp_func = [](StringEdge *a, StringEdge *b) { return a->score_ * a->identity_ > b->score_ * b->identity_; };
     
    for (auto &n : nodes_) {
        if (n.second->out_edges_.size() > 0) {
            const std::vector<StringEdge*>& edges = n.second->out_edges_;
            auto m = std::max_element(edges.begin(), edges.end(), best_cmp_func);
            assert(!(*m)->reduce_);
            best_edges.insert(*m);
            n.second->best_out_ = *m;
        }
        if (n.second->in_edges_.size() > 0) {
            const std::vector<StringEdge*>& edges = n.second->in_edges_;
            auto m = std::max_element(edges.begin(), edges.end(), best_cmp_func);
            assert(!(*m)->reduce_);
            best_edges.insert(*m);
            n.second->best_in_ = *m;
        }
    }

    for (auto &e : edges_) { // TODO check if the condition can be removed
        if (!e.second->reduce_) {
            if (best_edges.find(e.second) == best_edges.end()) {
                e.second->Reduce(StringEdge::NO_BEST);
                ReverseEdge(e.second)->Reduce(StringEdge::NO_BEST);
            }
        }
    }
}

void StringGraph::ResolveRepeatEdges() {
    std::unordered_set<StringEdge*> edges_to_reduce;
    std::unordered_set<StringNode*> nodes_to_test;
    for (auto &i : nodes_) {
        auto n = i.second;
        if (n->in_edges_.size() == 1 && n->out_edges_.size() == 1) {
            nodes_to_test.insert(n);
        }
    }

    for (auto n : nodes_to_test) {
        auto in_node = n->in_edges_[0]->in_node_;
        auto out_node = n->out_edges_[0]->out_node_;

        for (auto e : in_node->out_edges_) {
            //auto vv = e->in_node_;
            auto ww = e->out_node_;

            auto ww_out_nodes = ww->GetAllOutNodes();
            auto v_out_nodes = n->GetAllOutNodes();
            bool overlap = HasCommon(ww_out_nodes, v_out_nodes);

            int ww_in_count = ww->in_edges_.size();

            if (ww != n && !e->reduce_ && ww_in_count > 1 && !overlap) {
                edges_to_reduce.insert(e);
            }

        }

        for (auto e : out_node->in_edges_) {
            auto vv = e->in_node_;
            //auto ww = e->out_node_;

            auto vv_in_nodes = vv->GetAllInNodes();
            auto v_in_nodes = n->GetAllInNodes();
            bool overlap = HasCommon(vv_in_nodes, v_in_nodes);

            int vv_out_count = vv->out_edges_.size();

            if (vv != n && !e->reduce_ && vv_out_count > 1 && !overlap) {
                edges_to_reduce.insert(e);
            }

        }
    }
    for (auto e : edges_to_reduce) {
        e->Reduce(StringEdge::REMOVED);
    }
}

void StringGraph::IdentifySimplePaths() {
    std::unordered_set<StringEdge*> visited;

    for (auto &i : edges_) {
        StringEdge* e = i.second; // short name
        if (!e->reduce_ && visited.find(e) == visited.end()) {
            paths_.push_back(ExtendSimplePath(e, visited));
            auto vpath = Reverse(paths_.back());
            for (auto e : vpath) visited.insert(e);
            paths_.push_back(vpath);
        }
    }

//    assert(Assert_PathDual(paths_));  TOO SLOW
}

std::list<StringEdge*> StringGraph::ExtendSimplePath(StringEdge *e, std::unordered_set<StringEdge*> &visited) {
    assert(visited.find(e) == visited.end());

    std::list<StringEdge*> path;

    std::unordered_set<StringNode*> rnodes;

    visited.insert(e);
    path.push_back(e);
    rnodes.insert(ReverseNode(e->in_node_));
    rnodes.insert(ReverseNode(e->out_node_));

    StringEdge* curr = path.back();
    while (curr->out_node_->in_edges_.size() == 1 && curr->out_node_->out_edges_.size() == 1 && 
        visited.find(curr->out_node_->out_edges_.front()) == visited.end() &&
        rnodes.find(curr->out_node_->out_edges_.front()->out_node_) == rnodes.end()) {

        path.push_back(curr->out_node_->out_edges_.front());
        visited.insert(curr->out_node_->out_edges_.front());
        rnodes.insert(ReverseNode(curr->out_node_->out_edges_.front()->out_node_));
        curr = curr->out_node_->out_edges_.front();
    }

    curr = path.front();
    while (curr->in_node_->in_edges_.size() == 1 && curr->in_node_->out_edges_.size() == 1 &&
        visited.find(curr->in_node_->in_edges_.front()) == visited.end() &&
        rnodes.find(curr->in_node_->in_edges_.front()->in_node_) == rnodes.end()) {

        path.push_front(curr->in_node_->in_edges_.front());
        visited.insert(curr->in_node_->in_edges_.front());
        rnodes.insert(ReverseNode(curr->in_node_->in_edges_.front()->in_node_));

        curr = curr->in_node_->in_edges_.front();
    }

    return path;
}


template<typename T, typename U>
bool BelongTo(const T& src, const U& dst) {
    for (auto i : src) {
        if (dst.find(i) == dst.end()) return false;
    }
    return true;
}


std::vector<StringEdge*> StringGraph::ShortestPath(const StringNode* src, const StringNode *dst, std::unordered_set<StringEdge*> doable, int(*score)(StringEdge*)) {

    typedef std::tuple<const StringNode*, StringEdge*, int> WorkType;
    std::vector<WorkType> nodes;
    std::unordered_map<const StringNode*, WorkType> done;
    WorkType srcdst{ dst, nullptr, 0 };        // ���������յ���ͬ�����?

    auto cmp = [](const WorkType &a, const WorkType &b) { return std::get<2>(a) < std::get<2>(b); };

    std::make_heap(nodes.begin(), nodes.end(), cmp);

    nodes.push_back(std::make_tuple(src, (StringEdge*)nullptr, 0));
    std::push_heap(nodes.begin(), nodes.end());

    while (nodes.size() > 0) {
        WorkType i = nodes.front();
        std::pop_heap(nodes.begin(), nodes.end());
        nodes.pop_back();

        // ���������յ���ͬ�����?
        if (std::get<0>(i) == src && src == dst) {
            if ((std::get<1>(srcdst) == nullptr && std::get<1>(i) != nullptr) ||
                (std::get<1>(srcdst) != nullptr && std::get<1>(i) != nullptr && std::get<2>(i) < std::get<2>(srcdst)))
                
                srcdst = i;
        }

        if (done.find(std::get<0>(i)) == done.end()) {
            done[std::get<0>(i)] = i;
            
            for (auto e : std::get<0>(i)->out_edges_) {
                if (doable.find(e) != doable.end()) {
                    nodes.push_back(std::make_tuple(e->out_node_, e, std::get<2>(i) + score(e)));
                    std::push_heap(nodes.begin(), nodes.end());
                }
            }
        }
    }

    std::vector<StringEdge*> path;

    auto r = done.end();
    if (src == dst && std::get<1>(srcdst) != nullptr) {
        path.push_back(std::get<1>(srcdst));
        r = done.find(std::get<1>(srcdst)->in_node_);
    }
    else {
        r = done.find(dst);
    }

    if (r != done.end() && r->first != src) {
        while (r->first != src) {
            path.push_back(std::get<1>(r->second));
            r = done.find(std::get<1>(r->second)->in_node_);
        }
    }

    std::reverse(path.begin(), path.end());

    return path;
   
}


std::list<StringEdge*> StringGraph::Reverse(const std::list<StringEdge*>& path) {
    std::list<StringEdge*> vpath;
    for (auto &e : path) {
        vpath.push_front(ReverseEdge(e));
    }
    return vpath;
    
}

bool StringGraph::Assert_PathDual(const std::list<std::list<StringEdge*>> paths) {
    auto reverse_path = [&](const std::list<StringEdge*> &p) {
        std::list<StringEdge*> vp;
        for (auto &e : p) {
            vp.push_front(ReverseEdge(e));
        }
        return vp;
    };

    auto is_equal = [](const std::list<StringEdge*> &a, const std::list<StringEdge*> &b) {
        if (a.size() != b.size()) return false;

        auto ia = a.begin();
        auto ib = b.begin();
        for (; ia != a.end(); ++ia, ++ib) {
            if (*ia != *ib) return false;
        }
        return true;
    };
    
    std::unordered_set<const std::list<StringEdge*>*> done;

    for (auto &p1 : paths) {
        if (done.find(&p1) == done.end()) {
            done.insert(&p1);

            auto vp1 = reverse_path(p1);
            bool found = false;
            for (auto &p2 : paths) {
                if (done.find(&p2) == done.end()) {
                    if (is_equal(p2, vp1)) {
                        done.insert(&p2);
                        found = true;
                        break;
                    }
                }
            }
            if (!found) return false;
        }
    }
    return true;
}

void StringGraph::Dump() {

	for (const auto &n : edges_) {
        auto e = n.second;
        if (!e->reduce_)
		    printf("%d, %d, %d, %d, %d\n", e->in_node_->id_, e->out_node_->id_, e->reduce_, e->length_, e->score_);
	}

}

void StringGraph::DumpPaths() {
    for (const auto &path : paths_) {
        for (auto p : path) {
            printf("%d, ", p->in_node_->id_);
        }
        printf("%d\n", path.back()->out_node_->id_);
        
    }
}

void StringGraph::SaveChimerNode(const std::string &fname) {
    FILE* file = fopen(fname.c_str(), "w");
    if (file != NULL) {
        for (auto i : nodes_) {
            auto n = i.second;
            if (n->in_edges_.size() == 0 and n->out_edges_.size() == 0) {
                bool is_chimer_node = false;
                for (auto e : n->reduced_in_edges_) {
                    if (e->type_ == StringEdge::CHIMER) {
                        is_chimer_node = true;
                        break;
                    }
                }
                if (!is_chimer_node) {
                    for (auto e : n->reduced_in_edges_) {
                        if (e->type_ == StringEdge::CHIMER) {
                            is_chimer_node = true;
                            break;
                        }
                    }
                }

                if (is_chimer_node) {
                    std::string&& str = NodeIdString(n->id_);
                    fprintf(file, "%s\n", str.c_str());
                }
            }

        }
        fclose(file);
    }
}
void StringGraph::SaveEdges(const std::string &fname) {
    gzFile file = gzopen(fname.c_str(), "w");
    if (file != nullptr) {
        for (auto i : edges_) {
            StringEdge* e = i.second;
            const char* type = e->type_ == StringEdge::CHIMER ? "C" :
                e->type_ == StringEdge::REMOVED ? "R" :
                e->type_ == StringEdge::SPUR ? "S" :
                e->type_ == StringEdge::TRANSITIVE ? "TR" :
                e->type_ == StringEdge::NO_BEST ? "NO_BEST" :
                "G";
            

            gzprintf(file, "%s %s %s %5d %5d %5d %5.2f %s\n", \
                NodeIdString(e->in_node_->id_).c_str(),
                NodeIdString(e->out_node_->id_).c_str(),
                ReadIdString(e->read_).c_str(),
                e->start_, e->end_,
                e->score_,
                e->identity_, type);

        }
        gzclose(file);
    }
}

} // namespace fsa {
    

#include "path_graph.hpp"

#include <algorithm>
#include <cstdio>
#include <cassert>
#include <list>
#include <unordered_set>
#include <iterator> 
#include <stack>
#include "logger.hpp"
#include "file_io.hpp"

namespace fsa {

std::vector<PathNode*> PathNode::GetInNodes() {
    std::vector<PathNode*> nodes;
    for (auto e : in_edges_) {
        nodes.push_back(e->in_node_);
    }
    return nodes;
}


std::vector<PathNode*> PathNode::GetOutNodes() {
    std::vector<PathNode*> nodes;
    for (auto e : out_edges_) {
        nodes.push_back(e->out_node_);
    }
    return nodes;
}

PathEdge* PathNode::GetBestInEdge() const {
    for (auto e : in_edges_) {
        e->Contain(string_node_->GetBestInEdge());
        return e;
    }
    return nullptr;
}

PathEdge* PathNode::GetBestOutEdge() const {
    for (auto e : out_edges_) {
        e->Contain(string_node_->GetBestOutEdge());
        return e;
    }
    return nullptr;
}

std::string PathEdge::IdToString(ID id, const char *sep) {
    std::string str = StringGraph::NodeIdString(std::get<0>(id));
    str += sep;
    str += StringGraph::NodeIdString(std::get<1>(id));;
    str += sep;
    str += StringGraph::NodeIdString(std::get<3>(id));
    return str;
}


std::string SimplePathEdge::ToString() const {
    assert(path_.size() > 0);
    std::string path_or_edges;
    path_or_edges = StringGraph::NodeIdString(path_.front()->in_node_->Id());
    path_or_edges += "~";
    path_or_edges += StringGraph::NodeIdString(path_.front()->out_node_->Id());
    path_or_edges += "~";
    path_or_edges += StringGraph::NodeIdString(path_.back()->out_node_->Id());
    return path_or_edges;

}

std::string SimplePathEdge::ToDetailString() const {
    assert(path_.size() > 0);

    std::string detail = StringGraph::NodeIdString(path_.front()->in_node_->Id());

    for (auto p : path_) {
        detail += "~";
        detail += StringGraph::NodeIdString(p->out_node_->Id());
    }

    return detail;
}

std::string CompoundPathEdge::ToString() const {
    std::string path_or_edges;
    assert(simple_paths_.size() > 0);
    auto iter = simple_paths_.begin();
    path_or_edges += (*iter)->ToString();
    ++iter;
    for (; iter != simple_paths_.end(); ++iter) {
        path_or_edges += "|";
        path_or_edges += (*iter)->ToString();
    }
    return path_or_edges;
}

std::string CompoundPathEdge::ToDetailString() const {

    std::string detail;
    assert(simple_paths_.size() > 0);
    auto iter = simple_paths_.begin();
    detail += (*iter)->ToString();
    ++iter;
    for (; iter != simple_paths_.end(); ++iter) {
        detail += "|";
        detail += (*iter)->ToString();
    }
    return detail;
}

void PathNode::ReduceEdge(std::vector<PathEdge*> &src, std::vector<PathEdge*> &dst, const PathEdge *e) {
    auto it = std::find(src.begin(), src.end(), e);
    assert(it != src.end());
    dst.push_back(*it);
    src.erase(it);
}


PathGraph::~PathGraph() {
    for (auto& i : nodes_) {
        delete i.second;
    }
    for (auto& i : edges_) {
        delete i.second;
    }
}


void PathGraph::AddEdge(std::list<StringEdge*>& path) {
    assert(path.size() > 0);

    // TODO 循环链

    auto in = nodes_.find(path.front()->in_node_->id_);
    if (in == nodes_.end()) {
        PathNode::ID id = path.front()->in_node_->id_;
        auto r = nodes_.insert(std::make_pair(id, new PathNode(path.front()->in_node_)));
        // TODO
        assert(r.second);
        in = r.first;
    }

    auto out = nodes_.find(path.back()->out_node_->id_);
    if (out == nodes_.end()) {
        PathNode::ID id = path.back()->out_node_->id_;
        auto r = nodes_.insert(std::make_pair(id, new PathNode(path.back()->out_node_)));
        // TODO
        assert(r.second);
        out = r.first;
    }

    PathEdge *e = new SimplePathEdge(in->second, out->second, std::move(path));
    e->type_ = "simple";
    edges_[e->Id()] = e;

    in->second->out_edges_.push_back(e);
    out->second->in_edges_.push_back(e);

}


template<typename T, typename U>
bool BelongTo(const T& src, const U& dst) {
    for (const auto &i : src) {
        if (dst.find(i) == dst.end()) return false;
    }
    return true;
}
/**
 * 选择入度为0得点，找到ego_graph（限制边数）
 * 分析ego_graph的点，是否有来自非ego_graph的点
 * 如果有则找起点和该点的最短路径，当长度小于阈值删除整条路径，
 * 如果中间产生新的入度为0点的，则加入后续点中。
 */
void PathGraph::IdentifyPathSpur(int depth_threshold, int length_threshold) {

    std::unordered_set<PathNode*> candidates;

    for (auto n : nodes_) {
        if (n.second->InDegree() == 0) {
            candidates.insert(n.second);
        }
    }
    
    while (candidates.size() > 0) {
        PathNode* n = *candidates.begin();
        //candidates.erase(candidates.begin());
        bool found = false;
        assert(n->InDegree() == 0);

        std::list<PathNode*> &&ego_nodes = GetEgoNodes(n, depth_threshold, length_threshold*10);    // TODO 
        std::unordered_set<PathNode*> ego_nodes_set(ego_nodes.begin(), ego_nodes.end());


        for (auto b_node : ego_nodes) {
            if (b_node->InDegree() > 1) {
				std::vector<PathNode*> b_in_nodes = b_node->GetInNodes();
                assert(b_in_nodes.size() > 1);

				if (!BelongTo(b_in_nodes, ego_nodes_set)) {

					std::vector<PathEdge*>&& shortest = ShortestPath(n, b_node, ego_nodes_set);
                    std::vector<PathEdge*> vshortest = ReversePath(shortest);
                    int length = std::accumulate(shortest.begin(), shortest.end(), 0, [](int a, const PathEdge *b) { return a + b->Length(); });
                    int vlength = std::accumulate(vshortest.begin(), vshortest.end(), 0, [](int a, const PathEdge *b) { return a + b->Length(); });
					if (length < length_threshold || vlength < length_threshold) {

                        DUMPER["spur.txt"]("%d, %zd, %d, %zd\n", n->Id(), ego_nodes.size(), b_node->Id(), shortest.size()+1);
						for (auto e : shortest) {
                            e->MarkReduce("spur:2", true);
                            ReverseEdge(e)->MarkReduce("spur:2", true);
						}

                        // 添加新产生的候选节点
                        for (auto e : shortest) {
                            assert(std::find(ego_nodes.begin(), ego_nodes.end(), e->out_node_) != ego_nodes.end()); 
                            assert(std::find(ego_nodes.begin(), ego_nodes.end(), e->in_node_) != ego_nodes.end());
                            if (e->out_node_->InDegree() == 0) {
                                candidates.insert(e->out_node_);
                                DUMPER["spur.txt"]("--, %d\n", e->out_node_->Id());
                            }
                        }
                        found = true;
                        break;
					}
				}
            }
        }
        if (!found) candidates.erase(candidates.begin());
    }

}


/** 如果两个节点的有多条长度（子节点数）<=3的的路径，则只保留第一个（排序根据节点名称） */
void PathGraph::RemoveDuplicateSimplePath() {

    std::unordered_map<StringEdge::ID, std::vector<PathEdge*>, StringEdge::Hash, StringEdge::Compare> dup_edges;

    for (const auto &i : edges_) {
        if (i.second->type_ == "simple") {
            SimplePathEdge *e = static_cast<SimplePathEdge*>(i.second);
            if (e->path_.size() < 3) {
                StringEdge::ID id = StringEdge::ID{e->in_node_->id_, e->out_node_->id_};
                auto iter = dup_edges.find(id);
                if (iter != dup_edges.end()) {
                    iter->second.push_back(e);
                }
                else {
                    dup_edges[id] = std::vector<PathEdge*>{ e };
                }
            }
        }
    }

    std::unordered_set<PathEdge*> done;     // 判断对偶路径是否已经处理完成
    for (auto &i : dup_edges) {
        std::vector<PathEdge*> &edges = i.second;
        if (edges.size() > 1 && done.find(edges[0]) == done.end()) {
            done.insert(edges[0]);
            done.insert(ReverseEdge(edges[0]));
            for (size_t j = 1; j < edges.size(); ++j) {
                edges[j]->MarkReduce("simple_dup", true);
                ReverseEdge(edges[j])->MarkReduce("simple_dup", true);
                done.insert(edges[j]);
                done.insert(ReverseEdge(edges[j]));
            }
        }
    }
}

void PathGraph::GetEgoNodes1(PathNode* n, std::unordered_set<PathNode*>& nodes, int depth) {

    nodes.insert(n);
    std::list<PathNode*> candidates{ n };
    int d = 0;
    PathNode* end = n;

    while (d < depth && candidates.size() > 0) {
        auto curr = candidates.front();
        candidates.pop_front();

        for (auto e : curr->out_edges_) {
            if (nodes.find(e->out_node_) == nodes.end()) {
                candidates.push_back(e->out_node_);
                nodes.insert(e->out_node_);
            }

        }

        if (curr == end) {
            d++;
            end = candidates.back();
        }

    }
}

std::list<PathNode*> PathGraph::GetEgoNodes(PathNode* n, int depth_threshold) {
    std::list<PathNode*> nodes{ n };
    std::unordered_set<PathNode*> nodes_set{ n };
    int depth = 0;
    auto curr = nodes.begin();
    auto level_end = nodes.end();
    level_end--;

    while (depth < depth_threshold && curr != nodes.end()) {

        for (auto e: (*curr)->out_edges_) {
            if (nodes_set.find(e->out_node_) == nodes_set.end()) {
                nodes.push_back(e->out_node_);
                nodes_set.insert(e->out_node_);
            }
        }
        
        if (curr == level_end) {
            depth++;
            level_end = nodes.end();
            level_end--;
        }
        curr++;
    }


    return std::list<PathNode*>(nodes.begin(), curr); 
}

std::list<PathNode*> PathGraph::GetEgoNodes(PathNode* n, int depth_threshold, int length_threshold) {
    std::list<PathNode*> nodes{ n }; 
    std::list<int> lens{ 0 };

    int depth = 0;
    auto curr = nodes.begin();
    auto curr_len = lens.begin();

    auto level_end = nodes.end();
    level_end--;

    while (depth < depth_threshold && curr != nodes.end()) {

        for (auto e : (*curr)->out_edges_) {
            if (std::find(nodes.begin(), nodes.end(), e->out_node_) == nodes.end() && 
                *curr_len + e->Length() < length_threshold) {
                nodes.push_back(e->out_node_);
                lens.push_back(*curr_len + e->Length());
            }
        }

        if (curr == level_end) {
            depth++;
            level_end = nodes.end();
            level_end--;
        }
        curr++;
        curr_len++;
    }


    return std::list<PathNode*>(nodes.begin(), curr);
}

std::vector<PathEdge*> PathGraph::ShortestPath(const PathNode* src, const PathNode *dst,
    std::unordered_set<PathNode*> candnodes, int(*score)(PathEdge*)) {

    typedef std::tuple<const PathNode*, PathEdge*, int> WorkType;
    std::vector<WorkType> nodes;
    std::unordered_map<const PathNode*, WorkType> done;

    auto cmp = [](const WorkType &a, const WorkType &b) { return std::get<2>(a) < std::get<2>(b); };
    std::make_heap(nodes.begin(), nodes.end(), cmp);

    nodes.push_back(std::make_tuple(src, (PathEdge*)nullptr, 0));
    std::push_heap(nodes.begin(), nodes.end());

    while (nodes.size() > 0) {
        WorkType i = nodes.front();
        std::pop_heap(nodes.begin(), nodes.end());
        nodes.pop_back();

        if (done.find(std::get<0>(i)) == done.end()) {
            done[std::get<0>(i)] = i;
        
            for (auto e : std::get<0>(i)->out_edges_) {
                if (candnodes.find(e->out_node_) != candnodes.end()) {
                    nodes.push_back(std::make_tuple(e->out_node_, e, std::get<2>(i) + score(e)));
                    std::push_heap(nodes.begin(), nodes.end());
                }
            }
        }
    }

    std::vector<PathEdge*> path;


    auto r = done.find(dst);
    if (r != done.end() && r->first !=src) {
        while (r->first != src) {
            path.push_back(std::get<1>(r->second));
            r = done.find(std::get<1>(r->second)->in_node_);
        }
    }

    std::reverse(path.begin(), path.end());

    return path;

}

CompoundPathEdge* PathGraph::FindBundle(PathNode* start_node, int depth_cutoff, int width_cutoff, int length_cutoff) {
    PathNode* end_node = nullptr;

    std::list<PathEdge*> bundle_edges;      // 存放气泡的边
    std::unordered_map<PathNode*, std::pair<int, int>> visited; // length, score

    std::list<PathNode*> &&local_node_list = GetEgoNodes(start_node, depth_cutoff);
    std::unordered_set<PathNode*> local_nodes(local_node_list.begin(), local_node_list.end());
    std::unordered_set<PathNode*> tips;

    visited[start_node] = std::make_pair(0, 0);
    for (auto e : start_node->out_edges_) {
        tips.insert(e->out_node_);
        bundle_edges.push_back(e);
    }

    int depth = 0;
    double width = 1.0;
    int length = 0;

    bool loop_detect = false;
    bool meet_error = false;
    bool spur = false;

    do {
        std::unordered_map<PathNode*, PathEdge*> new_visited;     // 最新被访问节点，延后加入visited
        std::unordered_set<PathNode*> newtips, oldtips;     // 新产生的末梢节点和未处理的末梢节点

        for (auto n : tips) {
            //if (n->out_edges_.size() == 0) continue;        // dead end

            PathEdge *best_in_edge = nullptr;
            for (auto e : n->in_edges_) {
                // 检查入节点，分成三类：不在局部集合中、已经访问、没有访问
                // 如果所有入节点都已经访问，则找出分数最高的边。并且可以扩展它的出节点
                // 否则改节点延后处理

                if (local_nodes.find(e->in_node_) != local_nodes.end()) {

                    if (visited.find(e->in_node_) != visited.end()) {
                        if (best_in_edge == nullptr || best_in_edge->score_ < e->score_) {
                            best_in_edge = e;
                        }
                    }
                    else {
                        best_in_edge = nullptr;     // 
                        break;
                    }
                }
                else {
                    // 忽略这个入节点
                }
            }

            if (best_in_edge != nullptr) {

                assert(n == best_in_edge->out_node_);
                new_visited[n] = best_in_edge;

                // 如果气泡没有收敛，继续添加新的末梢节点
                if (tips.size() > 1) {
                    for (auto e : n->out_edges_) {
                        if (visited.find(e->out_node_) != visited.end() ||
                            new_visited.find(e->out_node_) != new_visited.end()) {
                            loop_detect = true;
                            break;
                        }

                        PathNode *revese_node = ReverseNode(e->out_node_);
                        if (local_nodes.find(e->out_node_) != local_nodes.end() && 
                            visited.find(revese_node) == visited.end() &&
                            new_visited.find(revese_node) == new_visited.end() ) {

                            if (tips.find(e->out_node_) == tips.end()) {
                                newtips.insert(e->out_node_);
                            }
                            bundle_edges.push_back(e);
                        }
                        else {
                            meet_error = true;
                            break;
                        }
                    }

                    if (n->out_edges_.size() == 0) {
                        spur = true;
                        break;
                    }
                }
                else {
                    end_node = n; // tips[0]
                }
            }
            else {
                oldtips.insert(n);
            }

        }

        for (auto &i : new_visited) {
            visited[i.first] = std::make_pair(
                visited[i.second->in_node_].first + i.second->length_,
                visited[i.second->in_node_].second + i.second->score_);

            // 更新当前长度
            if (length < visited[i.first].first) {
                length = visited[i.first].first;
            }

        }

        depth += 1;
        width = 1.0 * bundle_edges.size() / depth;


        if (tips.size() >= 1) {
            tips.clear();
            tips.insert(newtips.begin(), newtips.end());
            tips.insert(oldtips.begin(), oldtips.end());
        }



    } while (tips.size() >= 1 && tips.size() < 6 && !loop_detect && !meet_error && !spur && depth <= depth_cutoff && length <= length_cutoff && (depth <= 10 || width <= width_cutoff));

    if (end_node != nullptr && !loop_detect && !meet_error && !spur && depth <= depth_cutoff && length <= length_cutoff && (depth <= 10 || width <= width_cutoff)) {
        return new CompoundPathEdge(start_node, end_node, bundle_edges, visited[end_node].first, width, visited[end_node].second);
    }
    else {
        return nullptr;
    }
}


void PathGraph::ConstructCompoundPaths(size_t thread_size) {

    LOG(INFO)("Start FindBubble");
    //std::vector<CompoundPathEdge*>  compound_path0;
    //for (const auto &i : nodes_) {
    //    auto n = i.second;
    //    if (n->out_edges_.size() > 1) {
    //        auto rr =  FindBundle(n);
    //        if (rr != nullptr) {
    //            compound_path0.push_back(std::move(rr));
    //            DUMPER["bubbles.txt"]("%d, %d, %zd\n", rr->in_node_->Id(), rr->out_node_->Id(), rr->simple_paths_.size());
    //        }
    //    }
    //}

    auto work_func = [this](const std::array<decltype(nodes_)::const_iterator, 2> &range) {
        std::vector<CompoundPathEdge*> outputs;
        for (auto i = range[0]; i != range[1]; ++i) {
            auto n = i->second;
            if (n->out_edges_.size() > 1) {
                auto rr =  FindBundle(n);
                if (rr != nullptr) {
                    outputs.push_back(std::move(rr));
                    DUMPER["bubbles.txt"]("%d, %d, %zd\n", rr->in_node_->Id(), rr->out_node_->Id(), rr->simple_paths_.size());
                }
            }
        }
        return outputs;
    };

    std::vector<CompoundPathEdge*>  compound_path0 = MultiThreadRun(thread_size, nodes_, SplitConstIterater<decltype(nodes_)>, work_func, MoveCombineVector<std::vector<CompoundPathEdge*>>);
    LOG(INFO)("End FindBubble %zd", compound_path0.size());

    
    std::sort(compound_path0.begin(), compound_path0.end(), [](const CompoundPathEdge* a, const CompoundPathEdge *b) { return a->simple_paths_.size() > b->simple_paths_.size(); });

    std::unordered_set<PathEdge*> edge_to_cpath;
    std::vector<CompoundPathEdge*>  compound_path1;
    std::unordered_set<PathEdge::ID, PathEdge::Hash, PathEdge::Compare> path1_id;
    for (size_t i = 0; i < compound_path0.size(); ++i) {
        auto path = compound_path0[i];

        bool overlapped = false;
        std::list<PathEdge*> reverse_simple_paths;
        for (auto e : path->simple_paths_) {
            auto re = ReverseEdge(e);
            if (edge_to_cpath.find(e) == edge_to_cpath.end() && 
                edge_to_cpath.find(re) == edge_to_cpath.end()) {
                edge_to_cpath.insert(e);
                edge_to_cpath.insert(re);
                reverse_simple_paths.push_back(re);
            }
            else {
                overlapped = true;
                break;
            }
        }

        if (!overlapped) {
            compound_path1.push_back(path);
            fflush(stdout);
            CompoundPathEdge* reverse_path = new CompoundPathEdge(ReverseNode(path->out_node_), ReverseNode(path->in_node_), reverse_simple_paths, path->length_, path->width_, path->score_);

            compound_path1.push_back(reverse_path);
            path1_id.insert(path->Id());
            path1_id.insert(reverse_path->Id());
        }
        else {
            delete path;
            compound_path0[i] = nullptr;    // TODO 使用智能指针
        }
    }
    

    std::vector<CompoundPathEdge*>  compound_path2;
    // compound_paths_1 -> compound_paths_2  检查镜像节点是否在图中
    // 此步骤需要，FindBubble不能保证能够找镜像
    for (auto &p : compound_path1) {
        if (path1_id.find(PathEdge::ReverseId(p->Id())) != path1_id.end()) {
            compound_path2.push_back(p);
        }  else         {
            delete p;
            p = nullptr;

        }
    }

    // compound_paths_2 -> compound_paths_3 检查每一条edge是否属于compound_path唯一
    // 此检查鱼overlapped检查重复

    // compound_paths_3 -> compound_paths_4  检查镜像节点是否在图中
    // 无需检查必定在图中

    std::unordered_set<PathEdge*> edges_to_reduce;
    for (auto path : compound_path2) {
        for (auto p : path->simple_paths_) {
            edges_to_reduce.insert(p);
        }

        PathEdge* e = path;
        
        edges_[e->Id()] = e;

        e->in_node_->out_edges_.push_back(e);
        e->out_node_->in_edges_.push_back(e);
    }

    for (auto e : edges_to_reduce) {
        e->MarkReduce("contained");
    }
    
    return;
}

void PathGraph::MarkRepeatBridge(int length_threshold) {
    const int LENGTH_THRESHOLD = length_threshold;
    std::unordered_set<PathEdge*> removed;
    for (auto &i : edges_) {
        auto e = i.second;
        if (e->reduce_) continue;
        if (e->in_node_->in_edges_.size() == 1 && e->in_node_->out_edges_.size() >= 2) {
            std::vector<PathEdge*> edge = {e};
            int total_length = edge.back()->Length();
            int total_vlength = ReverseEdge(edge.back())->Length();

            while (total_length < LENGTH_THRESHOLD || total_vlength < LENGTH_THRESHOLD) {
                if (edge.back()->out_node_->in_edges_.size() >= 2 && edge.back()->out_node_->out_edges_.size() == 1) {
                    //removed.insert(edge.begin(), edge.end());
                    removed.insert(edge.front());
                    removed.insert(edge.back());
                    break;

                } else if (edge.back()->out_node_->in_edges_.size() == 1 && edge.back()->out_node_->out_edges_.size() == 1) {
                    edge.push_back(edge.back()->out_node_->out_edges_[0]);
                    total_length += edge.back()->Length();
                    total_vlength += ReverseEdge(edge.back())->Length();
                } else {

                    break;
                }

            }

        }
        /*
        if (e->in_node_->in_edges_.size() == 1 && e->in_node_->out_edges_.size() >= 2 && 
            e->out_node_->in_edges_.size() >= 2 && e->out_node_->out_edges_.size() == 1) {
            
            if (e->Length() < LENGTH_THRESHOLD) {
                removed.insert(e);
            }
        }
        */
    }

    for (auto e : removed) {
        if (!e->reduce_) {  // TODO 是否可以换成断言
            e->MarkReduce("repeat_bridge", true);
            ReverseEdge(e)->MarkReduce("repeat_bridge", true);
        }
                            
    }
}

void PathGraph::IdentifyPaths(const std::string &method) {

    std::list<std::list<PathEdge*>> paths;
    std::unordered_set<PathEdge*> visited;

    for (auto &i : edges_) {
        std::deque<PathEdge*> stack;
        stack.push_back(i.second);

        while (!stack.empty()) {
            PathEdge* e = stack.front();
            stack.pop_front();


            if (!e->reduce_ && visited.find(e) == visited.end()) {
                paths.push_back(ExtendPath(e, visited, method));

                for (auto ie : paths.back().back()->out_node_->out_edges_) {
                    if (!ie->reduce_ && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
                for (auto ie : paths.back().back()->out_node_->in_edges_) {
                    if (!ie->reduce_ && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
                for (auto ie : paths.back().front()->in_node_->out_edges_) {
                    if (!ie->reduce_ && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
                for (auto ie : paths.back().front()->in_node_->in_edges_) {
                    if (!ie->reduce_ && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
            }

        }
    }

    visited.clear();
    for (auto & path : paths) {
        bool overlapped = false;
        for (auto p : path) {
            if (p->reduce_ || visited.find(p) != visited.end() || 
                ReverseEdge(p)->reduce_ || visited.find(ReverseEdge(p)) != visited.end()) {
                overlapped = true;
                break;
            }
        }

        if (!overlapped) {
            paths_.push_back(path);
            std::list<PathEdge*> rpath;
            for (auto p : path) {
                visited.insert(p);
                PathEdge* rp = ReverseEdge(p);
                assert(rp != nullptr);
                visited.insert(rp);
                rpath.push_back(rp);
            }
            std::reverse(rpath.begin(), rpath.end());
            paths_.push_back(rpath);
        }
    }

 
}

template<typename TI, typename TO>
std::list<PathEdge*> PathGraph::ExtendPathWithMethod(PathEdge* e, std::unordered_set<PathEdge*> &visited, TI get_in_edge, TO get_out_edge) {
    assert(visited.find(e) == visited.end());

    std::list<PathEdge*> path;

    std::unordered_set<PathNode*> rnodes;

    visited.insert(e);
    visited.insert(ReverseEdge(e));
    path.push_back(e);
    rnodes.insert(ReverseNode(e->in_node_));
    rnodes.insert(ReverseNode(e->out_node_));


    PathEdge* next = get_out_edge(path.back(), visited);
    while (next != nullptr && rnodes.find(next->out_node_) == rnodes.end()) {
        assert(visited.find(next) == visited.end());
        
        path.push_back(next);
        visited.insert(next);
        visited.insert(ReverseEdge(next));
        rnodes.insert(ReverseNode(next->out_node_));
        next = get_out_edge(next, visited);
    }

    PathEdge* prev = get_in_edge(path.front(), visited);
    while (prev != nullptr && rnodes.find(prev->in_node_) == rnodes.end()) {
        assert(visited.find(prev) == visited.end());

        path.push_front(prev);
        visited.insert(prev);
        visited.insert(ReverseEdge(prev));
        rnodes.insert(ReverseNode(prev->in_node_));
        prev = get_in_edge(prev, visited);
    }
    return path;
}


std::list<PathEdge*> PathGraph::ExtendPath(PathEdge *e, std::unordered_set<PathEdge*> &visited, const std::string &method) {
    assert(visited.find(e) == visited.end());

    if (method == "no") {

        auto get_in_edge = [](const PathEdge *e, const std::unordered_set<PathEdge*> &visited) {
            if (e->in_node_->InDegree() == 1 && e->in_node_->OutDegree() == 1) {
                if (visited.find(e->in_node_->in_edges_.front()) == visited.end()) {
                    return e->in_node_->in_edges_.front();
                }
            }
            return (PathEdge*)nullptr;
        };

        auto get_out_edge = [](const PathEdge *e, const std::unordered_set<PathEdge*> &visited) {
            if (e->out_node_->InDegree() == 1 && e->out_node_->OutDegree() == 1) {
                if (visited.find(e->out_node_->out_edges_.front()) == visited.end()) {
                    return e->out_node_->out_edges_.front();
                }
            }
            return (PathEdge*)nullptr;
        };

        return ExtendPathWithMethod(e, visited, get_in_edge, get_out_edge);
    }
    else if (method == "best") {

        auto get_in_edge = [](const PathEdge* e, const std::unordered_set<PathEdge*> &visited) {
            if (e->in_node_->InDegree() == 1) {
                PathEdge* best_out = e->in_node_->GetBestOutEdge();
                if (best_out == e) {
                    if (visited.find(e->in_node_->in_edges_.front()) == visited.end())
                        return e->in_node_->in_edges_.front();
                }
            }
            return (PathEdge*)nullptr;
        };

        auto get_out_edge = [](const PathEdge* e, const std::unordered_set<PathEdge*> &visited) {
            if (e->out_node_->OutDegree() == 1) {
                PathEdge* best_in = e->out_node_->GetBestInEdge();
                if (best_in == e) {
                    if (visited.find(e->out_node_->out_edges_.front()) == visited.end())
                        return e->out_node_->out_edges_.front();
                }
            }
            return (PathEdge*)nullptr;
        };

        return ExtendPathWithMethod(e, visited, get_in_edge, get_out_edge);
    }
    else {
        LOG(FATAL)("Unknow --select-branch = %s", method.c_str());
        return std::list<PathEdge*>{e};
    }


}


void PathGraph::RemoveCrossEdges() {
    // 如果一个节点同时与反向互补多有
    std::list<PathEdge*> removed;

    for (auto &i : nodes_) {

        auto n = i.second;
        if (n->OutDegree() >= 2) {
            std::vector<PathNode*> out_nodes = n->GetOutNodes();
            std::unordered_set<PathNode*> cands;
            for (auto out : out_nodes) {
                PathNode* rout = ReverseNode(out);
                if (std::find(out_nodes.begin(), out_nodes.end(), rout) != out_nodes.end()) {
                    cands.insert(out);
                    cands.insert(rout);
                }
            }

            if (cands.size() > 0) {
                std::unordered_set<PathNode*> revesed;
                revesed.insert(ReverseNode(n));
                for (auto e : n->in_edges_) {
                    revesed.insert(ReverseNode(e->in_node_));
                    for (auto e1 : e->in_node_->in_edges_) {
                        revesed.insert(ReverseNode(e1->in_node_));

                    }
                }
                for (auto cand : cands) {
                    bool found = false;
                    for (auto e : cand->out_edges_) {
                        if (revesed.find(e->out_node_) == revesed.end()) {
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        for (auto e : n->out_edges_) {
                            if (e->out_node_ == cand) {
                                removed.push_back(e);
                            }
                        }
                    }
                }

            }

        }
    }

    for (auto e : removed) {
        e->MarkReduce("cross");
        ReverseEdge(e)->MarkReduce("cross");
    }
}

void PathGraph::SaveEdges(const std::string &fname) {
    gzFile file = gzopen(fname.c_str(), "w");
    if (file != NULL) {
        for (auto &i : edges_) {
            auto e = i.second;
            auto id = i.first;

            std::string &&detail = e->ToDetailString();


            gzprintf(file, "%s %s %s %s %d %d %s\n", 
                StringGraph::NodeIdString(e->in_node_->id_).c_str(), 
                StringGraph::NodeIdString(std::get<1>(id)).c_str(),
                StringGraph::NodeIdString(e->out_node_->id_).c_str(),
                e->type_.c_str(),
                e->Length(), e->Score(), detail.c_str());

        }
        gzclose(file);
    }
}


void PathGraph::SavePaths(const std::string &fname) {
    FILE *file = fopen(fname.c_str(), "w");
    if (file != NULL) {

        paths_.sort([](const std::list<PathEdge*> &a, const std::list<PathEdge*> &b) { 
            return std::accumulate(a.begin(), a.end(), 0, [](int a, PathEdge* b) { return a + b->score_; }) > 
                   std::accumulate(b.begin(), b.end(), 0, [](int a, PathEdge* b) { return a + b->score_; }); }
        );

        // TODO assert(paths_.size() % 2 == 0);

        int ctgid = 0;
        for (const auto &path : paths_) {

            int score = std::accumulate(path.begin(), path.end(), 0, [](int a, PathEdge* b) { return a + b->score_; });
            int length = std::accumulate(path.begin(), path.end(), 0, [](int a, PathEdge* b) { return a + b->length_; });
            std::string path_string = std::accumulate(path.begin(), path.end(), std::string(""), [](const std::string& a, const PathEdge* b) {
                return (a.empty() ? a : a + "|") + PathEdge::IdToString(b->Id(), "~"); });
      
            fprintf(file, "%06d%c %s %s %s %d %d %s\n",
                ctgid / 2,
                ((ctgid % 2) == 0 ? 'F':'R'),
                path.front()->in_node_ != path.back()->out_node_ ? "ctg_linear" : "ctg_circular",
                PathEdge::IdToString(path.front()->Id(), "~").c_str(),
                StringGraph::NodeIdString(path.back()->out_node_->id_).c_str(),
                length,
                score,
                path_string.c_str()
                );
            ctgid++;
        }
        fclose(file);
    }
}

void PathGraph::SaveCompoundPaths(const std::string &fname) {
    FILE *file = fopen(fname.c_str(), "w");
    if (file != NULL) {

        for (const auto &i : edges_) {
            if (i.second->type_ == "compound") {
                CompoundPathEdge* e = static_cast<CompoundPathEdge*>(i.second);

                fprintf(file, "%s %.01f %d %d %s\n",
                    PathEdge::IdToString(e->Id(), " ").c_str(),
                    e->width_,
                    e->length_,
                    e->score_,
                    e->ToString().c_str()
                );
            }
        }
        fclose(file);
    }
}

size_t PathGraph::PathLength(const std::list<PathEdge*> &path) {
    size_t len =  std::accumulate(path.begin(), path.end(), 0, [](int a, PathEdge* b) { return a + b->score_; });
    return len;
}

void PathGraph::Dump(const std::string& fname) const {
    FILE *file = fopen(fname.c_str(), "w");

    if (file != NULL) {
        for (const auto &i : edges_) {
            PathEdge *e = i.second;
            fprintf(file, "%d, %d, %d, %s\n", e->in_node_->id_, e->out_node_->id_, e->Score(), e->type_.c_str());

        }
        fclose(file);
    }
}

void PathGraph::DumpPaths() const {
    for (const auto &path : paths_) {
        for (auto p : path) {
            printf("%d, ", p->in_node_->id_);
            break;
        }
        printf("%d\n", path.back()->out_node_->id_);
    }   
}

} // namesapce fsa {
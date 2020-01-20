#ifndef PATH_GRAPH_HPP
#define PATH_GRAPH_HPP

#include "string_graph.hpp"
#include <algorithm>
#include <list>

namespace fsa {
class PathEdge;
class PathNode {
    friend class PathEdge;
    friend class PathGraph;
public:
    typedef int ID;

    //PathNode(ID id, StringNode *sn) : id_(id), string_node_(sn) {}
    PathNode(StringNode *n) : id_(n->Id()), string_node_(n) {}

    static ID ReverseId(ID id) {
        return StringNode::ReverseId(id);
    }

    ID Id() const { return id_;  }
    std::vector<PathNode*> GetInNodes();
    std::vector<PathNode*> GetOutNodes();

    void ReduceInEdge(const PathEdge *e) { ReduceEdge(in_edges_, reduced_in_edges_, e); }
    void ReduceOutEdge(const PathEdge *e) { ReduceEdge(out_edges_, reduced_out_edges_, e); }

    void ReduceEdge(std::vector<PathEdge*> &src, std::vector<PathEdge*> &dst, const PathEdge *e);

    size_t InDegree() const {  return in_edges_.size(); }
    size_t OutDegree() const { return out_edges_.size(); }

    PathEdge* GetBestInEdge() const;
    PathEdge* GetBestOutEdge() const;

    ID id_;
    StringNode* string_node_;
    std::vector<PathEdge*> out_edges_;
    std::vector<PathEdge*> in_edges_;
    std::vector<PathEdge*> reduced_out_edges_;
    std::vector<PathEdge*> reduced_in_edges_;
};


class PathEdge {
    friend class PathNode;
    friend class PathGraph;
public:
    using ID = std::array<StringNode::ID, 4>;
    using Hash = ArrayHash<StringNode::ID, 4>; 
    using Compare = ArrayCompare<StringNode::ID, 4>;
    enum ReduceType {
        REPEAT_BRIDGE,
        SPUR,
        SIMPLE_DUP,
        CONTAINED,
        CROSS    
    };
public:
    PathEdge(PathNode* in, PathNode* out) : in_node_(in), out_node_(out) {}
    virtual ~PathEdge() {}

    static ID ReverseId(ID id) {
        return ID{StringNode::ReverseId(id[3]), StringNode::ReverseId(id[2]),
            StringNode::ReverseId(id[1]), StringNode::ReverseId(id[0])};
    }
    virtual ID Id() const = 0;
    virtual std::string ToString() const = 0;
    virtual std::string ToDetailString() const = 0;
    virtual void ReduceStringEdge() = 0;

    static std::string IdToString(ID id, const char *sep=" ");
    int Length() const { return length_; }
    int Score() const { return score_; }
    void MarkReduce(const char* reduce_type, bool remve_string_edge=false) {
        if (!reduce_) {
            reduce_ = true;
            type_ = reduce_type;
            Reduce(remve_string_edge);
        }
    }

    void Reduce(bool remve_string_edge=false) {
        if (reduce_) {
            if (remve_string_edge) ReduceStringEdge();
            out_node_->ReduceInEdge(this);
            in_node_->ReduceOutEdge(this);
        }
    }

    virtual bool Contain(const StringEdge* e) const = 0;
    virtual std::unordered_set<Seq::Id> GetReads() = 0;

    bool reduce_{ false };
    PathNode* in_node_{ nullptr };
    PathNode* out_node_{ nullptr };
    int length_{ 0 };
    double width_{ 0.0 };
    int score_{ 0 };
    std::string type_{ "" };
    
};


class SimplePathEdge : public PathEdge {
public:
    SimplePathEdge(PathNode* in, PathNode* out, std::list<StringEdge*>&& path)
        : PathEdge(in, out), path_(path) {
        score_ = std::accumulate(path.begin(), path.end(), 0, [](int a, StringEdge* b) { return a + b->score_; });
        length_ = std::accumulate(path.begin(), path.end(), 0, [](int a, StringEdge* b) { return a + b->length_; });
        width_ = 1.0;
        type_ = "simple";
    }
    virtual ID Id() const { 
        return ID{path_.front()->in_node_->Id(), path_.front()->out_node_->Id(), 
            path_.back()->in_node_->Id(), path_.back()->out_node_->Id()}; 
    }

    virtual bool Contain(const StringEdge* e) const  {
        return std::find(path_.begin(), path_.end(), e) != path_.end();
    }

    virtual std::unordered_set<Seq::Id> GetReads() {
        std::unordered_set<Seq::Id> reads;
        for (const auto p : path_) {
            reads.insert(p->in_node_->ReadId());
            reads.insert(p->out_node_->ReadId());
        }
        return reads;
    }

    virtual void ReduceStringEdge() {
        for (auto e : path_) {
            e->Reduce();
        }
    }
    virtual std::string ToString() const;
    virtual std::string ToDetailString() const;
    std::list<StringEdge*> path_;
};

class CompoundPathEdge : public PathEdge {
public:
    CompoundPathEdge(PathNode *in, PathNode *out, std::list<PathEdge*> &path, int length, double width, int score)
        : PathEdge(in, out), simple_paths_(path)
    {
        type_ = "compound";
        length_ = length;
        score_ = score;
        width_ = width;
    }
    virtual ID Id() const { return ID{in_node_->id_, 0, 0, out_node_->id_}; }
    virtual std::string ToString() const;
    virtual std::string ToDetailString() const;

    virtual bool Contain(const StringEdge* e) const {
        for (auto p : simple_paths_) {
            if (p->Contain(e)) return true;
        }
        return false;
    }

    virtual std::unordered_set<Seq::Id> GetReads() {
        std::unordered_set<Seq::Id> reads;
        for (const auto p : simple_paths_) {
            auto r = p->GetReads();
            reads.insert(r.begin(), r.end());
        }
        return reads;
    }
    virtual void ReduceStringEdge() {
        for (auto p : simple_paths_) {
            p->ReduceStringEdge();
        }
    }


    std::list<PathEdge*> simple_paths_;

};

class PathGraph {
public:
    virtual ~PathGraph();
     
    void AddEdge(std::list<StringEdge*> &path);

    PathNode* ReverseNode(PathNode* n) {
        return nodes_[PathNode::ReverseId(n->id_)];
    }

    PathEdge* ReverseEdge(PathEdge* e) {
        return edges_[PathEdge::ReverseId(e->Id())];
    }

    std::vector<PathEdge*> ReversePath(std::vector<PathEdge*> &path) {
        std::vector<PathEdge*> vpath(path.size(), nullptr);
        std::transform(path.rbegin(), path.rend(), vpath.begin(), [&](PathEdge* e) { return ReverseEdge(e); });
        return vpath;
    }

    void IdentifyPathSpur(int depth_threshold=10, int length_threshold=50000);
	void RemoveDuplicateSimplePath();
    void GetEgoNodes1(PathNode* n, std::unordered_set<PathNode*>& nodes, int depth);
	std::list<PathNode*> GetEgoNodes(PathNode* n, int depth);
    std::list<PathNode*> GetEgoNodes(PathNode* n, int depth, int length);
    std::vector<PathEdge*> ShortestPath(const PathNode* src, const PathNode *dst, 
        std::unordered_set<PathNode*> candnodes, int(*score)(PathEdge*) = [](PathEdge*) {return 1; });

    CompoundPathEdge*  FindBundle(PathNode* s, int depth_cutoff = 48, int width_cutoff = 16, int length_cutoff = 500000);
    void ConstructCompoundPaths(size_t thread_size);

    void MarkRepeatBridge(int length_threshold=60000);
    void IdentifyPaths(const std::string &method="no");
    std::list<PathEdge*> ExtendPath(PathEdge* e, std::unordered_set<PathEdge*> &visited, const std::string &method);
    
    template<typename TI, typename TO>
    std::list<PathEdge*> ExtendPathWithMethod(PathEdge* e, std::unordered_set<PathEdge*> &visited, TI get_in_edge, TO get_out_edge);

    void RemoveCrossEdges();

    void SortPaths();
    std::list<std::list<PathEdge*>>& GetPaths() { return paths_; }
 
    void SaveEdges(const std::string &fname);
    void SavePaths(const std::string &fname);
    void SaveCompoundPaths(const std::string &fname);

    size_t PathLength(const std::list<PathEdge*> &path);
    
    void Dump(const std::string &fname) const;
    void DumpPaths() const;

protected:

    std::unordered_map<PathEdge::ID, PathEdge*, PathEdge::Hash, PathEdge::Compare> edges_;
    std::unordered_map<PathNode::ID, PathNode*> nodes_;
    std::list<std::list<PathEdge*>> paths_;
};



} // namespace fsa {

#endif // FSA_PATH_GRAPH_HPP  

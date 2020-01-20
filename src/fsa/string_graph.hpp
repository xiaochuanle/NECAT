#ifndef FSA_STRING_GRAPH_HPP
#define FSA_STRING_GRAPH_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <deque>
#include <numeric>
#include <functional>

#include "sequence.hpp"
#include "utility.hpp"

namespace fsa {
class Overlap;


class StringEdge;
class StringNode {
    friend class StringGraph;
    friend class StringEdge;
    friend class PathEdge;
    friend class PathGraph;
public:
    typedef Seq::EndId ID;
    StringNode(ID id) : id_(id) { }

    size_t InDegree() const { return in_edges_.size(); }
    size_t OutDegree() const { return out_edges_.size(); }
    std::vector<StringNode*> GetOutNodes();
    std::vector<StringNode*> GetInNodes();

    std::vector<StringNode*> GetAllOutNodes(); 
    std::vector<StringNode*> GetAllInNodes();

    void ReduceInEdge(const StringEdge *e) { ReduceEdge(in_edges_, reduced_in_edges_, e); }
    void ReduceOutEdge(const StringEdge *e) { ReduceEdge(out_edges_, reduced_out_edges_, e); }

    ID Id() const { return id_;  }

    static ID ReverseId(ID id) { return -id;  }

    int ReadId() const { return Seq::EndIdToId(id_);  }
    StringEdge* GetBestInEdge() const { return best_in_;  }
    StringEdge* GetBestOutEdge() const { return best_out_; }
protected:
    void ReduceEdge(std::vector<StringEdge*> &src, std::vector<StringEdge*> &dst, const StringEdge *e);
public:
    Seq::EndId id_{ 0 };
    std::vector<StringEdge*> out_edges_;
    std::vector<StringEdge*> in_edges_;
    std::vector<StringEdge*> reduced_out_edges_;
    std::vector<StringEdge*> reduced_in_edges_;

    int mark_;
    StringEdge* best_in_{ nullptr };
    StringEdge* best_out_{ nullptr };
};


class StringEdge {
    friend class StringGraph;
    friend class StringNode;

public:
    enum Type {
        ACTIVE,
        SPUR,
        CHIMER,
        CHIMER1,
        REMOVED,
        TRANSITIVE,
        NO_BEST,
    };
    using Hash = ArrayHash<int, 2>;
    using Compare = ArrayCompare<int, 2>;
    using ID = std::array<StringNode::ID, 2>;

    StringEdge(StringNode* in_node, StringNode *out_node) {
        reduce_ = false;
        out_node_ = out_node;
        in_node_ = in_node;
        length_ = 0;
        score_ = 0;
    }

    static ID ReverseID(ID id) {
        return ID{StringNode::ReverseId(id[1]), StringNode::ReverseId(id[0])};
    }

    ID Id() {
        return ID{in_node_->id_, out_node_->id_};
    }

    std::tuple<int, bool, int, int> GetSeqArea() const;
    void Reduce(Type t=REMOVED) {
        reduce_ = true;
        out_node_->ReduceInEdge(this);
        in_node_->ReduceOutEdge(this);
        type_ = t;
    }

    bool reduce_{ false };
    StringNode* out_node_{ nullptr };
    StringNode* in_node_{ nullptr };
    int read_{ 0 };
    int start_{ 0 };
    int end_{ 0 };
    int length_{ 0 };
    int score_{ 0 };
    double identity_{ 0 };
    Type type_{ ACTIVE };
};


class StringGraph {
public:
    virtual ~StringGraph();
public:
    static std::string NodeIdString(int id);
    static std::string ReadIdString(int id);

	static int ReverseNode(int id) {
		return -id;
	}
	static StringEdge::ID ReverseEdge(StringEdge::ID id) {
		return StringEdge::ID{ReverseNode(id[1]), ReverseNode(id[0])};
	}
	StringNode* ReverseNode(StringNode* n) {
		return nodes_[ReverseNode(n->id_)];
	}

	StringEdge* ReverseEdge(StringEdge* e) {
		return edges_[ReverseEdge(e->Id())];
	}
    
    std::list<StringEdge*> Reverse(const std::list<StringEdge*>& path);

    StringNode* GetNode(StringNode::ID id) {
        auto i = nodes_.find(id);
        return i != nodes_.end() ? i->second : nullptr;
    }

	void AddOverlap(const Overlap* overlap);
	void AddOverlaps(const std::deque<Overlap> &ovlps, int min_length, int min_aligned_lenght, float min_identity);
    bool FilterOverlap(const Overlap &ovlp, const std::unordered_set<StringNode::ID> &contained, int min_length, int min_aligned_length, float min_identity);
	void AddEdge(int in_node, int out_node, int len, int score, double identity, int read, int start, int end);

	void MarkTransitiveEdges();

	void MarkChimerEdges();
	void MarkSpurEdges();
	std::unordered_set<StringNode*> BfsNodes(StringNode* n, StringNode* exclude=nullptr, int depth=5);
	
	void MarkBestOverlap();
	void ResolveRepeatEdges();

	void IdentifySimplePaths();
    std::list<StringEdge*> ExtendSimplePath(StringEdge*n, std::unordered_set<StringEdge*> &visited);

        
    std::list<std::list<StringEdge*>>& GetPaths() {  return paths_; }

    std::vector<StringEdge*> ShortestPath(const StringNode* src, const StringNode *dst, std::unordered_set<StringEdge*> edges, int(*score)(StringEdge*) = [](StringEdge*) {return 1; });


    // To check whether path corresponds to a reversed path.
    bool Assert_PathDual(const std::list<std::list<StringEdge*>> paths);

	void Dump();
	void DumpPaths();
    void SaveChimerNode(const std::string &fname);
    void SaveEdges(const std::string &fname);

protected:
	std::unordered_map<StringNode::ID, StringNode*> nodes_;
	std::unordered_map<StringEdge::ID, StringEdge*, StringEdge::Hash, StringEdge::Compare> edges_;

    std::list<std::list<StringEdge*>> paths_;
};

} //namespace fsa {

#endif // FSA_STRING_GRAPH_HPP  

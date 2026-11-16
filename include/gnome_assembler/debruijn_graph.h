#pragma once

#include "gnome_assembler/kmer_counter.h"

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace gnome {

// A directed edge in the de Bruijn graph.
// Represents a k-mer: the edge goes from prefix (k-1 bases) to suffix (k-1 bases).
struct GraphEdge {
    std::string from;        // (k-1)-mer prefix node
    std::string to;          // (k-1)-mer suffix node
    uint32_t multiplicity;   // number of times this k-mer was observed (after filtering)
    bool used = false;       // for future Eulerian traversal (Milestone 3)
};

// Degree information for a single node.
struct NodeDegree {
    uint32_t in_degree = 0;
    uint32_t out_degree = 0;
};

// Classification of a node's role in the graph.
enum class NodeKind {
    Normal,     // in-degree == out-degree (balanced)
    StartNode,  // out-degree > in-degree by 1 (Euler path start)
    EndNode,    // in-degree > out-degree by 1 (Euler path end)
    Branch,     // out-degree > 1 (branching point)
    Tip,        // in-degree == 0 OR out-degree == 0 (but not start/end by one)
};

NodeKind classify_node(const NodeDegree& degree);

// The complete de Bruijn graph built from a filtered KmerTable.
struct DeBruijnGraph {
    int k = 0;                                               // k-mer length used to build this graph
    std::vector<GraphEdge> edges;                            // all directed edges
    std::unordered_map<std::string, NodeDegree> nodes;      // all nodes with their degrees

    std::size_t component_count = 0;   // number of weakly connected components
    std::size_t tip_count = 0;         // nodes with in=0 or out=0
    std::size_t branch_count = 0;      // nodes with out>1
    bool is_eulerian_path = false;     // exactly one start, one end, rest balanced
    bool is_eulerian_cycle = false;    // all nodes balanced and one component

    // Returns the edge index for a given (from, to) pair, or -1 if absent.
    int find_edge(const std::string& from, const std::string& to) const;
};

// Build a de Bruijn graph from a filtered KmerTable.
DeBruijnGraph build_graph(const KmerTable& kmers, int k);

// Recompute component_count, tip_count, branch_count, is_eulerian_path, is_eulerian_cycle.
void analyse_graph(DeBruijnGraph& graph);

// Write a human-readable deterministic dump of the graph to an output stream.
// Output is sorted for reproducibility: nodes alphabetically, edges by (from, to).
void dump_graph(const DeBruijnGraph& graph, std::ostream& out);

} // namespace gnome

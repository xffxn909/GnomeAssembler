#include "gnome_assembler/debruijn_graph.h"

#include <algorithm>
#include <queue>
#include <unordered_set>
#include <sstream>

namespace gnome {

// ---------------------------------------------------------------------------
// NodeKind classification
// ---------------------------------------------------------------------------

NodeKind classify_node(const NodeDegree& degree) {
    const uint32_t in = degree.in_degree;
    const uint32_t out = degree.out_degree;

    if (in == out) {
        return NodeKind::Normal;
    }
    if (in == 0 || out == 0) {
        return NodeKind::Tip;
    }
    if (out == in + 1) {
        return NodeKind::StartNode;
    }
    if (in == out + 1) {
        return NodeKind::EndNode;
    }
    if (out > 1) {
        return NodeKind::Branch;
    }
    return NodeKind::Normal;
}

// ---------------------------------------------------------------------------
// find_edge
// ---------------------------------------------------------------------------

int DeBruijnGraph::find_edge(const std::string& from, const std::string& to) const {
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        if (edges[static_cast<std::size_t>(i)].from == from &&
            edges[static_cast<std::size_t>(i)].to == to) {
            return i;
        }
    }
    return -1;
}

// ---------------------------------------------------------------------------
// build_graph
// ---------------------------------------------------------------------------

DeBruijnGraph build_graph(const KmerTable& kmers, int k) {
    DeBruijnGraph graph;
    graph.k = k;

    const int prefix_len = k - 1;

    for (const auto& [kmer, count] : kmers) {
        if (static_cast<int>(kmer.size()) != k) {
            continue;
        }
        const std::string from = kmer.substr(0, static_cast<std::size_t>(prefix_len));
        const std::string to   = kmer.substr(1, static_cast<std::size_t>(prefix_len));

        // Ensure nodes exist.
        graph.nodes[from]; // default-constructs NodeDegree{0,0} if absent
        graph.nodes[to];

        // Find an existing edge or create a new one.
        const int existing = graph.find_edge(from, to);
        if (existing >= 0) {
            graph.edges[static_cast<std::size_t>(existing)].multiplicity += count;
        } else {
            graph.edges.push_back({from, to, count, false});
        }

        // Update degrees.
        graph.nodes[from].out_degree += count;
        graph.nodes[to].in_degree   += count;
    }

    analyse_graph(graph);
    return graph;
}

// ---------------------------------------------------------------------------
// analyse_graph — weakly-connected components via BFS on undirected adjacency
// ---------------------------------------------------------------------------

void analyse_graph(DeBruijnGraph& graph) {
    // Build an adjacency list (undirected) over node labels.
    std::unordered_map<std::string, std::vector<std::string>> adj;
    for (const auto& e : graph.edges) {
        adj[e.from].push_back(e.to);
        adj[e.to].push_back(e.from);
    }

    // BFS to count weakly connected components.
    std::unordered_set<std::string> visited;
    std::size_t components = 0;
    for (const auto& [node, _] : graph.nodes) {
        if (visited.count(node) == 0) {
            ++components;
            std::queue<std::string> q;
            q.push(node);
            visited.insert(node);
            while (!q.empty()) {
                const std::string current = q.front();
                q.pop();
                for (const std::string& neighbour : adj[current]) {
                    if (visited.count(neighbour) == 0) {
                        visited.insert(neighbour);
                        q.push(neighbour);
                    }
                }
            }
        }
    }
    graph.component_count = components;

    // Count tips, branches, and check Eulerian conditions.
    std::size_t tips = 0;
    std::size_t branches = 0;
    std::size_t start_nodes = 0;
    std::size_t end_nodes = 0;
    std::size_t balanced_nodes = 0;
    std::size_t imbalanced_nodes = 0;

    for (const auto& [name, degree] : graph.nodes) {
        if (degree.in_degree == 0 || degree.out_degree == 0) {
            ++tips;
        }
        if (degree.out_degree > 1) {
            ++branches;
        }

        if (degree.in_degree == degree.out_degree) {
            ++balanced_nodes;
        } else if (degree.out_degree == degree.in_degree + 1) {
            ++start_nodes;
        } else if (degree.in_degree == degree.out_degree + 1) {
            ++end_nodes;
        } else {
            ++imbalanced_nodes;
        }
    }

    graph.tip_count    = tips;
    graph.branch_count = branches;

    // Eulerian cycle: all balanced, single component, at least one edge.
    graph.is_eulerian_cycle = (components == 1) && !graph.edges.empty() && (balanced_nodes == graph.nodes.size());

    // Eulerian path: exactly one start + one end, all other nodes balanced, single component.
    graph.is_eulerian_path =
        (components == 1) && !graph.edges.empty() &&
        start_nodes == 1 &&
        end_nodes   == 1 &&
        imbalanced_nodes == 0;
}

// ---------------------------------------------------------------------------
// dump_graph — deterministic human-readable output
// ---------------------------------------------------------------------------

void dump_graph(const DeBruijnGraph& graph, std::ostream& out) {
    out << "GnomeAssembler de Bruijn graph dump\n";
    out << "k = " << graph.k << "\n";
    out << "Nodes: " << graph.nodes.size() << "\n";
    out << "Edges: " << graph.edges.size() << "\n";
    out << "Components: " << graph.component_count << "\n";
    out << "Tips: " << graph.tip_count << "\n";
    out << "Branches: " << graph.branch_count << "\n";
    out << "Eulerian cycle: " << (graph.is_eulerian_cycle ? "yes" : "no") << "\n";
    out << "Eulerian path: "  << (graph.is_eulerian_path  ? "yes" : "no") << "\n";
    out << "\n";

    // Sorted node list for reproducibility.
    std::vector<std::string> node_names;
    node_names.reserve(graph.nodes.size());
    for (const auto& [name, _] : graph.nodes) {
        node_names.push_back(name);
    }
    std::sort(node_names.begin(), node_names.end());

    out << "--- Nodes ---\n";
    for (const std::string& name : node_names) {
        const NodeDegree& d = graph.nodes.at(name);
        const NodeKind kind = classify_node(d);
        const char* kind_str = "normal";
        switch (kind) {
            case NodeKind::Tip:       kind_str = "tip";       break;
            case NodeKind::Branch:    kind_str = "branch";    break;
            case NodeKind::StartNode: kind_str = "start";     break;
            case NodeKind::EndNode:   kind_str = "end";       break;
            case NodeKind::Normal:    break;
        }
        out << name << "  in=" << d.in_degree << "  out=" << d.out_degree
            << "  [" << kind_str << "]\n";
    }

    // Sorted edge list.
    std::vector<const GraphEdge*> edge_ptrs;
    edge_ptrs.reserve(graph.edges.size());
    for (const auto& e : graph.edges) {
        edge_ptrs.push_back(&e);
    }
    std::sort(edge_ptrs.begin(), edge_ptrs.end(), [](const GraphEdge* a, const GraphEdge* b) {
        if (a->from != b->from) return a->from < b->from;
        return a->to < b->to;
    });

    out << "\n--- Edges ---\n";
    for (const GraphEdge* e : edge_ptrs) {
        out << e->from << " -> " << e->to << "  x" << e->multiplicity << "\n";
    }
}

} // namespace gnome

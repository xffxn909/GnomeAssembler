#include "gnome_assembler/assembly.h"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gnome {

namespace {

std::string invalid_graph_message(const DeBruijnGraph& graph) {
    if (graph.edges.empty()) {
        return "Cannot assemble an empty graph: no retained k-mers are available.";
    }
    if (graph.component_count != 1) {
        return "Cannot assemble graph: it has " + std::to_string(graph.component_count) +
               " disconnected components; a single complete contig requires one component.";
    }

    std::size_t starts = 0;
    std::size_t ends = 0;
    for (const auto& [name, degree] : graph.nodes) {
        (void)name;
        if (degree.out_degree == degree.in_degree + 1) ++starts;
        else if (degree.in_degree == degree.out_degree + 1) ++ends;
    }
    return "Cannot assemble graph: it is not Eulerian (" + std::to_string(starts) +
           " start-like and " + std::to_string(ends) +
           " end-like nodes). Adjust k-mer filtering or inspect the graph dump.";
}

} // namespace

AssemblyResult assemble_graph(const DeBruijnGraph& graph) {
    AssemblyResult result;
    if (!graph.is_eulerian_path && !graph.is_eulerian_cycle) {
        result.error = invalid_graph_message(graph);
        return result;
    }

    std::string start;
    if (graph.is_eulerian_path) {
        for (const auto& [name, degree] : graph.nodes) {
            if (degree.out_degree == degree.in_degree + 1) {
                start = name;
                break;
            }
        }
    } else {
        for (const auto& [name, degree] : graph.nodes) {
            if (start.empty() || name < start) start = name;
        }
    }

    struct TraversalEdge { std::string to; std::size_t index; };
    std::unordered_map<std::string, std::vector<TraversalEdge>> adjacency;
    std::vector<uint32_t> remaining;
    remaining.reserve(graph.edges.size());
    std::size_t total_edges = 0;
    for (std::size_t index = 0; index < graph.edges.size(); ++index) {
        const GraphEdge& edge = graph.edges[index];
        adjacency[edge.from].push_back({edge.to, index});
        remaining.push_back(edge.multiplicity);
        total_edges += edge.multiplicity;
    }
    for (auto& [node, edges] : adjacency) {
        (void)node;
        std::sort(edges.begin(), edges.end(), [](const TraversalEdge& left,
                                                 const TraversalEdge& right) {
            if (left.to != right.to) return left.to > right.to;
            return left.index > right.index;
        });
    }

    // Each multiplicity unit is an independently consumable parallel edge.
    std::vector<std::string> stack;
    std::vector<std::string> reversed_vertices;
    stack.push_back(start);
    while (!stack.empty()) {
        const std::string& node = stack.back();
        auto& outgoing = adjacency[node];
        while (!outgoing.empty() && remaining[outgoing.back().index] == 0) {
            outgoing.pop_back();
        }
        if (outgoing.empty()) {
            reversed_vertices.push_back(node);
            stack.pop_back();
        } else {
            const TraversalEdge edge = outgoing.back();
            --remaining[edge.index];
            ++result.traversed_edges;
            stack.push_back(edge.to);
        }
    }

    if (result.traversed_edges != total_edges ||
        reversed_vertices.size() != total_edges + 1) {
        result.error = "Traversal stopped before every graph edge was consumed.";
        result.contigs.clear();
        return result;
    }

    std::reverse(reversed_vertices.begin(), reversed_vertices.end());
    AssemblyContig contig;
    contig.circular = graph.is_eulerian_cycle;
    contig.traversed_edges = total_edges;
    contig.sequence = reversed_vertices.front();
    for (std::size_t i = 1; i < reversed_vertices.size(); ++i) {
        contig.sequence.push_back(reversed_vertices[i].back());
    }
    // The final k-1 bases in a cycle only close the loop and must not be
    // reported as extra biological sequence.
    if (contig.circular) {
        contig.sequence.resize(total_edges);
    }
    result.contigs.push_back(std::move(contig));

    if (graph.branch_count != 0) {
        result.warnings.push_back(
            "The graph contains branches or repeated edges; this traversal is one deterministic Eulerian ordering, not proof of a unique assembly.");
    }
    return result;
}

} // namespace gnome

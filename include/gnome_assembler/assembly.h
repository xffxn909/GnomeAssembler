#pragma once

#include "gnome_assembler/debruijn_graph.h"

#include <cstddef>
#include <string>
#include <vector>

namespace gnome {

struct AssemblyContig {
    std::string sequence;
    bool circular = false;
    std::size_t traversed_edges = 0;
};

// Result of traversing a complete Eulerian de Bruijn graph.  A non-empty
// error means no contig was produced; warnings describe ambiguity without
// overstating the certainty of an otherwise valid traversal.
struct AssemblyResult {
    std::vector<AssemblyContig> contigs;
    std::vector<std::string> warnings;
    std::string error;
    std::size_t traversed_edges = 0;
};

// Traverse a single connected Eulerian graph with iterative Hierholzer's
// algorithm. Edge multiplicity is treated as parallel graph edges.
AssemblyResult assemble_graph(const DeBruijnGraph& graph);

} // namespace gnome

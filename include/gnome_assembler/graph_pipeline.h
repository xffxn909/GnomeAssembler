#pragma once

#include "gnome_assembler/cleaning_options.h"
#include "gnome_assembler/debruijn_graph.h"
#include "gnome_assembler/kmer_counter.h"
#include "gnome_assembler/kmer_options.h"

#include <string>

namespace gnome {

// Options for building the graph (cleaning + k-mer settings together).
struct GraphPipelineOptions {
    CleaningOptions cleaning;
    KmerOptions kmer;
    std::string output_prefix;  // may be empty; if set, a .graph-dump.txt is written
};

// Summary produced by the graph pipeline.
struct GraphPipelineResult {
    KmerStatistics kmer_stats;
    DeBruijnGraph graph;
    std::string summary;
    std::string error;
};

// Run the full pipeline: clean the FASTQ, extract & filter k-mers, build the graph.
// If output_prefix is set, writes <prefix>.graph-dump.txt alongside the cleaning outputs.
GraphPipelineResult run_graph_pipeline(const GraphPipelineOptions& options);

// Format a human-readable summary of the graph-build result.
std::string format_graph_summary(const KmerStatistics& kmer_stats,
                                 const DeBruijnGraph& graph,
                                 const KmerOptions& kmer_options);

} // namespace gnome

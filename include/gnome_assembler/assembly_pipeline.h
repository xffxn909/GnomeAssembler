#pragma once

#include "gnome_assembler/assembly.h"
#include "gnome_assembler/graph_pipeline.h"
#include "gnome_assembler/verification.h"
#include "gnome_assembler/blast_identification.h"

#include <string>

namespace gnome {

struct AssemblyPipelineResult {
    GraphPipelineResult graph_result;
    AssemblyResult assembly;
    VerificationReport verification;
    BlastResult blast;
    std::string summary;
    std::string error;
};

// Clean reads, build a graph, traverse it, and optionally write contigs to
// <output-prefix>.contigs.fasta.
AssemblyPipelineResult run_assembly_pipeline(const GraphPipelineOptions& options,
                                             const BlastOptions& blast_options = {});

std::string format_assembly_summary(const GraphPipelineResult& graph_result,
                                    const AssemblyResult& assembly);

} // namespace gnome

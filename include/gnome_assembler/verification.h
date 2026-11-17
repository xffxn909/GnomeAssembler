#pragma once

#include "gnome_assembler/assembly.h"
#include "gnome_assembler/debruijn_graph.h"

#include <cstddef>
#include <string>
#include <vector>

namespace gnome {

struct VerificationReport {
    std::size_t contig_count = 0;
    std::size_t total_contig_length = 0;
    std::size_t longest_contig_length = 0;
    std::size_t n50 = 0;
    std::size_t cleaned_reads = 0;
    std::size_t covered_reads = 0;
    std::size_t graph_edges = 0;
    std::size_t traversed_edges = 0;
    bool all_edges_used = false;
    bool assembly_complete = false;
    std::vector<std::string> warnings;
};

// Verify an assembly against its graph and the cleaned reads used to build it.
// A read is covered when it occurs in a contig or its reverse complement; for
// circular contigs, matches may span the reported sequence boundary.
VerificationReport verify_assembly(const DeBruijnGraph& graph,
                                   const AssemblyResult& assembly,
                                   const std::vector<std::string>& cleaned_reads);

std::string format_verification_text(const VerificationReport& report);
std::string format_verification_json(const VerificationReport& report);

} // namespace gnome

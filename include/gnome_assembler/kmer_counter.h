#pragma once

#include "gnome_assembler/kmer_options.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace gnome {

// A table of k-mer -> count, built from one or more sequences.
using KmerTable = std::unordered_map<std::string, uint32_t>;

// Count all k-mers of length k in a single sequence.
// K-mers that contain any character other than A, C, G, T are skipped.
void count_kmers(const std::string& sequence, int k, KmerTable& table);

// Remove all entries whose count is strictly below min_abundance.
void filter_kmers(KmerTable& table, uint32_t min_abundance);

// Summary statistics produced alongside a KmerTable.
struct KmerStatistics {
    std::size_t total_kmers_seen = 0;       // sum of all raw counts before filtering
    std::size_t distinct_kmers_before = 0;  // distinct entries before filtering
    std::size_t distinct_kmers_after = 0;   // distinct entries after filtering
};

// Build and filter a KmerTable from a list of sequences.
// sequences must already be uppercase and unambiguous is acceptable to skip.
KmerTable build_kmer_table(const std::vector<std::string>& sequences,
                           const KmerOptions& options,
                           KmerStatistics& statistics);

} // namespace gnome

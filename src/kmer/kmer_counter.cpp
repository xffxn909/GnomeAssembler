#include "gnome_assembler/kmer_counter.h"

#include <algorithm>
#include <cctype>

namespace gnome {

namespace {

bool is_unambiguous_base(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

bool kmer_is_unambiguous(const std::string& sequence, std::size_t start, int k) {
    for (int i = 0; i < k; ++i) {
        if (!is_unambiguous_base(sequence[start + static_cast<std::size_t>(i)])) {
            return false;
        }
    }
    return true;
}

} // namespace

void count_kmers(const std::string& sequence, int k, KmerTable& table) {
    if (k <= 0 || static_cast<std::size_t>(k) > sequence.size()) {
        return;
    }
    const std::size_t limit = sequence.size() - static_cast<std::size_t>(k);
    for (std::size_t i = 0; i <= limit; ++i) {
        if (kmer_is_unambiguous(sequence, i, k)) {
            table[sequence.substr(i, static_cast<std::size_t>(k))]++;
        }
    }
}

void filter_kmers(KmerTable& table, uint32_t min_abundance) {
    for (auto it = table.begin(); it != table.end(); ) {
        if (it->second < min_abundance) {
            it = table.erase(it);
        } else {
            ++it;
        }
    }
}

KmerTable build_kmer_table(const std::vector<std::string>& sequences,
                           const KmerOptions& options,
                           KmerStatistics& statistics) {
    KmerTable table;
    for (const std::string& seq : sequences) {
        count_kmers(seq, options.kmer_length, table);
    }

    // Gather stats before filtering.
    statistics.distinct_kmers_before = table.size();
    statistics.total_kmers_seen = 0;
    for (const auto& entry : table) {
        statistics.total_kmers_seen += entry.second;
    }

    filter_kmers(table, options.min_kmer_abundance);
    statistics.distinct_kmers_after = table.size();

    return table;
}

} // namespace gnome

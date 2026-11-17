#pragma once

#include <string>
#include <vector>

namespace gnome {

struct BlastOptions {
    std::string executable_path;
    std::string database_path;
    double minimum_identity = 70.0;
    double minimum_query_coverage = 50.0;
};

struct BlastHit {
    std::string query_id;
    std::string subject_id;
    std::string title;
    double identity_percent = 0.0;
    double query_coverage_percent = 0.0;
    double evalue = 0.0;
    double bit_score = 0.0;
    std::size_t alignment_length = 0;
    std::size_t query_length = 0;
};

struct BlastResult {
    bool requested = false;
    bool identified = false;
    std::string status;
    std::string error;
    std::vector<BlastHit> qualified_hits;
    BlastHit best_hit;
};

bool blast_options_requested(const BlastOptions& options);
BlastResult identify_contigs_with_blast(const std::string& contigs_fasta,
                                        const BlastOptions& options,
                                        const std::string& tabular_output_path);
std::vector<BlastHit> parse_blast_tabular(const std::string& tabular_text,
                                          std::string& error);
std::vector<BlastHit> filter_and_rank_blast_hits(const std::vector<BlastHit>& hits,
                                                 const BlastOptions& options);
std::string format_blast_report_text(const BlastResult& result, const BlastOptions& options);
std::string format_blast_report_json(const BlastResult& result, const BlastOptions& options);

} // namespace gnome

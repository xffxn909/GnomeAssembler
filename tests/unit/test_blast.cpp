#include "gnome_assembler/blast_identification.h"

#include <iostream>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void test_parser_and_ranked_qualification() {
    const std::string fixture =
        "contig_1\tweak\t99.1\t80\t100\t1000\t1e-20\t150\tWeaker hit\n"
        "contig_1\tstrong\t95.0\t90\t100\t900\t1e-10\t200\tBest bit score\n"
        "contig_1\tlow_coverage\t99.9\t20\t100\t1000\t1e-30\t300\tFiltered coverage\n";
    std::string error;
    const auto hits = gnome::parse_blast_tabular(fixture, error);
    check(error.empty(), "Valid BLAST fixture should parse.");
    check(hits.size() == 3, "All fixture hits should parse.");
    if (hits.size() == 3) {
        check(hits[0].query_coverage_percent == 80.0, "Coverage should use alignment/query length.");
    }

    gnome::BlastOptions thresholds;
    thresholds.minimum_identity = 90.0;
    thresholds.minimum_query_coverage = 50.0;
    const auto qualified = gnome::filter_and_rank_blast_hits(hits, thresholds);
    check(qualified.size() == 2, "Only hits meeting identity and coverage thresholds should remain.");
    if (!qualified.empty()) {
        check(qualified.front().subject_id == "strong", "Highest bit score should select the best qualified hit.");
    }

    // The public execution function verifies paths before spawning. The parser
    // and report format are separately deterministic and testable without BLAST+.
    gnome::BlastResult result;
    result.requested = true;
    result.identified = true;
    result.status = "Identified: strong.";
    result.best_hit = hits[1];
    gnome::BlastOptions options;
    const std::string report = gnome::format_blast_report_text(result, options);
    const std::string json = gnome::format_blast_report_json(result, options);
    check(report.find("Best subject: strong") != std::string::npos,
          "Text report should expose selected hit metrics.");
    check(json.find("\"subject_id\": \"strong\"") != std::string::npos,
          "JSON report should preserve selected hit identity.");
}

void test_failures_are_explicit() {
    gnome::BlastOptions options;
    options.executable_path = "missing-blastn.exe";
    options.database_path = "missing-db";
    const auto result = gnome::identify_contigs_with_blast("missing-contigs.fasta", options,
                                                            "missing-output.tsv");
    check(result.requested, "Configured BLAST paths should mark identification requested.");
    check(!result.identified && !result.error.empty(), "Unavailable executable should have an explicit failure.");

    std::string error;
    const auto malformed = gnome::parse_blast_tabular("too\tfew\tfields\n", error);
    check(malformed.empty() && !error.empty(), "Malformed tabular BLAST data should be rejected.");
}

} // namespace

int main() {
    test_parser_and_ranked_qualification();
    test_failures_are_explicit();
    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed.\n";
        return 1;
    }
    std::cout << "All Milestone 5 BLAST tests passed.\n";
    return 0;
}

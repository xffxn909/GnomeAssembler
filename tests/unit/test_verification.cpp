#include "gnome_assembler/assembly.h"
#include "gnome_assembler/debruijn_graph.h"
#include "gnome_assembler/verification.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void test_complete_linear_assembly() {
    gnome::KmerTable kmers{{"ATC", 1}, {"TCG", 1}, {"CGG", 1}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto assembly = gnome::assemble_graph(graph);
    const std::vector<std::string> reads{"ATCG", "CGG", "CCG"}; // CCG is reverse-complement of CGG.
    const auto report = gnome::verify_assembly(graph, assembly, reads);

    check(report.assembly_complete, "All graph edges and reads should establish completeness.");
    check(report.all_edges_used, "Every retained graph edge should be marked used.");
    check(report.graph_edges == 3 && report.traversed_edges == 3, "Edge totals should agree.");
    check(report.contig_count == 1, "One contig should be reported.");
    check(report.total_contig_length == 5 && report.longest_contig_length == 5,
          "Contig length metrics should be correct.");
    check(report.n50 == 5, "N50 of one 5bp contig should be 5.");
    check(report.covered_reads == 3 && report.cleaned_reads == 3,
          "Forward and reverse-complement reads should be covered.");
}

void test_circular_boundary_coverage() {
    gnome::KmerTable kmers{{"ATC", 1}, {"TCA", 1}, {"CAT", 1}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto assembly = gnome::assemble_graph(graph);
    // TCA crosses the deterministic ATC reporting boundary. ATCAT spans the
    // boundary more than once and must also be recognized on a circular contig.
    const auto report = gnome::verify_assembly(graph, assembly, {"TCA", "ATCAT"});
    check(report.assembly_complete, "Reads crossing a circular boundary should be covered.");
    check(report.covered_reads == 2, "Circular boundary reads should count as covered.");
}

void test_incomplete_assembly_warning_and_json() {
    gnome::KmerTable kmers{{"ACG", 1}, {"ACT", 1}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto assembly = gnome::assemble_graph(graph);
    const auto report = gnome::verify_assembly(graph, assembly, {"ACG"});
    const std::string json = gnome::format_verification_json(report);

    check(!report.assembly_complete, "Failed traversal must not be reported complete.");
    check(!report.all_edges_used, "Failed traversal must report unused edges.");
    check(!report.warnings.empty(), "Incomplete assembly must explain the problem.");
    check(json.find("\"assembly_complete\": false") != std::string::npos,
          "JSON report should contain completion state.");
    check(json.find("\"warnings\": [") != std::string::npos,
          "JSON report should contain warnings.");
}

} // namespace

int main() {
    test_complete_linear_assembly();
    test_circular_boundary_coverage();
    test_incomplete_assembly_warning_and_json();
    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed.\n";
        return 1;
    }
    std::cout << "All Milestone 4 verification tests passed.\n";
    return 0;
}

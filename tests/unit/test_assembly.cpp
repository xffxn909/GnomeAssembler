#include "gnome_assembler/assembly.h"
#include "gnome_assembler/debruijn_graph.h"

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

void test_linear_contig() {
    // ATC -> TCG -> CGG spells ATCGG.
    gnome::KmerTable kmers{{"ATC", 1}, {"TCG", 1}, {"CGG", 1}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto result = gnome::assemble_graph(graph);

    check(result.error.empty(), "A simple Eulerian path should assemble.");
    check(result.contigs.size() == 1, "A linear graph should produce one contig.");
    if (!result.contigs.empty()) {
        check(result.contigs.front().sequence == "ATCGG", "Linear contig sequence should be ATCGG.");
        check(!result.contigs.front().circular, "Linear contig must not be marked circular.");
        check(result.contigs.front().traversed_edges == 3, "Linear traversal must use all three edges.");
    }
    check(result.traversed_edges == 3, "Result must report all three consumed edges.");
}

void test_circular_contig() {
    // AT -> TC -> CA -> AT. The reported circle contains each edge once;
    // any rotation is biologically equivalent, while deterministic traversal
    // starts at the alphabetically first node and yields ATC.
    gnome::KmerTable kmers{{"ATC", 1}, {"TCA", 1}, {"CAT", 1}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto result = gnome::assemble_graph(graph);

    check(result.error.empty(), "A simple Eulerian cycle should assemble.");
    check(result.contigs.size() == 1, "A circular graph should produce one contig.");
    if (!result.contigs.empty()) {
        check(result.contigs.front().sequence == "ATC", "Circular contig should omit duplicated closing bases.");
        check(result.contigs.front().circular, "Circular contig must be marked circular.");
        check(result.contigs.front().traversed_edges == 3, "Circular traversal must use all three edges.");
    }
}

void test_parallel_edges() {
    // Multiplicity represents parallel edges, each of which must be used.
    gnome::KmerTable kmers{{"ATA", 2}, {"TAT", 2}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto result = gnome::assemble_graph(graph);

    check(result.error.empty(), "A balanced graph with parallel edges should assemble.");
    check(result.traversed_edges == 4, "Traversal must consume every multiplicity unit.");
    if (!result.contigs.empty()) {
        check(result.contigs.front().sequence.size() == 4,
              "Circular contig length should equal the number of traversed edges.");
    }
}

void test_invalid_graph_diagnostic() {
    // AC -> CG and AC -> CT creates degree imbalance greater than one.
    gnome::KmerTable kmers{{"ACG", 1}, {"ACT", 1}};
    const auto graph = gnome::build_graph(kmers, 3);
    const auto result = gnome::assemble_graph(graph);

    check(!result.error.empty(), "A non-Eulerian graph must report an error.");
    check(result.contigs.empty(), "A non-Eulerian graph must not claim a contig.");
    check(result.error.find("not Eulerian") != std::string::npos,
          "Diagnostic should explicitly identify the Eulerian failure.");
}

} // namespace

int main() {
    test_linear_contig();
    test_circular_contig();
    test_parallel_edges();
    test_invalid_graph_diagnostic();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed.\n";
        return 1;
    }
    std::cout << "All Milestone 3 assembly tests passed.\n";
    return 0;
}

#include "gnome_assembler/debruijn_graph.h"
#include "gnome_assembler/kmer_counter.h"
#include "gnome_assembler/kmer_options.h"

#include <iostream>
#include <sstream>
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

// ---------------------------------------------------------------------------
// K-mer counting tests (hand-calculated)
// ---------------------------------------------------------------------------

void test_kmer_extraction() {
    // "ACGT" with k=3 produces: ACG, CGT
    gnome::KmerTable table;
    gnome::count_kmers("ACGT", 3, table);
    check(table.size() == 2, "k=3 on ACGT should produce 2 distinct k-mers.");
    check(table["ACG"] == 1, "ACG count should be 1.");
    check(table["CGT"] == 1, "CGT count should be 1.");
}

void test_kmer_skips_ambiguous() {
    // "ACNGT" with k=3: ACN skipped (N), CNG skipped, NGT skipped — all skipped
    gnome::KmerTable table;
    gnome::count_kmers("ACNGT", 3, table);
    check(table.empty(), "All k-mers containing N should be skipped.");
}

void test_kmer_counts_repeats() {
    // "AAAA" with k=2: AA appears 3 times
    gnome::KmerTable table;
    gnome::count_kmers("AAAA", 2, table);
    check(table.size() == 1, "Only one distinct 2-mer in AAAA.");
    check(table["AA"] == 3, "AA should appear 3 times in AAAA.");
}

void test_kmer_too_short() {
    // sequence shorter than k — nothing should be added
    gnome::KmerTable table;
    gnome::count_kmers("AC", 5, table);
    check(table.empty(), "No k-mers when sequence is shorter than k.");
}

void test_kmer_filtering() {
    gnome::KmerTable table;
    table["AAA"] = 1;
    table["CCC"] = 3;
    table["GGG"] = 2;
    gnome::filter_kmers(table, 2);
    check(table.size() == 2, "Only k-mers with count >= 2 should remain.");
    check(table.count("AAA") == 0, "AAA (count=1) should be filtered out.");
    check(table.count("CCC") == 1, "CCC (count=3) should remain.");
    check(table.count("GGG") == 1, "GGG (count=2) should remain (boundary).");
}

void test_build_kmer_table_statistics() {
    // Two reads that share the k-mer "ACG":
    //   "ACGT" -> ACG, CGT
    //   "ACGC" -> ACG, CGC
    // Raw counts: ACG=2, CGT=1, CGC=1
    // With min_abundance=2: only ACG survives.
    const std::vector<std::string> seqs = {"ACGT", "ACGC"};
    gnome::KmerOptions opts;
    opts.kmer_length = 3;
    opts.min_kmer_abundance = 2;
    gnome::KmerStatistics stats;
    const gnome::KmerTable table = gnome::build_kmer_table(seqs, opts, stats);

    check(stats.total_kmers_seen == 4, "Total k-mers seen should be 4.");
    check(stats.distinct_kmers_before == 3, "Three distinct raw k-mers.");
    check(stats.distinct_kmers_after == 1, "One k-mer survives abundance filter.");
    check(table.size() == 1, "Only ACG should remain.");
    check(table.count("ACG") == 1, "ACG should be in the filtered table.");
}

// ---------------------------------------------------------------------------
// Graph construction tests (hand-calculated)
//
// Hand example 1 (with multiplicity):
//   Reads: "ATCG", "TCGA" with k=3
//   k-mers: ATC, TCG, CGA (TCG appears in both reads -> count=2)
//   Edges:
//     ATC: AT -> TC x1
//     TCG: TC -> CG x2
//     CGA: CG -> GA x1
//   Nodes: AT(in=0,out=1), TC(in=1,out=2), CG(in=2,out=1), GA(in=1,out=0)
// ---------------------------------------------------------------------------

void test_graph_linear() {
    gnome::KmerTable kmers;
    kmers["ATC"] = 1;
    kmers["TCG"] = 2;
    kmers["CGA"] = 1;

    const gnome::DeBruijnGraph graph = gnome::build_graph(kmers, 3);

    check(graph.k == 3, "Graph should record k=3.");
    check(graph.edges.size() == 3, "Three edges expected.");
    check(graph.nodes.size() == 4, "Four nodes: AT, TC, CG, GA.");

    // Check degrees
    check(graph.nodes.at("AT").in_degree == 0, "AT: in=0.");
    check(graph.nodes.at("AT").out_degree == 1, "AT: out=1.");
    check(graph.nodes.at("TC").in_degree == 1, "TC: in=1.");
    check(graph.nodes.at("TC").out_degree == 2, "TC: out=2 (TCG has multiplicity 2).");
    check(graph.nodes.at("CG").in_degree == 2, "CG: in=2.");
    check(graph.nodes.at("CG").out_degree == 1, "CG: out=1.");
    check(graph.nodes.at("GA").in_degree == 1, "GA: in=1.");
    check(graph.nodes.at("GA").out_degree == 0, "GA: out=0.");

    // Edge multiplicities
    const int atc_idx = graph.find_edge("AT", "TC");
    check(atc_idx >= 0, "Edge AT->TC should exist.");
    if (atc_idx >= 0) {
        check(graph.edges[static_cast<std::size_t>(atc_idx)].multiplicity == 1,
              "AT->TC multiplicity should be 1.");
    }
    const int tcg_idx = graph.find_edge("TC", "CG");
    check(tcg_idx >= 0, "Edge TC->CG should exist.");
    if (tcg_idx >= 0) {
        check(graph.edges[static_cast<std::size_t>(tcg_idx)].multiplicity == 2,
              "TC->CG multiplicity should be 2 (k-mer count=2).");
    }

    // Graph analysis
    check(graph.component_count == 1, "One connected component.");
    check(graph.tip_count == 2, "AT and GA are tips.");
    check(graph.branch_count == 1, "TC is a branch (out_degree=2).");
    check(!graph.is_eulerian_cycle, "Linear path with multiplicity is not an Eulerian cycle.");
}

void test_graph_eulerian_path() {
    // Linear path without repeats: ATC, TCG, CGG with k=3
    // Edges:
    //   ATC: AT -> TC x1
    //   TCG: TC -> CG x1
    //   CGG: CG -> GG x1
    // Nodes:
    //   AT: in=0, out=1 (start: out=in+1)
    //   TC: in=1, out=1 (balanced)
    //   CG: in=1, out=1 (balanced)
    //   GG: in=1, out=0 (end: in=out+1)
    gnome::KmerTable kmers;
    kmers["ATC"] = 1;
    kmers["TCG"] = 1;
    kmers["CGG"] = 1;

    const gnome::DeBruijnGraph graph = gnome::build_graph(kmers, 3);
    check(graph.nodes.size() == 4, "Four nodes in linear path.");
    check(graph.edges.size() == 3, "Three edges in linear path.");
    check(graph.is_eulerian_path, "Simple linear path must be an Eulerian path.");
    check(!graph.is_eulerian_cycle, "Linear path is not an Eulerian cycle.");
    check(graph.tip_count == 2, "Start and end nodes are tips.");
    check(graph.branch_count == 0, "No branching in simple path.");
}

void test_graph_circular() {
    // Circular graph: ATC, TCA, CAT with k=3 (each appearing once).
    // Edges:
    //   ATC: AT -> TC
    //   TCA: TC -> CA
    //   CAT: CA -> AT
    // All nodes: in=1, out=1 -> Eulerian cycle.
    gnome::KmerTable kmers;
    kmers["ATC"] = 1;
    kmers["TCA"] = 1;
    kmers["CAT"] = 1;

    const gnome::DeBruijnGraph graph = gnome::build_graph(kmers, 3);

    check(graph.nodes.size() == 3, "Three nodes in circular graph.");
    check(graph.edges.size() == 3, "Three edges.");
    for (const auto& [name, degree] : graph.nodes) {
        check(degree.in_degree == 1 && degree.out_degree == 1,
              "Every node in a triangle should have in=1, out=1.");
    }
    check(graph.is_eulerian_cycle, "Triangle should be an Eulerian cycle.");
    check(!graph.is_eulerian_path, "Eulerian cycle should not also be flagged as a path.");
    check(graph.tip_count == 0, "No tips in circular graph.");
    check(graph.branch_count == 0, "No branches in simple circle.");
}

void test_graph_dump_deterministic() {
    // Build a tiny graph and check that dump output is sorted.
    gnome::KmerTable kmers;
    kmers["CGT"] = 1;
    kmers["ACG"] = 2;
    const gnome::DeBruijnGraph graph = gnome::build_graph(kmers, 3);

    std::ostringstream out1, out2;
    gnome::dump_graph(graph, out1);
    gnome::dump_graph(graph, out2);
    check(out1.str() == out2.str(), "Graph dump must be deterministic (identical on two calls).");
    check(!out1.str().empty(), "Graph dump must not be empty.");
    // Check that "AC -> CG" appears in the dump (from ACG).
    check(out1.str().find("AC -> CG") != std::string::npos,
          "Dump should contain edge AC -> CG.");
}

void test_node_kind_classification() {
    // Balanced
    check(gnome::classify_node({2, 2}) == gnome::NodeKind::Normal,
          "Balanced node should be Normal.");
    check(gnome::classify_node({1, 1}) == gnome::NodeKind::Normal,
          "in=1, out=1 should be Normal.");

    // Tip (dead end: in=0 or out=0)
    check(gnome::classify_node({0, 1}) == gnome::NodeKind::Tip,
          "in=0 node should be a Tip.");
    check(gnome::classify_node({1, 0}) == gnome::NodeKind::Tip,
          "out=0 node should be a Tip.");

    // Start node: out = in + 1 (with in > 0)
    check(gnome::classify_node({1, 2}) == gnome::NodeKind::StartNode,
          "out=in+1 with in>0 should be a StartNode.");

    // End node: in = out + 1 (with out > 0)
    check(gnome::classify_node({2, 1}) == gnome::NodeKind::EndNode,
          "in=out+1 with out>0 should be an EndNode.");

    // Branch: out > 1 (when not balanced and not simple start/end)
    check(gnome::classify_node({1, 3}) == gnome::NodeKind::Branch,
          "out=3, in=1 should be a Branch.");
}

} // namespace

int main() {
    test_kmer_extraction();
    test_kmer_skips_ambiguous();
    test_kmer_counts_repeats();
    test_kmer_too_short();
    test_kmer_filtering();
    test_build_kmer_table_statistics();
    test_graph_linear();
    test_graph_eulerian_path();
    test_graph_circular();
    test_graph_dump_deterministic();
    test_node_kind_classification();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed.\n";
        return 1;
    }
    std::cout << "All Milestone 2 k-mer and graph tests passed.\n";
    return 0;
}

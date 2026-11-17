#include "gnome_assembler/graph_pipeline.h"

#include "gnome_assembler/cleaner.h"
#include "gnome_assembler/fastq_reader.h"
#include "gnome_assembler/output_writer.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace gnome {

namespace {

std::filesystem::path resolve_graph_prefix(const std::string& requested,
                                           std::error_code& ec) {
    std::filesystem::path prefix = requested;
    if (std::filesystem::is_directory(prefix, ec)) {
        if (!ec) {
            prefix /= "gnome_graph";
        }
    } else if (ec == std::errc::no_such_file_or_directory) {
        ec.clear();
    }
    return prefix;
}

} // namespace

std::string format_graph_summary(const KmerStatistics& kmer_stats,
                                 const DeBruijnGraph& graph,
                                 const KmerOptions& kmer_options) {
    std::ostringstream out;
    out << "GnomeAssembler graph summary\n"
        << "k-mer length: "      << kmer_options.kmer_length       << "\n"
        << "Min abundance: "     << kmer_options.min_kmer_abundance << "\n"
        << "Total k-mers seen: " << kmer_stats.total_kmers_seen     << "\n"
        << "Distinct (raw): "    << kmer_stats.distinct_kmers_before << "\n"
        << "Distinct (kept): "   << kmer_stats.distinct_kmers_after  << "\n"
        << "Nodes: "             << graph.nodes.size()               << "\n"
        << "Edges: "             << graph.edges.size()               << "\n"
        << "Components: "        << graph.component_count            << "\n"
        << "Tips: "              << graph.tip_count                  << "\n"
        << "Branches: "          << graph.branch_count               << "\n"
        << "Eulerian cycle: "    << (graph.is_eulerian_cycle ? "yes" : "no") << "\n"
        << "Eulerian path: "     << (graph.is_eulerian_path  ? "yes" : "no") << "\n";
    return out.str();
}

GraphPipelineResult run_graph_pipeline(const GraphPipelineOptions& options) {
    GraphPipelineResult result;

    // ----------------------------------------------------------------
    // 1. Validate basic inputs.
    // ----------------------------------------------------------------
    if (options.cleaning.reads_path.empty()) {
        result.error = "An input FASTQ file is required.";
        return result;
    }
    if (options.kmer.kmer_length < 2) {
        result.error = "k-mer length must be at least 2.";
        return result;
    }
    if (options.kmer.min_kmer_abundance < 1) {
        result.error = "Minimum k-mer abundance must be at least 1.";
        return result;
    }

    // ----------------------------------------------------------------
    // 2. Stream the FASTQ, clean each read, collect retained sequences.
    // ----------------------------------------------------------------
    std::ifstream input(options.cleaning.reads_path);
    if (!input) {
        result.error = "Could not open FASTQ input: " + options.cleaning.reads_path;
        return result;
    }

    std::vector<std::string> sequences;
    try {
        FastqReader reader(input);
        FastqRecord record;
        while (reader.next(record)) {
            const CleaningResult cleaned = clean_record(record, options.cleaning);
            if (cleaned.retained) {
                sequences.push_back(cleaned.record.sequence);
            }
        }
    } catch (const std::exception& e) {
        result.error = std::string("FASTQ parse error: ") + e.what();
        return result;
    }

    if (sequences.empty()) {
        result.error = "No reads survived cleaning. Try relaxing the quality or length thresholds.";
        return result;
    }
    result.retained_reads = sequences;

    // ----------------------------------------------------------------
    // 3. Build k-mer table and graph.
    // ----------------------------------------------------------------
    KmerStatistics kmer_stats;
    const KmerTable kmer_table = build_kmer_table(sequences, options.kmer, kmer_stats);
    result.kmer_stats = kmer_stats;

    if (kmer_table.empty()) {
        result.error = "No k-mers survived abundance filtering. "
                       "Try a lower --min-kmer-abundance or a shorter k.";
        return result;
    }

    result.graph = build_graph(kmer_table, options.kmer.kmer_length);

    // ----------------------------------------------------------------
    // 4. Format summary.
    // ----------------------------------------------------------------
    result.summary = format_graph_summary(kmer_stats, result.graph, options.kmer);

    // ----------------------------------------------------------------
    // 5. Optionally write the graph dump to disk.
    // ----------------------------------------------------------------
    if (!options.output_prefix.empty()) {
        std::error_code ec;
        const std::filesystem::path prefix = resolve_graph_prefix(options.output_prefix, ec);
        if (!ec) {
            const std::filesystem::path parent = prefix.parent_path();
            if (!parent.empty()) {
                std::filesystem::create_directories(parent, ec);
            }
        }
        if (ec) {
            result.error = "Cannot write graph dump: " + ec.message();
            return result;
        }
        const std::string dump_path = prefix.string() + ".graph-dump.txt";
        std::ofstream dump_file(dump_path);
        if (!dump_file) {
            result.error = "Could not create graph dump file: " + dump_path;
            return result;
        }
        dump_graph(result.graph, dump_file);
        result.summary += "Graph dump: " + dump_path + "\n";
    }

    return result;
}

} // namespace gnome

#include "gnome_assembler/assembly_pipeline.h"

#include <filesystem>
#include <fstream>
#include <sstream>

namespace gnome {

std::string format_assembly_summary(const GraphPipelineResult& graph_result,
                                    const AssemblyResult& assembly) {
    std::ostringstream out;
    out << graph_result.summary;
    if (!assembly.error.empty()) return out.str();

    out << "Contigs: " << assembly.contigs.size() << "\n";
    out << "Edges traversed: " << assembly.traversed_edges << "\n";
    for (std::size_t i = 0; i < assembly.contigs.size(); ++i) {
        const AssemblyContig& contig = assembly.contigs[i];
        out << "Contig " << (i + 1) << ": " << contig.sequence.size()
            << " bp (" << (contig.circular ? "circular" : "linear") << ")\n";
    }
    for (const std::string& warning : assembly.warnings) {
        out << "Warning: " << warning << "\n";
    }
    return out.str();
}

AssemblyPipelineResult run_assembly_pipeline(const GraphPipelineOptions& options) {
    AssemblyPipelineResult result;
    result.graph_result = run_graph_pipeline(options);
    if (!result.graph_result.error.empty()) {
        result.error = result.graph_result.error;
        result.summary = result.graph_result.summary;
        return result;
    }

    result.assembly = assemble_graph(result.graph_result.graph);
    result.summary = format_assembly_summary(result.graph_result, result.assembly);
    result.verification = verify_assembly(result.graph_result.graph, result.assembly,
                                          result.graph_result.retained_reads);
    result.summary += "\n" + format_verification_text(result.verification);

    if (!options.output_prefix.empty()) {
        std::filesystem::path prefix(options.output_prefix);
        std::error_code ec;
        if (std::filesystem::is_directory(prefix, ec) && !ec) prefix /= "gnome_graph";
        const std::string text_path = prefix.string() + ".assembly-report.txt";
        const std::string json_path = prefix.string() + ".assembly-report.json";
        std::ofstream text_output(text_path);
        std::ofstream json_output(json_path);
        if (!text_output || !json_output) {
            result.error = "Could not create assembly report files beside: " + prefix.string();
            return result;
        }
        text_output << format_verification_text(result.verification);
        json_output << format_verification_json(result.verification);
        result.summary += "Text report: " + text_path + "\n";
        result.summary += "JSON report: " + json_path + "\n";
    }
    if (!result.assembly.error.empty()) {
        result.error = result.assembly.error;
        return result;
    }

    if (!options.output_prefix.empty()) {
        std::filesystem::path prefix(options.output_prefix);
        std::error_code ec;
        if (std::filesystem::is_directory(prefix, ec) && !ec) prefix /= "gnome_graph";
        const std::string path = prefix.string() + ".contigs.fasta";
        std::ofstream output(path);
        if (!output) {
            result.error = "Could not create contig FASTA file: " + path;
            return result;
        }
        for (std::size_t i = 0; i < result.assembly.contigs.size(); ++i) {
            const AssemblyContig& contig = result.assembly.contigs[i];
            output << ">contig_" << (i + 1) << " "
                   << (contig.circular ? "circular" : "linear")
                   << " edges=" << contig.traversed_edges << "\n"
                   << contig.sequence << "\n";
        }
        result.summary += "Contigs FASTA: " + path + "\n";
    }
    return result;
}

} // namespace gnome

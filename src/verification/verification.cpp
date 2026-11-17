#include "gnome_assembler/verification.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

namespace gnome {

namespace {

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return base;
    }
}

std::string reverse_complement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        result.push_back(complement(*it));
    }
    return result;
}

bool contig_contains(const AssemblyContig& contig, const std::string& read) {
    if (read.empty() || contig.sequence.empty()) return false;
    if (!contig.circular && read.size() > contig.sequence.size()) return false;
    std::string searchable = contig.sequence;
    if (contig.circular && read.size() > 1) {
        const std::size_t extension_length = read.size() - 1;
        while (searchable.size() < contig.sequence.size() + extension_length) {
            const std::size_t remaining = contig.sequence.size() + extension_length - searchable.size();
            searchable += contig.sequence.substr(0, remaining);
        }
    }
    return searchable.find(read) != std::string::npos;
}

bool read_is_covered(const std::vector<AssemblyContig>& contigs, const std::string& read) {
    const std::string reverse = reverse_complement(read);
    for (const AssemblyContig& contig : contigs) {
        if (contig_contains(contig, read) || contig_contains(contig, reverse)) return true;
    }
    return false;
}

std::string json_escape(const std::string& value) {
    std::string escaped;
    for (char character : value) {
        switch (character) {
            case '\\': escaped += "\\\\"; break;
            case '"': escaped += "\\\""; break;
            case '\n': escaped += "\\n"; break;
            case '\r': escaped += "\\r"; break;
            case '\t': escaped += "\\t"; break;
            default: escaped += character; break;
        }
    }
    return escaped;
}

} // namespace

VerificationReport verify_assembly(const DeBruijnGraph& graph,
                                   const AssemblyResult& assembly,
                                   const std::vector<std::string>& cleaned_reads) {
    VerificationReport report;
    report.contig_count = assembly.contigs.size();
    report.cleaned_reads = cleaned_reads.size();
    report.traversed_edges = assembly.traversed_edges;
    for (const GraphEdge& edge : graph.edges) report.graph_edges += edge.multiplicity;

    std::vector<std::size_t> lengths;
    lengths.reserve(assembly.contigs.size());
    for (const AssemblyContig& contig : assembly.contigs) {
        const std::size_t length = contig.sequence.size();
        report.total_contig_length += length;
        report.longest_contig_length = std::max(report.longest_contig_length, length);
        lengths.push_back(length);
    }
    std::sort(lengths.begin(), lengths.end(), std::greater<std::size_t>());
    std::size_t accumulated = 0;
    for (std::size_t length : lengths) {
        accumulated += length;
        if (accumulated * 2 >= report.total_contig_length) {
            report.n50 = length;
            break;
        }
    }

    for (const std::string& read : cleaned_reads) {
        if (read_is_covered(assembly.contigs, read)) ++report.covered_reads;
    }
    report.all_edges_used = !graph.edges.empty() &&
                            report.graph_edges == report.traversed_edges;
    report.assembly_complete = assembly.error.empty() && !assembly.contigs.empty() &&
                               report.all_edges_used &&
                               report.covered_reads == report.cleaned_reads;

    if (!assembly.error.empty()) report.warnings.push_back("Assembly did not complete: " + assembly.error);
    if (!report.all_edges_used) report.warnings.push_back("Not every retained graph edge was traversed.");
    if (report.covered_reads != report.cleaned_reads) {
        report.warnings.push_back("Some cleaned reads are not represented by the reported contigs.");
    }
    for (const std::string& warning : assembly.warnings) report.warnings.push_back(warning);
    if (report.assembly_complete && !assembly.warnings.empty()) {
        report.warnings.push_back("All retained edges and reads are represented, but graph ambiguity prevents claiming a unique assembly.");
    }
    return report;
}

std::string format_verification_text(const VerificationReport& report) {
    std::ostringstream out;
    const double coverage = report.cleaned_reads == 0 ? 0.0 :
        100.0 * static_cast<double>(report.covered_reads) / static_cast<double>(report.cleaned_reads);
    out << "GnomeAssembler verification report\n"
        << "Assembly complete: " << (report.assembly_complete ? "yes" : "no") << "\n"
        << "All retained edges used: " << (report.all_edges_used ? "yes" : "no") << "\n"
        << "Graph edges (with multiplicity): " << report.graph_edges << "\n"
        << "Edges traversed: " << report.traversed_edges << "\n"
        << "Contigs: " << report.contig_count << "\n"
        << "Total contig length: " << report.total_contig_length << " bp\n"
        << "Longest contig: " << report.longest_contig_length << " bp\n"
        << "N50: " << report.n50 << " bp\n"
        << "Cleaned reads: " << report.cleaned_reads << "\n"
        << "Covered reads: " << report.covered_reads << " (" << std::fixed << std::setprecision(2)
        << coverage << "%)\n";
    for (const std::string& warning : report.warnings) out << "Warning: " << warning << "\n";
    return out.str();
}

std::string format_verification_json(const VerificationReport& report) {
    std::ostringstream out;
    const double coverage = report.cleaned_reads == 0 ? 0.0 :
        100.0 * static_cast<double>(report.covered_reads) / static_cast<double>(report.cleaned_reads);
    out << "{\n"
        << "  \"assembly_complete\": " << (report.assembly_complete ? "true" : "false") << ",\n"
        << "  \"all_retained_edges_used\": " << (report.all_edges_used ? "true" : "false") << ",\n"
        << "  \"graph_edges\": " << report.graph_edges << ",\n"
        << "  \"traversed_edges\": " << report.traversed_edges << ",\n"
        << "  \"contig_count\": " << report.contig_count << ",\n"
        << "  \"total_contig_length\": " << report.total_contig_length << ",\n"
        << "  \"longest_contig_length\": " << report.longest_contig_length << ",\n"
        << "  \"n50\": " << report.n50 << ",\n"
        << "  \"cleaned_reads\": " << report.cleaned_reads << ",\n"
        << "  \"covered_reads\": " << report.covered_reads << ",\n"
        << "  \"covered_read_percent\": " << std::fixed << std::setprecision(2) << coverage << ",\n"
        << "  \"warnings\": [";
    for (std::size_t i = 0; i < report.warnings.size(); ++i) {
        if (i != 0) out << ", ";
        out << "\"" << json_escape(report.warnings[i]) << "\"";
    }
    out << "]\n}\n";
    return out.str();
}

} // namespace gnome

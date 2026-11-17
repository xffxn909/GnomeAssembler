#include "gnome_assembler/blast_identification.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#ifdef _WIN32
#include <process.h>
#endif

namespace gnome {

namespace {

bool database_exists(const std::string& prefix) {
    namespace fs = std::filesystem;
    if (prefix.empty()) return false;
    std::error_code ec;
    if (fs::is_regular_file(prefix, ec)) return true;
    for (const char* suffix : {".ndb", ".nin", ".nhr", ".nsq"}) {
        if (fs::exists(prefix + suffix, ec)) return true;
        ec.clear();
    }
    return false;
}

bool better_hit(const BlastHit& left, const BlastHit& right) {
    if (left.bit_score != right.bit_score) return left.bit_score > right.bit_score;
    if (left.evalue != right.evalue) return left.evalue < right.evalue;
    if (left.identity_percent != right.identity_percent) return left.identity_percent > right.identity_percent;
    if (left.query_coverage_percent != right.query_coverage_percent) {
        return left.query_coverage_percent > right.query_coverage_percent;
    }
    if (left.subject_id != right.subject_id) return left.subject_id < right.subject_id;
    return left.query_id < right.query_id;
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

void write_hit_json(std::ostream& out, const BlastHit& hit) {
    out << "{\"query_id\": \"" << json_escape(hit.query_id)
        << "\", \"subject_id\": \"" << json_escape(hit.subject_id)
        << "\", \"title\": \"" << json_escape(hit.title)
        << "\", \"identity_percent\": " << hit.identity_percent
        << ", \"query_coverage_percent\": " << hit.query_coverage_percent
        << ", \"evalue\": " << hit.evalue
        << ", \"bit_score\": " << hit.bit_score << "}";
}

} // namespace

bool blast_options_requested(const BlastOptions& options) {
    return !options.executable_path.empty() || !options.database_path.empty();
}

std::vector<BlastHit> parse_blast_tabular(const std::string& tabular_text, std::string& error) {
    std::vector<BlastHit> hits;
    std::istringstream input(tabular_text);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const std::size_t end = line.find('\t', start);
            fields.push_back(line.substr(start, end == std::string::npos ? std::string::npos : end - start));
            if (end == std::string::npos) break;
            start = end + 1;
        }
        if (fields.size() < 9) {
            error = "Malformed BLAST tabular result on line " + std::to_string(line_number) + ".";
            return {};
        }
        try {
            BlastHit hit;
            hit.query_id = fields[0];
            hit.subject_id = fields[1];
            hit.identity_percent = std::stod(fields[2]);
            hit.alignment_length = static_cast<std::size_t>(std::stoull(fields[3]));
            hit.query_length = static_cast<std::size_t>(std::stoull(fields[4]));
            hit.evalue = std::stod(fields[6]);
            hit.bit_score = std::stod(fields[7]);
            hit.title = fields[8];
            if (hit.query_length == 0) {
                error = "BLAST result has a zero query length on line " + std::to_string(line_number) + ".";
                return {};
            }
            hit.query_coverage_percent = 100.0 * static_cast<double>(hit.alignment_length) /
                                         static_cast<double>(hit.query_length);
            hits.push_back(std::move(hit));
        } catch (const std::exception&) {
            error = "Malformed numeric BLAST value on line " + std::to_string(line_number) + ".";
            return {};
        }
    }
    return hits;
}

std::vector<BlastHit> filter_and_rank_blast_hits(const std::vector<BlastHit>& hits,
                                                 const BlastOptions& options) {
    std::vector<BlastHit> qualified;
    for (const BlastHit& hit : hits) {
        if (hit.identity_percent >= options.minimum_identity &&
            hit.query_coverage_percent >= options.minimum_query_coverage) {
            qualified.push_back(hit);
        }
    }
    std::sort(qualified.begin(), qualified.end(), better_hit);
    return qualified;
}

BlastResult identify_contigs_with_blast(const std::string& contigs_fasta,
                                        const BlastOptions& options,
                                        const std::string& tabular_output_path) {
    BlastResult result;
    result.requested = blast_options_requested(options);
    if (!result.requested) {
        result.status = "Identification not requested.";
        return result;
    }
    if (options.executable_path.empty() || options.database_path.empty()) {
        result.error = "BLAST identification requires both a blastn executable and a database prefix.";
        result.status = "Not identified: BLAST settings are incomplete.";
        return result;
    }
    if (!std::filesystem::is_regular_file(options.executable_path)) {
        result.error = "Could not find blastn executable: " + options.executable_path;
        result.status = "Not identified: blastn executable unavailable.";
        return result;
    }
    if (!database_exists(options.database_path)) {
        result.error = "Could not find a formatted BLAST database at prefix: " + options.database_path;
        result.status = "Not identified: formatted database unavailable.";
        return result;
    }
    if (!std::filesystem::is_regular_file(contigs_fasta)) {
        result.error = "Contig FASTA was not created: " + contigs_fasta;
        result.status = "Not identified: contigs unavailable.";
        return result;
    }

    const std::vector<std::string> arguments{
        options.executable_path, "-query", contigs_fasta, "-db", options.database_path,
        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore stitle",
        "-out", tabular_output_path
    };
#ifdef _WIN32
    std::vector<const char*> argv;
    argv.reserve(arguments.size() + 1);
    for (const std::string& argument : arguments) argv.push_back(argument.c_str());
    argv.push_back(nullptr);
    const int exit_code = _spawnv(_P_WAIT, options.executable_path.c_str(), argv.data());
    if (exit_code != 0) {
        result.error = "blastn failed with exit code " + std::to_string(exit_code) + ".";
        result.status = "Not identified: blastn execution failed.";
        return result;
    }
#else
    result.error = "Local BLAST execution is currently supported on Windows builds only.";
    result.status = "Not identified: BLAST execution unavailable on this platform.";
    return result;
#endif

    std::ifstream output(tabular_output_path);
    if (!output) {
        result.error = "blastn did not create the expected tabular output: " + tabular_output_path;
        result.status = "Not identified: BLAST output unavailable.";
        return result;
    }
    std::ostringstream text;
    text << output.rdbuf();
    std::string parse_error;
    const std::vector<BlastHit> all_hits = parse_blast_tabular(text.str(), parse_error);
    if (!parse_error.empty()) {
        result.error = parse_error;
        result.status = "Not identified: BLAST output could not be parsed.";
        return result;
    }
    result.qualified_hits = filter_and_rank_blast_hits(all_hits, options);
    if (result.qualified_hits.empty()) {
        result.status = "Not identified: no hit met the configured identity and coverage thresholds.";
        return result;
    }
    result.identified = true;
    result.best_hit = result.qualified_hits.front();
    result.status = "Identified: " + result.best_hit.subject_id + ".";
    return result;
}

std::string format_blast_report_text(const BlastResult& result, const BlastOptions& options) {
    std::ostringstream out;
    out << "GnomeAssembler BLAST+ identification report\n"
        << "Status: " << result.status << "\n"
        << "blastn: " << (options.executable_path.empty() ? "not configured" : options.executable_path) << "\n"
        << "Database: " << (options.database_path.empty() ? "not configured" : options.database_path) << "\n"
        << "Minimum identity: " << options.minimum_identity << "%\n"
        << "Minimum query coverage: " << options.minimum_query_coverage << "%\n";
    if (result.identified) {
        const BlastHit& hit = result.best_hit;
        out << "Best subject: " << hit.subject_id << "\n"
            << "Title: " << hit.title << "\n"
            << "Query: " << hit.query_id << "\n"
            << "Identity: " << hit.identity_percent << "%\n"
            << "Query coverage: " << hit.query_coverage_percent << "%\n"
            << "E-value: " << hit.evalue << "\n"
            << "Bit score: " << hit.bit_score << "\n";
    }
    if (!result.error.empty()) out << "Detail: " << result.error << "\n";
    return out.str();
}

std::string format_blast_report_json(const BlastResult& result, const BlastOptions& options) {
    std::ostringstream out;
    out << "{\n  \"requested\": " << (result.requested ? "true" : "false")
        << ",\n  \"identified\": " << (result.identified ? "true" : "false")
        << ",\n  \"status\": \"" << json_escape(result.status) << "\""
        << ",\n  \"error\": \"" << json_escape(result.error) << "\""
        << ",\n  \"blastn\": \"" << json_escape(options.executable_path) << "\""
        << ",\n  \"database\": \"" << json_escape(options.database_path) << "\""
        << ",\n  \"minimum_identity\": " << options.minimum_identity
        << ",\n  \"minimum_query_coverage\": " << options.minimum_query_coverage
        << ",\n  \"best_hit\": ";
    if (result.identified) write_hit_json(out, result.best_hit); else out << "null";
    out << "\n}\n";
    return out.str();
}

} // namespace gnome

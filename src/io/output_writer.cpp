#include "gnome_assembler/output_writer.h"

#include "gnome_assembler/fastq_reader.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace gnome {

namespace {

std::filesystem::path resolve_output_prefix(const std::filesystem::path& requested_path,
                                            std::error_code& error) {
    std::filesystem::path prefix = requested_path;
    if (std::filesystem::is_directory(prefix, error)) {
        if (error) {
            return {};
        }
        prefix /= "cleaned_reads";
    }
    if (error == std::errc::no_such_file_or_directory) {
        error.clear();
    }
    return prefix;
}

std::string format_directory_error(const std::filesystem::path& directory,
                                   const std::error_code& error) {
    return "Cannot use output folder '" + directory.string() + "': " + error.message() +
           ". Choose a folder you can write to, such as a folder under Documents.";
}

} // namespace

std::string format_cleaning_summary(const CleaningStatistics& statistics,
                                    const CleaningOptions& options) {
    std::ostringstream summary;
    summary << "GnomeAssembler cleaning summary\n"
            << "Input reads: " << statistics.input_reads << "\n"
            << "Retained reads: " << statistics.retained_reads << "\n"
            << "Trimmed reads: " << statistics.trimmed_reads << "\n"
            << "Discarded reads: " << statistics.discarded_reads << "\n"
            << "Input bases: " << statistics.input_bases << "\n"
            << "Retained bases: " << statistics.retained_bases << "\n"
            << "Minimum quality: " << options.minimum_quality << "\n"
            << "Minimum read length: " << options.minimum_read_length << "\n";
    return summary.str();
}

PipelineResult run_cleaning_to_files(const CleaningOptions& options) {
    PipelineResult result;
    if (options.reads_path.empty() || options.output_prefix.empty()) {
        result.error = "An input FASTQ file and output prefix are required.";
        return result;
    }
    if (options.minimum_quality > 93) {
        result.error = "Minimum quality must be from 0 to 93 for Phred+33 input.";
        return result;
    }

    std::error_code path_error;
    const std::filesystem::path prefix = resolve_output_prefix(options.output_prefix, path_error);
    if (path_error) {
        result.error = format_directory_error(options.output_prefix, path_error);
        return result;
    }
    const std::filesystem::path parent = prefix.parent_path();
    std::error_code directory_error;
    if (!parent.empty()) {
        std::filesystem::create_directories(parent, directory_error);
        if (directory_error) {
            result.error = format_directory_error(parent, directory_error);
            return result;
        }
    }

    std::ifstream input(options.reads_path);
    if (!input) {
        result.error = "Could not open FASTQ input: " + options.reads_path;
        return result;
    }

    const auto fastq_path = prefix.string() + ".cleaned.fastq";
    const auto fasta_path = prefix.string() + ".cleaned.fasta";
    const auto summary_path = prefix.string() + ".cleaning-summary.txt";
    std::ofstream cleaned_fastq(fastq_path);
    std::ofstream cleaned_fasta(fasta_path);
    if (!cleaned_fastq || !cleaned_fasta) {
        result.error = "Could not create cleaned output files with prefix: " + prefix.string();
        return result;
    }

    CleaningStatistics statistics;
    try {
        FastqReader reader(input);
        FastqRecord record;
        while (reader.next(record)) {
            const CleaningResult cleaned = clean_record(record, options);
            update_statistics(statistics, cleaned);
            if (cleaned.retained) {
                cleaned_fastq << '@' << cleaned.record.name << '\n' << cleaned.record.sequence << "\n+\n"
                              << cleaned.record.quality << '\n';
                cleaned_fasta << '>' << cleaned.record.name << '\n' << cleaned.record.sequence << '\n';
            }
        }
    } catch (const std::exception& error) {
        result.error = error.what();
        result.summary = format_cleaning_summary(statistics, options);
        return result;
    }

    result.summary = format_cleaning_summary(statistics, options);
    result.summary += "Cleaned FASTQ: " + fastq_path + "\n";
    result.summary += "Cleaned FASTA: " + fasta_path + "\n";
    result.summary += "Summary file: " + summary_path + "\n";
    std::ofstream summary_file(summary_path);
    if (!summary_file) {
        result.error = "Could not create summary file: " + summary_path;
        return result;
    }
    summary_file << result.summary;
    return result;
}

} // namespace gnome

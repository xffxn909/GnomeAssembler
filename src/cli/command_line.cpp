#include "gnome_assembler/command_line.h"

#include <charconv>
#include <string>

namespace gnome {

namespace {

bool parse_nonnegative_integer(const std::string& text, int& value) {
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size() && value >= 0;
}

bool is_cleaning_option(const std::string& arg) {
    return arg == "--reads" || arg == "--output" ||
           arg == "--min-quality" || arg == "--min-read-length" ||
           arg == "--adapter-prefix" || arg == "--adapter-suffix";
}

bool is_graph_option(const std::string& arg) {
    return is_cleaning_option(arg) ||
           arg == "--kmer-length" || arg == "--min-kmer-abundance";
}

bool apply_cleaning_option(const std::string& arg, const std::string& value,
                            CleaningOptions& opts, std::string& error) {
    if (arg == "--reads")          { opts.reads_path = value; return true; }
    if (arg == "--output")         { opts.output_prefix = value; return true; }
    if (arg == "--adapter-prefix") { opts.adapter_prefix = value; return true; }
    if (arg == "--adapter-suffix") { opts.adapter_suffix = value; return true; }
    int number = 0;
    if (!parse_nonnegative_integer(value, number)) {
        error = "Option '" + arg + "' needs a non-negative integer.";
        return false;
    }
    if (arg == "--min-quality")      opts.minimum_quality = number;
    else if (arg == "--min-read-length") opts.minimum_read_length = static_cast<std::size_t>(number);
    return true;
}

} // namespace

ParsedCommand parse_command_line(int argc, char* argv[]) {
    ParsedCommand command;
    if (argc == 1) {
        return command;
    }

    int index = 1;
    const std::string first = argv[index];
    if (first == "--help" || first == "-h") {
        command.show_help = true;
        return command;
    }

    if (first == "clean") {
        command.has_cleaning_request = true;
        ++index;
        while (index < argc) {
            const std::string arg = argv[index++];
            if (!is_cleaning_option(arg)) {
                command.error = "Unknown option '" + arg + "'.";
                return command;
            }
            if (index == argc) {
                command.error = "Option '" + arg + "' needs a value.";
                return command;
            }
            const std::string value = argv[index++];
            if (!apply_cleaning_option(arg, value, command.options, command.error)) {
                return command;
            }
        }
        if (command.options.reads_path.empty() || command.options.output_prefix.empty()) {
            command.error = "Both --reads and --output are required for 'clean'.";
        }
        return command;
    }

    if (first == "build-graph") {
        command.has_graph_request = true;
        ++index;
        while (index < argc) {
            const std::string arg = argv[index++];
            if (!is_graph_option(arg)) {
                command.error = "Unknown option '" + arg + "'.";
                return command;
            }
            if (index == argc) {
                command.error = "Option '" + arg + "' needs a value.";
                return command;
            }
            const std::string value = argv[index++];
            if (arg == "--kmer-length" || arg == "--min-kmer-abundance") {
                int number = 0;
                if (!parse_nonnegative_integer(value, number) || number < 1) {
                    command.error = "Option '" + arg + "' needs a positive integer.";
                    return command;
                }
                if (arg == "--kmer-length")
                    command.graph_options.kmer.kmer_length = number;
                else
                    command.graph_options.kmer.min_kmer_abundance = static_cast<uint32_t>(number);
            } else {
                if (!apply_cleaning_option(arg, value,
                                           command.graph_options.cleaning, command.error)) {
                    return command;
                }
                // Also mirror reads/output into output_prefix for the graph dump.
                if (arg == "--reads")
                    command.graph_options.output_prefix =
                        command.graph_options.cleaning.output_prefix;
                if (arg == "--output")
                    command.graph_options.output_prefix = value;
            }
        }
        if (command.graph_options.cleaning.reads_path.empty()) {
            command.error = "--reads is required for 'build-graph'.";
        }
        return command;
    }

    command.error = "Expected 'clean', 'build-graph', or no command to open the GUI.";
    return command;
}

std::string command_help() {
    return
        "GnomeAssembler - de novo virus genome assembler\n\n"
        "Open the GUI:\n"
        "  gnome_assembler\n\n"
        "Milestone 1 - clean FASTQ reads:\n"
        "  gnome_assembler clean --reads INPUT.fastq --output PREFIX [options]\n\n"
        "Milestone 2 - build de Bruijn graph:\n"
        "  gnome_assembler build-graph --reads CLEANED.fastq --output PREFIX [options]\n\n"
        "Shared cleaning options:\n"
        "  --min-quality NUMBER       Minimum Phred+33 quality at the 3' end (default: 20)\n"
        "  --min-read-length NUMBER   Discard shorter reads after trimming (default: 50)\n"
        "  --adapter-prefix DNA       Remove this exact adapter from the 5' end\n"
        "  --adapter-suffix DNA       Remove this exact adapter from the 3' end\n\n"
        "Graph options:\n"
        "  --kmer-length NUMBER       k-mer length (default: 21)\n"
        "  --min-kmer-abundance NUMBER Discard k-mers below this count (default: 2)\n";
}

} // namespace gnome

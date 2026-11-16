#include "gnome_assembler/command_line.h"

#include <charconv>

namespace gnome {

namespace {

bool parse_nonnegative_integer(const std::string& text, int& value) {
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size() && value >= 0;
}

bool needs_value(const std::string& argument) {
    return argument == "--reads" || argument == "--output" ||
           argument == "--min-quality" || argument == "--min-read-length" ||
           argument == "--adapter-prefix" || argument == "--adapter-suffix";
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
    if (first != "clean" && first != "assemble") {
        command.error = "Expected 'clean', 'assemble', or no command to open the GUI.";
        return command;
    }
    command.has_cleaning_request = true;
    ++index;

    while (index < argc) {
        const std::string argument = argv[index++];
        if (!needs_value(argument)) {
            command.error = "Unknown option '" + argument + "'.";
            return command;
        }
        if (index == argc) {
            command.error = "Option '" + argument + "' needs a value.";
            return command;
        }
        const std::string value = argv[index++];
        if (argument == "--reads") command.options.reads_path = value;
        else if (argument == "--output") command.options.output_prefix = value;
        else if (argument == "--adapter-prefix") command.options.adapter_prefix = value;
        else if (argument == "--adapter-suffix") command.options.adapter_suffix = value;
        else {
            int number = 0;
            if (!parse_nonnegative_integer(value, number)) {
                command.error = "Option '" + argument + "' needs a non-negative integer.";
                return command;
            }
            if (argument == "--min-quality") command.options.minimum_quality = number;
            else command.options.minimum_read_length = static_cast<std::size_t>(number);
        }
    }

    if (command.options.reads_path.empty() || command.options.output_prefix.empty()) {
        command.error = "Both --reads and --output are required.";
    }
    return command;
}

std::string command_help() {
    return
        "GnomeAssembler Milestone 1 - FASTQ cleaning\n\n"
        "Open the GUI:\n"
        "  gnome_assembler\n\n"
        "Run cleaning from the command line:\n"
        "  gnome_assembler clean --reads INPUT.fastq --output OUTPUT_PREFIX [options]\n\n"
        "Options:\n"
        "  --min-quality NUMBER       Minimum Phred+33 quality at the 3' end (default: 20)\n"
        "  --min-read-length NUMBER   Discard shorter reads after trimming (default: 50)\n"
        "  --adapter-prefix DNA       Remove this exact adapter from the 5' end\n"
        "  --adapter-suffix DNA       Remove this exact adapter from the 3' end\n";
}

} // namespace gnome

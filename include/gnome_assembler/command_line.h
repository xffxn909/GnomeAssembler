#pragma once

#include "gnome_assembler/cleaning_options.h"

#include <string>

namespace gnome {

struct ParsedCommand {
    CleaningOptions options;
    bool has_cleaning_request = false;
    bool show_help = false;
    std::string error;
};

ParsedCommand parse_command_line(int argc, char* argv[]);
std::string command_help();

} // namespace gnome

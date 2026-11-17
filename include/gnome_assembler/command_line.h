#pragma once

#include "gnome_assembler/cleaning_options.h"
#include "gnome_assembler/graph_pipeline.h"
#include "gnome_assembler/kmer_options.h"

#include <string>

namespace gnome {

struct ParsedCommand {
    CleaningOptions options;
    GraphPipelineOptions graph_options;
    bool has_cleaning_request = false;
    bool has_graph_request = false;
    bool has_assembly_request = false;
    bool show_help = false;
    std::string error;
};

ParsedCommand parse_command_line(int argc, char* argv[]);
std::string command_help();

} // namespace gnome

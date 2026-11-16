#pragma once

#include "gnome_assembler/cleaning_options.h"
#include "gnome_assembler/cleaner.h"

#include <string>

namespace gnome {

struct PipelineResult {
    std::string summary;
    std::string error;
};

PipelineResult run_cleaning_to_files(const CleaningOptions& options);
std::string format_cleaning_summary(const CleaningStatistics& statistics,
                                    const CleaningOptions& options);

} // namespace gnome

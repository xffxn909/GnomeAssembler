#pragma once

#include "gnome_assembler/cleaning_options.h"
#include "gnome_assembler/fastq_record.h"

#include <cstddef>

namespace gnome {

struct CleaningStatistics {
    std::size_t input_reads = 0;
    std::size_t retained_reads = 0;
    std::size_t trimmed_reads = 0;
    std::size_t discarded_reads = 0;
    std::size_t input_bases = 0;
    std::size_t retained_bases = 0;
};

struct CleaningResult {
    bool retained = false;
    bool was_trimmed = false;
    std::size_t input_bases = 0;
    FastqRecord record;
};

int phred_score(char quality_character);
CleaningResult clean_record(FastqRecord record, const CleaningOptions& options);
void update_statistics(CleaningStatistics& statistics, const CleaningResult& result);

} // namespace gnome

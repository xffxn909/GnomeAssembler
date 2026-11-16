#pragma once

#include "gnome_assembler/fastq_record.h"

#include <istream>
#include <string>

namespace gnome {

class FastqReader {
public:
    explicit FastqReader(std::istream& input);
    bool next(FastqRecord& record);

private:
    std::istream& input_;
    std::size_t record_number_ = 0;
};

bool is_valid_dna_sequence(const std::string& sequence);

} // namespace gnome

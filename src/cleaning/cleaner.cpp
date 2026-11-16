#include "gnome_assembler/cleaner.h"

#include <algorithm>
#include <cctype>
#include <utility>

namespace gnome {

namespace {

bool starts_with(const std::string& value, const std::string& prefix) {
    return value.size() >= prefix.size() && value.compare(0, prefix.size(), prefix) == 0;
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string uppercase(std::string value) {
    for (char& character : value) {
        character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
    }
    return value;
}

void erase_prefix(FastqRecord& record, std::size_t length) {
    record.sequence.erase(0, length);
    record.quality.erase(0, length);
}

void erase_suffix(FastqRecord& record, std::size_t length) {
    record.sequence.erase(record.sequence.size() - length);
    record.quality.erase(record.quality.size() - length);
}

} // namespace

int phred_score(char quality_character) {
    return static_cast<unsigned char>(quality_character) - 33;
}

CleaningResult clean_record(FastqRecord record, const CleaningOptions& options) {
    CleaningResult result;
    result.input_bases = record.sequence.size();
    result.record = std::move(record);
    result.record.sequence = uppercase(std::move(result.record.sequence));

    const std::string adapter_prefix = uppercase(options.adapter_prefix);
    const std::string adapter_suffix = uppercase(options.adapter_suffix);
    if (!adapter_prefix.empty() && starts_with(result.record.sequence, adapter_prefix)) {
        erase_prefix(result.record, adapter_prefix.size());
        result.was_trimmed = true;
    }
    if (!adapter_suffix.empty() && ends_with(result.record.sequence, adapter_suffix)) {
        erase_suffix(result.record, adapter_suffix.size());
        result.was_trimmed = true;
    }

    std::size_t retained_length = result.record.quality.size();
    while (retained_length > 0 &&
           phred_score(result.record.quality[retained_length - 1]) < options.minimum_quality) {
        --retained_length;
    }
    if (retained_length != result.record.quality.size()) {
        result.record.sequence.resize(retained_length);
        result.record.quality.resize(retained_length);
        result.was_trimmed = true;
    }

    result.retained = result.record.sequence.size() >= options.minimum_read_length;
    return result;
}

void update_statistics(CleaningStatistics& statistics, const CleaningResult& result) {
    ++statistics.input_reads;
    statistics.input_bases += result.input_bases;
    if (result.was_trimmed) {
        ++statistics.trimmed_reads;
    }
    if (result.retained) {
        ++statistics.retained_reads;
        statistics.retained_bases += result.record.sequence.size();
    } else {
        ++statistics.discarded_reads;
    }
}

} // namespace gnome

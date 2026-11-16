#include "gnome_assembler/fastq_reader.h"
#include "gnome_assembler/cleaner.h"

#include <algorithm>
#include <stdexcept>

namespace gnome {

namespace {

std::string read_required_line(std::istream& input, std::size_t record_number,
                               const char* description) {
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("FASTQ record " + std::to_string(record_number) +
                                 " is missing its " + description + " line.");
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    return line;
}

} // namespace

FastqReader::FastqReader(std::istream& input) : input_(input) {}

bool FastqReader::next(FastqRecord& record) {
    std::string header;
    if (!std::getline(input_, header)) {
        return false;
    }
    if (!header.empty() && header.back() == '\r') {
        header.pop_back();
    }

    ++record_number_;
    if (header.empty() || header.front() != '@') {
        throw std::runtime_error("FASTQ record " + std::to_string(record_number_) +
                                 " must begin with '@'.");
    }

    record.name = header.substr(1);
    record.sequence = read_required_line(input_, record_number_, "sequence");
    const std::string separator = read_required_line(input_, record_number_, "separator");
    record.quality = read_required_line(input_, record_number_, "quality");
    record.source_record = record_number_;

    if (separator.empty() || separator.front() != '+') {
        throw std::runtime_error("FASTQ record " + std::to_string(record_number_) +
                                 " must have a '+' separator line.");
    }
    if (record.sequence.size() != record.quality.size()) {
        throw std::runtime_error("FASTQ record " + std::to_string(record_number_) +
                                 " has different sequence and quality lengths.");
    }
    if (!is_valid_dna_sequence(record.sequence)) {
        throw std::runtime_error("FASTQ record " + std::to_string(record_number_) +
                                 " contains an unsupported nucleotide.");
    }
    for (const char character : record.quality) {
        if (phred_score(character) < 0 || phred_score(character) > 93) {
            throw std::runtime_error("FASTQ record " + std::to_string(record_number_) +
                                     " contains a non-Phred+33 quality character.");
        }
    }
    return true;
}

bool is_valid_dna_sequence(const std::string& sequence) {
    return std::all_of(sequence.begin(), sequence.end(), [](unsigned char base) {
        return std::string("ACGTNacgtn").find(static_cast<char>(base)) != std::string::npos;
    });
}

} // namespace gnome

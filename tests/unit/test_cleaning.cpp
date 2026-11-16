#include "gnome_assembler/cleaner.h"
#include "gnome_assembler/fastq_reader.h"
#include "gnome_assembler/output_writer.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void check_throws(const std::function<void()>& action, const std::string& expected_text) {
    try {
        action();
        check(false, "Expected an exception containing: " + expected_text);
    } catch (const std::runtime_error& error) {
        check(std::string(error.what()).find(expected_text) != std::string::npos,
              "Exception should contain: " + expected_text);
    }
}

gnome::FastqRecord record(std::string sequence, std::string quality) {
    return {"read", std::move(sequence), std::move(quality), 1};
}

void test_fastq_reader() {
    std::istringstream input("@read1\nacgtN\n+\nI!5I!\n");
    gnome::FastqReader reader(input);
    gnome::FastqRecord parsed;
    check(reader.next(parsed), "A valid FASTQ record should be read.");
    check(parsed.name == "read1", "Header name should exclude @.");
    check(parsed.sequence == "acgtN", "Sequence should be preserved by reader.");
    check(!reader.next(parsed), "Reader should end after the final record.");

    std::istringstream no_header("read1\nACGT\n+\nIIII\n");
    gnome::FastqReader invalid_header(no_header);
    check_throws([&] { gnome::FastqRecord value; invalid_header.next(value); }, "begin with '@'");

    std::istringstream wrong_lengths("@read1\nACGT\n+\nIII\n");
    gnome::FastqReader invalid_lengths(wrong_lengths);
    check_throws([&] { gnome::FastqRecord value; invalid_lengths.next(value); }, "different sequence");

    std::istringstream invalid_base("@read1\nACGX\n+\nIIII\n");
    gnome::FastqReader invalid_nucleotide(invalid_base);
    check_throws([&] { gnome::FastqRecord value; invalid_nucleotide.next(value); }, "unsupported nucleotide");
}

void test_quality_and_adapters() {
    gnome::CleaningOptions options;
    options.minimum_quality = 20;
    options.minimum_read_length = 3;
    options.adapter_prefix = "aa";
    options.adapter_suffix = "tt";
    const auto cleaned = gnome::clean_record(record("aaCGTtt", "IIIII!!"), options);
    check(cleaned.retained, "Read remaining at the minimum length should be retained.");
    check(cleaned.was_trimmed, "Adapter removal should mark the read trimmed.");
    check(cleaned.record.sequence == "CGT", "Both adapters should be removed case-insensitively.");
    check(cleaned.record.quality == "III", "Quality must stay aligned with bases.");
    check(gnome::phred_score('5') == 20, "Phred+33 quality boundary should be correct.");
    check(gnome::phred_score('4') == 19, "One quality point below the boundary should be correct.");
}

void test_quality_trimming_and_filtering() {
    gnome::CleaningOptions options;
    options.minimum_quality = 20;
    options.minimum_read_length = 3;
    const auto retained = gnome::clean_record(record("ACGT", "II54"), options);
    check(retained.retained, "A three-base read should remain after low-quality end trimming.");
    check(retained.record.sequence == "ACG", "Only low-quality 3' bases should be trimmed.");

    const auto discarded = gnome::clean_record(record("ACGT", "II4!"), options);
    check(!discarded.retained, "A read below the minimum length should be discarded.");
    check(discarded.record.sequence == "AC", "Trimming should remove every trailing low-quality base.");

    gnome::CleaningStatistics statistics;
    gnome::update_statistics(statistics, retained);
    gnome::update_statistics(statistics, discarded);
    check(statistics.input_reads == 2 && statistics.retained_reads == 1 &&
              statistics.discarded_reads == 1 && statistics.trimmed_reads == 2,
          "Cleaning summary counts should be reproducible.");
    check(statistics.input_bases == 8 && statistics.retained_bases == 3,
          "Base totals should use input and retained lengths.");
}

void test_output_folder_is_accepted() {
    const std::filesystem::path output_folder =
        std::filesystem::temp_directory_path() / "gnome_assembler_milestone1_test";
    std::filesystem::create_directories(output_folder);
    const std::filesystem::path input_path = output_folder / "input.fastq";
    {
        std::ofstream input(input_path);
        input << "@read\nACGT\n+\nIIII\n";
    }

    gnome::CleaningOptions options;
    options.reads_path = input_path.string();
    options.output_prefix = output_folder.string();
    options.minimum_read_length = 4;
    const auto result = gnome::run_cleaning_to_files(options);
    check(result.error.empty(), "An existing output folder should be accepted.");
    check(std::filesystem::exists(output_folder / "cleaned_reads.cleaned.fastq"),
          "The cleaned FASTQ should be written inside the selected folder.");

    options.output_prefix = (output_folder / "new_output_folder" / "sample").string();
    const auto new_folder_result = gnome::run_cleaning_to_files(options);
    check(new_folder_result.error.empty(), "A new output folder should be created when needed.");
    check(std::filesystem::exists(output_folder / "new_output_folder" / "sample.cleaned.fastq"),
          "The cleaned FASTQ should be written to a newly created folder.");

}

} // namespace

int main() {
    test_fastq_reader();
    test_quality_and_adapters();
    test_quality_trimming_and_filtering();
    test_output_folder_is_accepted();
    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed.\n";
        return 1;
    }
    std::cout << "All cleaning tests passed.\n";
    return 0;
}

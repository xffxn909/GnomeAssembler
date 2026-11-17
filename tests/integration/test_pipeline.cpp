#include "gnome_assembler/assembly_pipeline.h"

#include <filesystem>
#include <iostream>

#ifndef GNOME_SOURCE_DIR
#error "GNOME_SOURCE_DIR must be configured for the integration test"
#endif

namespace {

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

} // namespace

int main() {
    namespace fs = std::filesystem;
    const fs::path source(GNOME_SOURCE_DIR);
    const fs::path input = source / "tests" / "fixtures" / "simulated_circular_virus.fastq";
    const fs::path output = fs::temp_directory_path() / "gnome_assembler_integration" / "simulated_virus";
    std::error_code ec;
    fs::remove_all(output.parent_path(), ec);
    fs::create_directories(output.parent_path(), ec);
    check(!ec, "Temporary integration output folder should be creatable.");

    gnome::GraphPipelineOptions options;
    options.cleaning.reads_path = input.string();
    options.cleaning.minimum_quality = 0;
    options.cleaning.minimum_read_length = 0;
    options.kmer.kmer_length = 7;
    options.kmer.min_kmer_abundance = 1;
    options.output_prefix = output.string();

    const auto result = gnome::run_assembly_pipeline(options);
    check(result.error.empty(), "Simulated circular-virus pipeline should finish.");
    check(result.verification.assembly_complete, "Verification should establish complete circular assembly.");
    check(result.verification.covered_reads == 1, "The simulated read should be represented.");
    check(result.assembly.contigs.size() == 1 && result.assembly.contigs.front().circular,
          "The simulated virus should assemble to one circular contig.");
    check(fs::exists(output.string() + ".contigs.fasta"), "Contig FASTA should be written.");
    check(fs::exists(output.string() + ".assembly-report.txt"), "Text report should be written.");
    check(fs::exists(output.string() + ".assembly-report.json"), "JSON report should be written.");

    fs::remove_all(output.parent_path(), ec);
    if (failures != 0) return 1;
    std::cout << "Simulated circular-virus integration test passed.\n";
    return 0;
}

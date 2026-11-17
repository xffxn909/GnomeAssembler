# GnomeAssembler

A Windows-focused C++ de novo virus genome assembler based on a de Bruijn graph. It provides a simple Qt 6 Widgets desktop interface and a command-line mode. It will read quality-scored FASTQ reads, clean them, construct and traverse a graph, verify the resulting assembly, and identify the virus with a locally installed BLAST+ database.

The project name and GitHub repository retain the requested `GnomeAssembler` spelling; the program's biological scope is genome assembly.

## Project status

- **Milestone 1**: FASTQ quality trimming, adapter removal, filtering, and summary reporting via GUI and CLI.
- **Milestone 2**: k-mer counting, abundance filtering, de Bruijn graph construction, deterministic graph dumps, and asynchronous background execution in the GUI.
- **Milestone 3**: deterministic iterative Hierholzer traversal, linear and circular contig spelling, explicit invalid-graph diagnostics, contig FASTA output, and GUI assembly results.
- **Milestone 4**: edge-use and read-to-contig verification, completeness status, N50 metrics, text/JSON assembly reports, and an output-folder action in the GUI.

Later assembly work is described in [the implementation plan](docs/IMPLEMENTATION_PLAN.md).

## Technology choices

- C++17, CMake 3.20 or newer, and Qt 6 Widgets
- Windows development with Visual Studio Build Tools or Visual Studio
- FASTQ as the initial input format, preserving per-base quality scores for trimming
- BLAST+ `blastn` launched by the program to identify the finished assembly against a configured local nucleotide database

## Build and open the GUI

```powershell
cmake -S . -B build -DCMAKE_PREFIX_PATH=C:\Qt\6.12.0\mingw_64
cmake --build build
.\build\Debug\gnome_assembler.exe
```

For a MinGW/Ninja single-configuration build, the executable is normally `build\gnome_assembler.exe`. Opening it without arguments shows the FASTQ-cleaning window.

## Clean FASTQ reads from the command line

```powershell
gnome_assembler clean `
  --reads data\sample.fastq `
  --min-quality 20 `
  --min-read-length 50 `
  --output results\sample
```

This writes `sample.cleaned.fastq`, `sample.cleaned.fasta`, and `sample.cleaning-summary.txt`. Input reads, reference data, BLAST databases, and generated results are intentionally not tracked by Git.

In the GUI, **Folder...** accepts an existing writable folder and writes files named `cleaned_reads.cleaned.*` inside it. Selecting an input FASTQ also suggests a `GnomeAssemblerResults` folder beside that input. You may instead type a full output prefix when you want a custom filename.

## Build de Bruijn graph from the command line

```powershell
gnome_assembler build-graph `
  --reads data\sample.fastq `
  --kmer-length 21 `
  --min-kmer-abundance 2 `
  --output results\sample
```

This cleans the input reads, extracts k-mers, filters low-abundance k-mers, constructs the de Bruijn graph, and writes `sample.graph-dump.txt`.

## Assemble contigs from the command line

```powershell
gnome_assembler assemble `
  --reads data\sample.fastq `
  --kmer-length 21 `
  --min-kmer-abundance 2 `
  --output results\sample
```

This runs the same cleaning and graph stages, traverses a single valid Eulerian path or cycle, and writes `sample.contigs.fasta`, `sample.assembly-report.txt`, and `sample.assembly-report.json`. The report checks that every retained graph edge was traversed, counts cleaned reads represented by a contig or its reverse complement, and reports contig count, total length, longest contig, N50, coverage, and ambiguity warnings. Circular contigs are labeled as circular and do not imply a biological start base. If the graph is disconnected or has invalid degree balance, the report explains why completeness was not established.

## Documentation

- [Complete implementation plan](docs/IMPLEMENTATION_PLAN.md)

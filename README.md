# GnomeAssembler

A Windows-focused C++ de novo virus genome assembler based on a de Bruijn graph. It will read quality-scored FASTQ reads, clean them, construct and traverse a graph, verify the resulting assembly, and identify the virus with a locally installed BLAST+ database.

The project name and GitHub repository retain the requested `GnomeAssembler` spelling; the program's biological scope is genome assembly.

## Project status

The repository foundation and implementation plan are in place. The executable currently builds as a placeholder while the pipeline is implemented in the milestones described in [the implementation plan](docs/IMPLEMENTATION_PLAN.md).

## Technology choices

- C++17 and CMake 3.20 or newer
- Windows development with Visual Studio Build Tools or Visual Studio
- FASTQ as the initial input format, preserving per-base quality scores for trimming
- BLAST+ `blastn` launched by the program to identify the finished assembly against a configured local nucleotide database

## Build the current foundation

```powershell
cmake -S . -B build
cmake --build build --config Debug
.\build\Debug\gnome_assembler.exe
```

For single-configuration generators, the executable may instead be at `build\gnome_assembler.exe`.

## Intended command shape

```powershell
gnome_assembler assemble `
  --reads data\sample.fastq `
  --kmer-size 31 `
  --min-quality 20 `
  --min-read-length 50 `
  --blast-db data\blast_db\viruses `
  --output results\sample
```

The exact command-line interface will be delivered and tested in Milestone 1. Input reads, reference data, BLAST databases, and generated results are intentionally not tracked by Git.

## Documentation

- [Complete implementation plan](docs/IMPLEMENTATION_PLAN.md)

# GnomeAssembler user guide

## Scope and limitations

GnomeAssembler is a teaching-oriented, de novo assembler for small, single-virus FASTQ datasets. It uses a de Bruijn graph and does not use a reference sequence to construct contigs. It supports Phred+33 FASTQ input, one k-mer length per run, and a single complete Eulerian component. It is not appropriate for metagenomes, diploid/polyploid genomes, large eukaryotic genomes, or assemblies with unresolved repeats. A successful report means all retained graph edges and cleaned reads were represented; it does not prove biological correctness or uniqueness when ambiguity warnings are present.

## Build on Windows

Install CMake 3.20+, a C++17 compiler, and Qt 6 Widgets. For a Qt MinGW installation:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=C:\Qt\6.12.0\mingw_64
cmake --build build
ctest --test-dir build --output-on-failure
```

For Visual Studio, omit `-G Ninja` and build the desired configuration:

```powershell
cmake -S . -B build-vs -DCMAKE_PREFIX_PATH=C:\Qt\6.12.0\msvc2022_64
cmake --build build-vs --config Debug
ctest --test-dir build-vs -C Debug --output-on-failure
```

Run `gnome_assembler` with no arguments to open the GUI. Use `gnome_assembler --help` for the complete CLI option list.

## Reproducible CLI workflow

Record the input FASTQ checksum, program commit, command line, k-mer settings, and output reports. An example assembly is:

```powershell
.\gnome_assembler.exe assemble `
  --reads data\sample.fastq `
  --output results\sample `
  --min-quality 20 `
  --min-read-length 50 `
  --kmer-length 21 `
  --min-kmer-abundance 2
```

This writes a deterministic graph dump, contig FASTA, and text/JSON assembly reports for the same input and options. The JSON report is intended for workflow automation. Circular contig start positions are deterministic for reporting but are not biological assertions.

## Optional local BLAST+ identification

Install NCBI BLAST+ separately and build or obtain a local, preformatted nucleotide virus database. GnomeAssembler never downloads executables or databases. Supply the exact path to `blastn.exe` and the database prefix:

```powershell
.\gnome_assembler.exe assemble `
  --reads data\sample.fastq `
  --output results\sample `
  --kmer-length 21 `
  --min-kmer-abundance 2 `
  --blastn C:\blast\bin\blastn.exe `
  --blast-db D:\virus_db\viruses `
  --min-blast-identity 90 `
  --min-blast-coverage 80
```

The database prefix must point to a formatted database (for example, files with `.ndb` or legacy `.nin/.nhr/.nsq` companions). The program saves raw tabular BLAST results plus text/JSON identification reports. Qualified hits are ranked by bit score, then e-value, identity, and query coverage. Failed or unavailable identification is reported as `not identified` and does not discard assembly results.

## GUI workflow

Choose FASTQ input and output prefix, select cleaning and k-mer settings, and optionally fill in the local BLAST fields. **Assemble Contigs** runs the shared pipeline in a worker thread. The run log contains assembly verification and BLAST status. **Open Output Folder** opens the report location after a run.

## Test data provenance

Tracked FASTQ fixtures are synthetic and exist only to exercise parser, graph, circular traversal, verification, and integration behavior. `simulated_circular_virus.fastq` is a small artificial circular-read case; it is not phiX174 or a biological benchmark. Users must provide and retain provenance for their own sequencing data and local BLAST database.

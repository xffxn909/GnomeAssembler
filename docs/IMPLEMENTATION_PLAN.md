# Implementation Plan

## 1. Scope and decisions

GnomeAssembler is a Windows-first C++17 application for **de novo assembly of small virus genomes**. It has a simple Qt 6 Widgets desktop GUI for normal use and a CLI for reproducible scripted runs. Both interfaces call the same Qt-independent C++ pipeline.

The assembler will use a de Bruijn graph and will not use a reference to construct an assembly. It initially accepts standard FASTQ with Phred+33 qualities. The target is small, single-genome viral data, including circular genomes such as phiX174; it is not a production metagenomic, diploid, or large-eukaryotic assembler.

After assembly verification, the program will optionally run local BLAST+ `blastn` against a user-selected preformatted virus database. It never silently downloads a database.

## 2. Technology decisions

- C++17, CMake 3.20+, and Qt 6 Widgets.
- Qt Widgets is selected because it is established, simple to learn, CMake-friendly, and suitable for a clear Windows desktop interface.
- Core biological and file-processing code must not depend on Qt. GUI and CLI only convert user input to typed options and display the core result.
- The GUI launches when the executable has no pipeline command. CLI commands remain available for automation.

## 3. Deliverables

1. `gnome_assembler.exe`, built with CMake on Windows and Qt 6.
2. A shared FASTQ cleaning, k-mer, graph, assembly, verification, reporting, and optional BLAST pipeline.
3. A GUI to select files, adjust supported settings, run the pipeline, and view status, warnings, and output locations.
4. A CLI that exposes the same pipeline settings and useful exit codes.
5. Automated unit and integration tests, including a simulated phiX174-style circular dataset.
6. Documentation for build, GUI, CLI, BLAST setup, limitations, and reproducibility.

## 4. Repository layout

```text
GnomeAssembler/
├── CMakeLists.txt
├── README.md
├── docs/IMPLEMENTATION_PLAN.md
├── include/gnome_assembler/       # shared public core interfaces
├── src/
│   ├── main.cpp
│   ├── cli/                       # command parsing only
│   ├── gui/                       # Qt Widgets presentation only
│   ├── io/
│   ├── cleaning/
│   ├── kmer/
│   ├── graph/
│   ├── assembly/
│   ├── verification/
│   ├── blast/
│   └── reporting/
├── tests/unit/
├── tests/integration/
├── tests/fixtures/
└── data/                          # ignored local reads, databases, results
```

## 5. Architecture

### 5.1 GUI and CLI

- The initial window is intentionally small: FASTQ input picker, output-prefix/folder picker, minimum-quality and minimum-read-length controls, optional exact 5'/3' adapter fields, Run button, and read-only log.
- The GUI shows plain-language errors and all generated output paths. It never reimplements cleaning or assembly logic.
- CLI invocation is `gnome_assembler clean --reads INPUT.fastq --output PREFIX [options]`.
- Milestone 2 moves long-running calls off the GUI thread, disables Run while active, and provides progress/status. This keeps the Milestone 1 design easy to understand while preventing later UI freezes.

### 5.2 Input and cleaning

- Stream FASTQ records rather than loading a full file.
- Validate four-line records, sequence/quality length equality, DNA alphabet, and Phred+33 quality range; include the record number in errors.
- Normalize retained sequence to uppercase DNA bases.
- Remove configured exact adapter prefixes/suffixes, trim low-quality 3' ends at a configurable threshold, and discard reads below the configured minimum surviving length.
- Write cleaned FASTQ, cleaned FASTA, and a deterministic text summary of input, retained, trimmed, discarded, and base counts.

### 5.3 K-mers and graph

- Extract fixed-length k-mers from retained reads, skipping ambiguous-base k-mers initially.
- Count with `std::unordered_map<std::string, uint32_t>` and discard counts below a configurable abundance threshold.
- Represent each retained k-mer as a directed edge from `(k-1)` prefix to `(k-1)` suffix. Maintain adjacency, degrees, multiplicity, and edge usage.
- Identify connected components, tips, branches, degree imbalance, and disconnected graphs.

### 5.4 Assembly

- Check Eulerian path/cycle conditions and traverse valid components with iterative Hierholzer traversal.
- Spell contigs from traversal edges and prove each retained edge was consumed once.
- Treat circular contigs explicitly and canonicalize only for comparison/reporting; do not assert a biological start point.
- Report repeat/branch ambiguity and emit maximal valid contigs rather than falsely claiming uniqueness.

### 5.5 Verification and reporting

- Count cleaned reads that occur exactly in each contig or its reverse complement.
- Calculate contig count, total length, longest contig, N50, covered-read percentage, edge-use status, and ambiguity warnings.
- For known fixtures, compare circular sequences modulo rotation and reverse complement.
- Write human-readable text and machine-readable JSON reports. The GUI displays the concise report and offers an output-folder action.

### 5.6 BLAST+ identification

- Require a user-configured local `blastn` executable and local preformatted nucleotide virus database.
- Write contigs to FASTA, invoke `blastn` with a safe argument list, and parse deterministic tabular output.
- Select the strongest qualified hit by bit score, then e-value, identity, and coverage.
- Keep assembly outputs if identification cannot run, and show an explicit `not identified` status.

## 6. Milestones and acceptance criteria

### Milestone 1 — Foundation, FASTQ cleaning, and simple GUI

- Implement typed configuration, CLI parsing, FASTQ reader, Phred conversion, adapter removal, 3' quality trimming, filtering, output writers, and summary logging.
- Implement the Qt 6 Widgets cleaning window with the same supported options and result log.
- Add tests for valid and malformed FASTQ, quality boundaries, adapters, minimum length, and summary statistics.
- **Acceptance:** a fixture produces reproducible cleaned FASTQ/FASTA and summary through CLI; the GUI selects the same input, invokes the shared pipeline, and visibly reports the same outcome and output paths.

### Milestone 2 — K-mers and graph construction

- Implement k-mer counting/filtering and graph nodes, edges, multiplicities, and degrees.
- Add deterministic graph dumps and tests for hand-worked examples.
- Add k-mer settings, background execution, and progress to the GUI.
- **Acceptance:** hand-worked reads produce the exact expected k-mer table and graph.

### Milestone 3 — Traversal and contigs — complete

- Implement iterative Hierholzer traversal, contig spelling, Euler path/cycle handling, and clear invalid-graph diagnostics.
- Show contigs and ambiguity warnings in the GUI using the shared report data.
- **Acceptance:** synthetic linear and circular datasets reconstruct expected sequences, allowing equivalent circular rotations.

Implementation notes: traversal uses an iterative, deterministic Hierholzer algorithm; k-mer multiplicity is traversed as parallel edges, so all counted edge occurrences are consumed. A valid cycle is emitted without its duplicated closing `(k-1)` suffix and is explicitly marked circular. Disconnected or non-Eulerian graphs stop with diagnostics rather than an incomplete contig. The CLI command is `assemble`; it writes a contig FASTA alongside the graph dump. The GUI's **Assemble Contigs** button uses the same pipeline and shows contig lengths, topology, and ambiguity warnings in its log.

### Milestone 4 — Verification and reporting — complete

- Implement edge-use checks, read/contig coverage, ambiguity warnings, N50, JSON, and text summaries.
- Add simulated viral integration tests starting with phiX174.
- Add a concise GUI report panel and output-folder action.
- **Acceptance:** the report establishes completeness or explains a non-unique result.

Implementation notes: verification checks the number of retained graph edges (including multiplicity) against the traversal count, and searches every cleaned read and its reverse complement in the contigs. Circular contigs are searched across their reporting boundary. The pipeline writes deterministic text and JSON reports with completion status, edge-use status, read coverage, contig metrics, N50, and warnings. The GUI log displays the same concise text report and provides an output-folder action.

### Milestone 5 — BLAST+ identification — complete

- Implement executable/database discovery, FASTA output, process invocation, result parsing, thresholds, and failure messages.
- Add GUI selectors for BLAST/database settings and visible identification status.
- **Acceptance:** a known-virus contig produces an auditable identification report with selected-hit metrics.

Implementation notes: `assemble` accepts optional local `--blastn` and `--blast-db` settings plus identity and coverage thresholds. The program validates the executable and formatted database prefix before invoking Windows `blastn` with a fixed argument vector and tabular output format. It retains all assembly outputs on a BLAST failure, reports `not identified`, and writes raw tabular results plus text/JSON identification reports. Saved-tabular parser tests cover thresholding and deterministic best-hit selection.

### Milestone 6 — Hardening and submission — complete

- Document Qt build configuration, GUI/CLI workflow, BLAST setup, data provenance, limitations, and reproducibility.
- Test Debug and Release builds and assess runtime/memory on the target viral dataset.
- **Acceptance:** a new Windows user can build, run the sample workflow in either interface, inspect the report, and reproduce the result.

Implementation notes: the user guide documents Qt/CMake builds for Ninja and Visual Studio, CLI/GUI workflows, local BLAST configuration, data provenance, limitations, and reproducibility records. Automated tests cover all core milestones and an end-to-end synthetic circular-virus run that writes contig and verification artifacts. Debug and Release configurations are built and tested as part of the release verification.

## 7. Test strategy

| Layer | Test | Expected result |
| --- | --- | --- |
| FASTQ | valid, malformed, and low-quality fixtures | valid reads parse; invalid records are diagnosed; trimming is correct |
| Cleaning | quality threshold, adapters, minimum length, statistics | bases/qualities stay aligned and counts are exact |
| GUI | manual Milestone 1 smoke test; GUI-to-core test where practical | chosen values reach the shared core and log is readable |
| K-mers | hand-calculated reads | exact counts and filtering |
| Graph | linear, cyclic, branching examples | expected degrees and warnings |
| Traversal | Euler path/cycle and non-Euler graph | complete traversal or explicit failure |
| Assembly | generated viral reads | reference match allowing circular rotation/reverse complement |
| Verification | omitted edge, duplicate edge, missing reads | each defect is detected |
| BLAST | saved tabular fixtures and optional local run | correct best hit and threshold handling |

## 8. Complexity targets

With `B` retained bases, `K` retained distinct k-mers, and `E` retained graph edges: FASTQ processing and k-mer extraction are expected `O(B)`; graph construction and Hierholzer traversal are `O(E)`; memory is initially `O(K + E)`. GUI work must not alter these core targets.

## 9. Risks and controls

| Risk | Control |
| --- | --- |
| Low-quality reads create false branches | configurable trimming/filtering and discarded-data statistics |
| Repeats prevent a unique genome | warn explicitly and preserve contigs |
| Incorrect k fragments/collapses graph | expose k and document/test trade-offs |
| Missing/outdated BLAST database | require explicit local path and log it |
| Circular start differs from reference | rotation- and reverse-complement-aware comparison |
| External tool fails | retain assembly/report output and state failure clearly |
| GUI and CLI diverge | keep all pipeline logic in Qt-independent shared code |
| GUI freezes on larger inputs | use a worker thread and progress UI from Milestone 2 |
| Qt is unavailable at build time | document the Qt CMake prefix path and surface CMake's clear configuration error |

## 10. Project status

All planned milestones are implemented. Future work should be versioned as new requirements rather than changing the completed milestone scope.

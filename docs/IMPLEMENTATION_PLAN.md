# Implementation Plan

## 1. Scope and decisions

This is a Windows-first C++17 command-line application for **de novo assembly of small virus genomes**. It will use a de Bruijn graph and will not use a reference sequence to construct the assembly. BLAST+ is part of the application workflow: after verification, the program will run `blastn` against a user-configured **local** virus database and include the best hit in its report.

Initial input is standard FASTQ with Phred+33 quality encoding. FASTQ is preferred over FASTA at this stage because it retains the quality scores required for the proposal's cleaning stage. A later FASTA input mode may be added only after the FASTQ pipeline is verified.

The initial target is small, single-genome viral data (including circular genomes such as phiX174). It is not intended to be a production assembler for mixed metagenomic, diploid, or large eukaryotic genomes.

## 2. Deliverables

1. `gnome_assembler.exe`, built with CMake on Windows.
2. A FASTQ cleaning, k-mer counting, graph construction, assembly, verification, and reporting pipeline.
3. A BLAST+ integration that writes the assembled contig to FASTA, invokes `blastn`, parses tabular output, and reports the strongest qualified match.
4. Automated unit and integration tests, including a simulated phiX174-style circular-genome dataset.
5. Example command lines and an explanation of inputs, outputs, limitations, and reproducibility.

## 3. Proposed repository layout

```text
GnomeAssembler/
├── CMakeLists.txt
├── README.md
├── docs/
│   └── IMPLEMENTATION_PLAN.md
├── src/
│   ├── main.cpp
│   ├── cli/
│   ├── io/
│   ├── cleaning/
│   ├── kmer/
│   ├── graph/
│   ├── assembly/
│   ├── verification/
│   ├── blast/
│   └── reporting/
├── include/gnome_assembler/
├── tests/
│   ├── unit/
│   ├── integration/
│   └── fixtures/
├── scripts/
└── data/                         # ignored: local reads, databases, results
```

## 4. Architecture

### 4.1 Input and cleaning

- Stream FASTQ records rather than loading the entire file by default.
- Validate sequence/quality length equality and nucleotide alphabet; report malformed record locations.
- Remove configured adapter suffixes/prefixes.
- Trim low-quality ends using a configurable Phred threshold.
- Drop reads below the configured minimum surviving length.
- Record input, retained, trimmed, discarded, and base-count statistics.

### 4.2 K-mers

- Normalize each retained read to uppercase DNA bases.
- Extract fixed-length k-mers, skipping k-mers with ambiguous bases in the first version.
- Store counts in `std::unordered_map<std::string, uint32_t>`.
- Remove k-mers below a configurable abundance threshold before graph creation, reducing isolated sequencing errors.
- Expose `k` and minimum abundance through the CLI and record them in every report.

### 4.3 De Bruijn graph

- Represent each retained k-mer as a directed edge from its `(k-1)` prefix to its `(k-1)` suffix.
- Maintain adjacency lists, per-node in/out degree, and per-edge multiplicity/usage state.
- Build connected components and flag disconnected graphs, tips, branches, and degree imbalance.
- Preserve enough edge data to prove traversal coverage during verification.

### 4.4 Assembly

- Identify whether each usable component has Eulerian-path or Eulerian-cycle degree conditions.
- Traverse valid components with iterative Hierholzer traversal, consuming each retained edge exactly once.
- Spell contigs from node labels and traversal edges.
- Treat circular assemblies explicitly: canonicalize a circular sequence only for comparison/reporting, without claiming a biological origin position.
- Report branch/repeat ambiguity rather than silently claiming uniqueness. Initial behavior is to emit maximal valid contigs and warnings; sophisticated repeat resolution is outside the first release.

### 4.5 Verification and metrics

- Confirm degree conditions before assembly and that each expected graph edge was consumed exactly once.
- Count original cleaned reads that occur exactly in a contig or its reverse complement; report coverage rather than falsely treating imperfect/error reads as failures.
- Calculate contig count, total assembled length, longest contig, N50, covered-read percentage, and ambiguous/repeat warnings.
- For known test references, compare linear sequences directly and circular sequences modulo rotation and reverse complement.
- Report contamination signals: unassembled high-abundance components and BLAST hits inconsistent with the leading result.

### 4.6 BLAST+ identification

- Require `blastn` plus a local, preformatted nucleotide virus database selected with `--blast-db`; do not download databases silently.
- Write assembled contigs to an output FASTA file.
- Invoke `blastn` with a safe argument list and tabular output fields covering accession, title, percent identity, alignment length, e-value, bit score, query coverage, and taxonomy fields when the database provides them.
- Check process exit status and parse the best hit deterministically: strongest bit score, then e-value, identity, and coverage.
- Produce a clear `not identified` result when BLAST+ is missing, the database is absent, no hit meets configured thresholds, or execution fails. Assembly output remains available in all cases.

## 5. Milestones and acceptance criteria

### Milestone 1 — Foundation and FASTQ cleaning

- Implement the command-line parser, typed configuration, logging, FASTQ reader, quality conversion, trimming, and filtering.
- Add unit tests for malformed FASTQ, quality boundaries, adapters, and minimum length.
- **Acceptance:** a small FASTQ fixture produces a reproducible cleaning summary and cleaned FASTQ/FASTA output.

### Milestone 2 — K-mers and graph construction

- Implement k-mer counting/filtering and graph nodes/edges/degrees.
- Add deterministic graph dumps for debugging and unit tests for counts, multiplicity, and degree calculations.
- **Acceptance:** hand-worked read examples produce the exact expected k-mer table and graph.

### Milestone 3 — Traversal and contigs

- Implement iterative Hierholzer traversal and contig spelling.
- Handle Eulerian paths, cycles, disconnected components, and invalid degree patterns with useful diagnostics.
- **Acceptance:** synthetic linear and circular datasets reconstruct expected sequences, including equivalent circular rotations.

### Milestone 4 — Verification and reporting

- Implement edge-use checks, read/contig coverage, ambiguity warnings, N50, and machine-readable JSON plus human-readable text summaries.
- Add integration tests using simulated short reads from a known small viral reference, beginning with phiX174.
- **Acceptance:** the report establishes whether the fixture assembled completely and explains any non-unique result.

### Milestone 5 — BLAST+ identification

- Implement executable/database discovery, FASTA output, process invocation, result parsing, thresholds, and failure messages.
- Test the parser with saved BLAST tabular fixtures and add an optional end-to-end test when a local database is available.
- **Acceptance:** a known-virus contig produces an auditable identification report containing the selected hit and its alignment metrics.

### Milestone 6 — Hardening and project submission

- Document installation, local BLAST database setup, data provenance, command examples, limitations, and reproducibility steps.
- Run tests under Debug and Release builds; assess memory/runtime on the chosen viral dataset.
- **Acceptance:** a new Windows user can build, run the sample workflow, inspect the verification report, and reproduce the stated result.

## 6. Test strategy

| Layer | Test | Expected result |
| --- | --- | --- |
| FASTQ | valid, malformed, and low-quality fixtures | valid reads parsed; invalid records diagnosed; correct trimming |
| K-mers | hand-calculated reads | exact counts and filtering |
| Graph | linear, cyclic, branching examples | expected degrees and warnings |
| Traversal | Euler path/cycle and non-Euler graph | complete traversal or explicit failure |
| Assembly | generated viral reads | sequence matches reference, allowing circular rotation/reverse complement |
| Verification | omitted edge, duplicate edge, missing reads | check detects each defect |
| BLAST | saved outfmt-6 results and optional local run | correct top hit/threshold handling |

## 7. Complexity targets

Let `B` be retained input bases, `K` the number of retained distinct k-mers, and `E` the retained graph-edge count. FASTQ processing and k-mer extraction should be `O(B)` expected time. Hash-map insertion is expected constant time per k-mer. Graph construction is `O(E)`, and Hierholzer traversal is `O(E)`. Memory is dominated by k-mer/graph storage, approximately `O(K + E)` in the initial transparent representation.

## 8. Risks and controls

| Risk | Control |
| --- | --- |
| Low-quality reads introduce false branches | configurable trimming and k-mer abundance filtering; report discarded data |
| Repeats prevent a unique genome | warn explicitly; preserve contigs and do not overstate correctness |
| Incorrect `k` fragments or collapses the graph | expose `k`, document trade-offs, and compare selected values in tests |
| BLAST database is missing/out-of-date | require explicit local path; log database path and failure status |
| Circular genome start differs from reference | use rotation- and reverse-complement-aware comparison |
| External tool failure | retain assembly/report outputs and return a clear nonzero status only for the identification stage |

## 9. Immediate next task

Implement Milestone 1: define the CLI contract, add the FASTQ record/model types, and write the quality-trimming tests before adding graph code.

# nf-core/nanorepertoire: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 2.0.0 - 2026-07-23

Second major release focused on repertoire analysis robustness, reporting refactor, ARM64 and container compatibility, and sustainability reporting.

### `Added`

- Replaced deprecated `getcdr3` path with `nanocdr-x` integration using `run_nanocdrx.py`, preserving pipeline I/O compatibility.
- Added CO2 footprint reporting with `nf-co2footprint` plugin and documentation for generated trace artifacts.
- Added support for remote CSV test input from bucket-based test data.
- Added `<sample>_cdr3_summary.tsv` to `NANOCDRX`, counting the cluster representatives excluded by length, excluded for non-standard residues, sent to nanoCDR-X, and with or without a CDR3.

### `Changed`

- Refactored reporting by introducing `repertoire_report` subworkflow and hybrid reporting integration.
- Parameterized CD-HIT and renamed report stage to `AGGREGATE_STATS` for clearer module semantics.
- Updated pipeline ownership references from `nf-core` to `lescailab` metadata where applicable.
- Refined output routing for pipeline info and CO2-related artifacts to align with official plugin guidance.
- `run_nanocdrx.py` writes every cluster representative to `<sample>_cdr3.tsv` exactly once, with the new statuses `cdr3-too-long`, `excluded-length` and `excluded-nonstandard`. Representatives containing non-standard residues are now excluded instead of having those residues deleted, which shortened their CDR3.
- `run_nanocdrx.py` calls `predict_cdrs` once per sample instead of once per 500 sequences.
- `NANOTRANSLATE` is faster and uses less memory: the start and end motifs are built from their amino acids instead of being listed one by one, each motif is searched once, reads are processed in groups instead of being loaded all at once, and the script uses the CPUs assigned by Nextflow instead of all those of the machine. The output is unchanged (verified byte by byte on two llama samples and an alpaca subset).

### `Fixed`

- Fixed `NANOCDRX` silently discarding whole 500-sequence batches: a single sequence longer than 150 residues, the input length of the nanoCDR-X model, made `predict_cdrs` fail on its batch, and the failure was ignored. The length filter is now 70-150 residues, and a failing or incomplete `predict_cdrs` run stops the step with an error.
- Fixed CDR3s longer than 50 residues being labelled `non-unique`.
- Fixed the nanoCDR-X version recorded in `versions.yml`, which was read from a `--version` option that `predict_cdrs` does not have.
- Fixed `NANOTRANSLATE` discarding merged reads that hold the VHH on the reverse strand, as half of the reads of a non-directional library do: when the start or end motif is missing on the forward strand, the reverse complement is searched too. The translation log reports how many reads were reverse-complemented.
- Fixed translations containing undetermined residues (`X`, from codons with an `N`) being kept: they are now discarded and counted in the translation log, like those with a stop codon. Reads with many `X` formed spurious singleton clusters.
- Fixed `FLASH` writing the merged reads in a different order on every run when it used more than one thread, which changed the read numbering of `RENAME` and the CD-HIT clusters between runs of the same data. FLASH now runs with `-t 1`, the setting its manual gives for deterministic output.
- Fixed the CD-HIT clusters depending on the number of threads: with more than one thread CD-HIT assigns some sequences to different clusters depending on the thread count, and slightly differently on every run. `CDHIT_CDHIT` now runs with one thread, the CD-HIT default. The `-T 0` and `-M 16000` options set in `conf/modules.config` had no effect, because the module appends its own `-T` and `-M`, and were removed.
- Fixed `NANOTRANSLATE` discarding VHH whose first three residues or final `TVSS` differ by one amino acid from the start and end motifs (for example `QVR`, `QGQ`, `TVSA`, `IVSS`): reads with no exact motif on either strand are searched again in the six reading frames with anchors that tolerate one substitution, confirmed by residues conserved in the IMGT alpaca germline (the FR1 cysteine, the second cysteine and the `WG.G` motif of framework 4). The rescued reads are counted in the translation log and listed in a `.rescued` file next to `.notfound` and `.outframe`.
- Fixed an 11-nt end motif in `NANOTRANSLATE` (`CTGTTAGCTCA`, a truncated encoding of `TVSS`): it never matched its intended `ACTGTTAGCTCA` and could match elsewhere in a read, leaving the VHH region out of frame. All end motifs are now the 12-nt encodings of `TVSS`.
- The aggregated tables and the Quarto report list the samples sorted by name, instead of in the order in which the samples finished clustering.
- Fixed the software versions report, which listed only FastQC and the reports: the versions of the `FASTQ_TO_FASTA` and `FASTA_CLUSTERING` steps were never collected. Also fixed the seqtk version recorded by `RENAME` and removed the pandas and plotly versions from `NANOREPERTOIRE_REPORT`, which does not use them.
- Improved Docker compatibility on Apple Silicon and fixed ARM64 Wave container handling for report rendering.
- Fixed resource-limit guard logic (`check_max`) and conda channel checks for Nextflow `25.10.4` + micromamba.
- Addressed nf-core compliance and linting issues across templates/modules/subworkflows.

### `Dependencies`

- Updated nf-core template baseline to `nf-core/tools 3.4.1`.
- Updated module set and subworkflows to newer nf-core revisions.
- Added `nf-co2footprint@1.0.0-beta` plugin.
- Updated Plotly-related container setup for reporting components.
- Aligned the container images with the conda environments: cutadapt 4.9, seqtk 1.4 and biopython 1.78 for every container engine, so that all profiles run the same tool versions.

### `Deprecated`

- Deprecated legacy `getcdr3`-driven CDR3 extraction path in favor of `nanocdr-x`.
- Removed the unused `GETCDR3` module and its test.
- Removed `bin/getcdr3.py`, used only by the removed `GETCDR3` module.

## 1.0.0 - 2025-10-17 Alpaca Lypse

Initial release of nf-core/nanorepertoire, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

No fixes in this initial release

### `Dependencies`

- Uses `nf-schema@2.1.x` plugin for parameter schema validation
- Standardizes container registries on `quay.io`
- Requires Nextflow ≥ 24.10.5

### `Deprecated`

No deprecations in this release

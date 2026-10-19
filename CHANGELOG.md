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
- Added a sequence accounting table to the interactive report, and `report/sequence_accounting.tsv`: for each sample, the number of read pairs, merged reads, translated sequences, clustered sequences, clusters, annotated representatives and distinct CDR3s, with the sequences lost between two steps and the reason. The counts are read from the Cutadapt reports, the translation logs, the cluster sizes and the nanoCDR-X summaries, and the report flags neighbouring counts that do not agree.
- Added `CLUSTER_MOSTFREQUENT` and `bin/cluster_most_frequent.py`, which choose the most frequent sequence of each CD-HIT cluster (ties go to the member CD-HIT lists first) and write it, under the identifier of the CD-HIT representative, to `cluster_most_frequent/<sample>_most_frequent.fasta`, the input of `NANOCDRX`, with a per-cluster table `<sample>_most_frequent.tsv`.

### `Changed`

- Refactored reporting by introducing `repertoire_report` subworkflow and hybrid reporting integration.
- Parameterized CD-HIT and renamed report stage to `AGGREGATE_STATS` for clearer module semantics.
- `CDHIT_CDHIT`, the slowest step of the pipeline, uses up to 16 CPUs, capped by `check_max` and the resource limits of the profile: with `-g 1 -det 1` the clusters no longer depend on the number of threads.
- `CDHIT_CDHIT` runs with `-g 1`, which puts each sequence in the cluster of its most similar representative instead of the first one found above the identity threshold, as in Deschaght et al. 2017. The clusters and their representatives are unchanged; the cluster sizes change (on the llama `before` sample, Shannon 6.372 with the previous setting, 6.644 with `-g 1`).
- Updated pipeline ownership references from `nf-core` to `lescailab` metadata where applicable.
- `run_nanocdrx.py` writes every cluster representative to `<sample>_cdr3.tsv` exactly once, with the new statuses `cdr3-too-long`, `excluded-length` and `excluded-nonstandard`. Representatives containing non-standard residues are now excluded instead of having those residues deleted, which shortened their CDR3.
- `run_nanocdrx.py` calls `predict_cdrs` once per sample instead of once per 500 sequences.
- `NANOTRANSLATE` is faster and uses less memory: the start and end motifs are built from their amino acids instead of being listed one by one, each motif is searched once, reads are processed in groups instead of being loaded all at once, and the script uses the CPUs assigned by Nextflow instead of all those of the machine. The output is unchanged (verified byte by byte on two llama samples and an alpaca subset).
- Redesigned the interactive report page (`nanorepertoire_report.html`) without changing its data: each sample keeps one colour across all figures, from a palette checked for colour blindness (a single colour above eight samples, without legends); light and dark themes; the Top-10 clusters are one ranked list per sample with the share of its clustered sequences; the identity histograms share their bins; the CDR3 boundary offsets have one panel per terminus. Plotly.js is embedded in the page when the plotly package is installed, as in the module container, so the report opens without a network connection; `plotly=6.7.0` is added to the conda environment of `NANOREPERTOIRE_REPORT`.
- The interactive report cites every tool and method it names (Cutadapt, FLASH, nanoCDR-X, IMGT, Green Algorithms), no longer attributes its cluster size classes to Deschaght et al. 2017, and says where the execution and CO₂ reports are; `CITATIONS.md` lists nf-co2footprint and Green Algorithms.

### `Fixed`

- Fixed the CDR3 of large clusters being read on an atypical sequence: nanoCDR-X annotated the CD-HIT representative, which is the longest sequence of the cluster, and in large clusters the longest sequence is often an atypical one, for example with an extra stretch after framework 4. In the largest cluster of the llama `after` sample with the previous CD-HIT setting (121,504 sequences) the representative is 150 residues long, while 97.4% of the members are 128 residues long and 91.7% carry the CDR3 `AADQTPLRYSYSTRVDAYDY`; nanoCDR-X read `AADQTPLRYYYSTRVDAYDYWGQGTQVT` on the representative. The CDR3 is now annotated on the most frequent sequence of each cluster (`CLUSTER_MOSTFREQUENT`): on the two llama samples, with CD-HIT `-g 1 -det 1`, the clusters whose annotated sequence has a CDR3 that neither runs into framework 4 nor is missing or shorter than four residues hold 100% and 99.8% of the sequences, against 36.2% and 74.9% with the representative.
- Fixed `NANOCDRX` silently discarding whole 500-sequence batches: a single sequence longer than 150 residues, the input length of the nanoCDR-X model, made `predict_cdrs` fail on its batch, and the failure was ignored. The length filter is now 70-150 residues, and a failing or incomplete `predict_cdrs` run stops the step with an error.
- Fixed CDR3s longer than 50 residues being labelled `non-unique`.
- Fixed the nanoCDR-X version recorded in `versions.yml`, which was read from a `--version` option that `predict_cdrs` does not have.
- Fixed `NANOTRANSLATE` discarding merged reads that hold the VHH on the reverse strand, as half of the reads of a non-directional library do: when the start or end motif is missing on the forward strand, the reverse complement is searched too. The translation log reports how many reads were reverse-complemented.
- Fixed translations containing undetermined residues (`X`, from codons with an `N`) being kept: they are now discarded and counted in the translation log, like those with a stop codon. Reads with many `X` formed spurious singleton clusters.
- Fixed `FLASH` writing the merged reads in a different order on every run when it used more than one thread, which changed the read numbering of `RENAME` and the CD-HIT clusters between runs of the same data. FLASH now runs with `-t 1`, the setting its manual gives for deterministic output.
- Fixed the CD-HIT clusters depending on the number of threads and on the memory limit: with more than one thread CD-HIT assigns some sequences to different clusters depending on the thread count, and slightly differently on every run. `CDHIT_CDHIT` now uses CD-HIT 4.8.1 built from the reproducible branch of [lescailab/cdhit](https://github.com/lescailab/cdhit/tree/reproducible) (package `lescailab::cd-hit`), which adds the `-det` option, and runs with `-g 1 -det 1`: on the two llama samples the clusters are identical byte by byte with 1, 2, 4, 8 and 16 threads and with different memory limits, on macOS and on Linux. The `-T 0` and `-M 16000` options set in `conf/modules.config` had no effect, because the module appends its own `-T` and `-M`, and were removed.
- Fixed `NANOTRANSLATE` cutting the VHH at the first `TVSS` of the read when that motif also occurs inside the VHH (in a CDR or framework 3): the end is now the first in-frame `TVSS` that starts 7 residues after the tryptophan of `WG.G`, its position in the IMGT germline; when no `TVSS` is in that position the first one is kept, as before.
- Fixed `NANOTRANSLATE` discarding VHH whose first three residues or final `TVSS` differ by one amino acid from the start and end motifs (for example `QVR`, `QGQ`, `TVSA`, `IVSS`): reads with no exact motif on either strand are searched again in the six reading frames with anchors that tolerate one substitution, confirmed by residues conserved in the IMGT alpaca germline (the FR1 cysteine, the second cysteine and the `WG.G` motif of framework 4). The rescued reads are counted in the translation log and listed in a `.rescued` file next to `.notfound` and `.outframe`.
- Fixed an 11-nt end motif in `NANOTRANSLATE` (`CTGTTAGCTCA`, a truncated encoding of `TVSS`): it never matched its intended `ACTGTTAGCTCA` and could match elsewhere in a read, leaving the VHH region out of frame. All end motifs are now the 12-nt encodings of `TVSS`.
- The aggregated tables and the Quarto report list the samples sorted by name, instead of in the order in which the samples finished clustering.
- Fixed the nf-co2footprint reports being written in the launch directory: they are now written in `pipeline_info/`, next to the Nextflow execution reports, and listed in `docs/output.md`.
- Fixed the software versions report, which listed only FastQC and the reports: the versions of the `FASTQ_TO_FASTA` and `FASTA_CLUSTERING` steps were never collected. Also fixed the seqtk version recorded by `RENAME` and removed the pandas and plotly versions from `NANOREPERTOIRE_REPORT`, which does not use them.
- Improved Docker compatibility on Apple Silicon and fixed ARM64 Wave container handling for report rendering.
- Fixed resource-limit guard logic (`check_max`) and conda channel checks for Nextflow `25.10.4` + micromamba.
- Addressed nf-core compliance and linting issues across templates/modules/subworkflows.
- Fixed `AGGREGATE_STATS` writing `null` as the sample ID in `sampledata.tsv`: it read the key `sampleID`, which is never set (the sample name is `meta.id`), so the sample metadata never joined the clusters and individual, immunisation and boost were `NA` in every row of `cdr3_boost_overview_table.tsv`. Samplesheet columns left empty are now written as `NA` instead of `[]`. The alignment and tree step of the Quarto report, which the fix made run on every sample with metadata, now runs only with `--calculate_tree`, which already controlled the drawing of the trees.
- Fixed the methods text of both reports, which described FLASH with a 300 bp maximum overlap (the pipeline uses the FLASH defaults) and the translation with the old fixed motifs `ATGGCACAG` and `ACCGTCTCCTCA`; it now describes the current start and end motifs, the framework-4 end rule and the rescue with tolerant anchors. The Quarto report quotes the CD-HIT identity threshold and word size of the run instead of a fixed 0.9, and no longer prints debugging lines, a sentence left from development and two unfinished sections.
- Fixed the cluster size plot of the Quarto report, which drew only part of the clusters: its axes were fixed at 0-500 members and 0-5,000 clusters, and ggplot2 silently drops what falls outside, so the bar of the clusters with 1-8 members and every cluster above 500 members were never drawn (542 of 4,194 clusters drawn for the June llama `after` sample). The plot now shows, on logarithmic axes, how many clusters have each size, and draws every cluster.
- Fixed the Quarto report (`aggregate_stats/analysis_report.html`) showing no figures and no styling: the page linked them from an `analysis_report_files/` folder that the pipeline does not publish. The page now embeds its figures, styles and scripts (`embed-resources: true`).
- Fixed `NANOTRANSLATE` writing empty translations and fragments as translated VHH. When a read's first end motif came before its first start motif, the region between them was empty and an empty sequence was written; the end motif is now searched after the start, so such reads are tried on the other strand and with the tolerant anchors like any read without an end. Translations shorter than 70 residues (primer dimers such as `MAQVTVSS`, and fragments ending on a `TVSS` inside the VHH or starting at an `LQL` of framework 3) are now discarded and counted in the translation log. CD-HIT used to drop the empty and 8-residue translations, which is why a few translated sequences were missing from the clusters.
- Fixed the "Unique CDR3 Sequences – Distinct paratopes" card of the interactive report, which summed the distinct CDR3s of each sample and so counted a CDR3 found in several samples once per sample: it now counts the distinct CDR3s of all samples, and is labelled "Distinct CDR3 sequences". The card "Sequences with a called CDR3" counted cluster representatives, not sequences, and is now labelled "Representatives with a CDR3". The descriptions of the CDR3 panels say that they count each distinct CDR3 of a sample once, whatever the size of its cluster.
- Fixed the interpretation of cluster sizes in the reports. The interactive report described clusters with at least 5 members as B-cell lineages that had undergone antigen-driven somatic expansion, clusters with at least 1000 members as prime candidates for high-affinity nanobodies and the largest clusters as the most affinity-matured binders, and the Quarto report stated that it showed which lineages were affinity matured. The reports now list the possible causes of a large cluster (clonal expansion, a selection step such as phage display panning, transcript abundance, PCR bias, repeated sequencing of the same molecule), which they cannot tell apart, and state that cluster size does not measure binding affinity. Statements on CDR3 length classes, loop shapes and amino acid composition that no cited source supported are removed. The citation of Spinelli et al. 2022, whose DOI does not resolve, is removed, and the CDR3 length comparison now reports the values of Mitchell & Colwell 2018; the citations of Desmyter et al. 2002 and Muyldermans 2013, given with DOIs that do not resolve, are removed. The class of clusters with fewer than 5 members is labelled "1–4 members" instead of "Singletons (<5)".

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

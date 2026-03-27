<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-nanorepertoire_logo_dark.png">
    <img alt="nf-core/nanorepertoire" src="docs/images/nf-core-nanorepertoire_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-core/nanorepertoire)
[![GitHub Actions CI Status](https://github.com/nf-core/nanorepertoire/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanorepertoire/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17379842-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17379842)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.4.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/nanorepertoire)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nanorepertoire-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nanorepertoire)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction


**lescailab/nanorepertoire** is a Nextflow-based bioinformatics pipeline designed to characterise *nanobody repertoires* from raw sequencing data (FASTQ files).
It performs quality control, read preprocessing, translation, clustering, and CDR3 extraction, producing comprehensive reports that describe nanobody diversity and repertoire composition.

The workflow is divided into three main subworkflows:

### 1. `fastq_to_fasta`
Processes raw sequencing data to produce translated FASTA sequences:
- **FastQC** – quality control of raw reads
- **Cutadapt** – adapter trimming
- **FLASH** – paired-end read merging
- **Rename** – standardized renaming of merged reads
- **Nanotranslate** – translation from nucleotide to amino acid sequences

### 2. `fasta_clustering`
Clusters and analyses translated nanobody sequences:
- **CD-HIT** – clustering of identical or highly similar sequences (configurable identity threshold)
- **ReadCDHIT** – extraction and summarization of cluster statistics
- **getCDR3** – extraction of CDR3 regions from translated nanobodies
- **Nanocdr-x** – CDR3 extraction and analysis tool

### 3. `reporting`
Generates integrated reports and visual outputs:
- **MultiQC** – aggregation of QC metrics
- **Aggregate Stats (R-based)** – production of summary tables (`.csv`, `.tsv`), serialized objects (`.RData`), and a static analysis report.
- **Nanorepertoire Report (Python-based)** – generation of a modern, interactive HTML dashboard.

## Parameters

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `--cdhit_identity` | `0.9` | Sequence identity threshold for CD-HIT clustering. |
| `--cluster_size_threshold` | `1` | Minimum sequences per cluster for inclusion in reports. |
| `--calculate_tree` | `false` | Enable/disable phylogenetic tree calculation. |
| `--adapterfile` | `null` | Custom adapter file for Cutadapt. |

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

1. Prepare a samplesheet (`samplesheet.csv`):
```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

2. Run the pipeline:
```bash
nextflow run lescailab/nanorepertoire \
  -profile docker \
  --input samplesheet.csv \
  --outdir results/
```

## Pipeline output

The pipeline produces:
- **Translated Sequences**: Amino acid sequences in FASTA format.
- **Cluster Stats**: Detailed statistics on nanobody clusters.
- **QC Reports**: `multiqc_report.html` for technical quality metrics.
- **Hybrid Reports**:
    - `analysis_report.html`: Static scientific summary (R/Quarto).
    - `nanorepertoire_report.html`: Interactive, premium dashboard (Python/Plotly).
    - `nanobodies_report.RData`: Serialized R object for downstream analysis.

## Credits

**lescailab/nanorepertoire** was originally written by Francesco Lescai.
Key contributors: Davide Bagordo.

## Citations

If you use lescailab/nanorepertoire, please cite: [10.5281/zenodo.17379842](https://doi.org/10.5281/zenodo.17379842)

An extensive list of references for the tools used can be found in [`CITATIONS.md`](CITATIONS.md).
for proper execution.

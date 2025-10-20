<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-nanorepertoire_logo_dark.png">
    <img alt="nf-core/nanorepertoire" src="docs/images/nf-core-nanorepertoire_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-core/nanorepertoire)
[![GitHub Actions CI Status](https://github.com/nf-core/nanorepertoire/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanorepertoire/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
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
- **CD-HIT** – clustering of identical or highly similar sequences
- **ReadCDHIT** – extraction and summarization of cluster statistics
- **getCDR3** – extraction of CDR3 regions from translated nanobodies
- **MAFFT** – multiple sequence alignment of CDR3 clusters

### 3. `nanoreport`
Generates integrated reports and visual outputs:
- **MultiQC** – aggregation of QC metrics
- **Report (R-based)** – production of summary tables (`.tsv`), serialized objects (`.RData`), and an HTML report summarizing repertoire metrics

A visual representation of the three main subworkflows:


<p align="center">
  <img src="docs/images/nanorepertoire_workflow.png" width="800" alt="Nanorepertoire workflow overview">
</p>

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.


First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).



Now, you can run the pipeline using:

**Example command**
```bash
nextflow run lescailab/nanorepertoire \
  -profile docker \
  --input samplesheet.csv \
  --outdir results/
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/nanorepertoire/usage) and the [parameter documentation](https://nf-co.re/nanorepertoire/parameters).

## Pipeline output

The pipeline produces:

- Translated nanobody FASTA sequences

- Clustered and aligned CDR3 repertoires

- Cluster statistics and summary tables

- Quality control reports (multiqc.html)

- Final integrated report (report.html, .RData, .tsv)


## Credits

lescailab/nanorepertoire was originally written by Francesco Lescai.

We thank the following people for their extensive assistance in the development of this pipeline:
- Davide Bagordo


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#nanorepertoire` channel](https://nfcore.slack.com/channels/nanorepertoire) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations


If you use lescailab/nanorepertoire for your analysis, please cite it using the following doi: [10.5281/zenodo.17379842](https://doi.org/10.5281/zenodo.17379842)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


## Limitations

It is important to note that the pipeline currently does not pass all GitHub CI tests due to incorrect container versioning affecting some modules.

We recommend using Nextflow version 24.10.4 build 5934 for proper execution.

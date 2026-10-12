#!/usr/bin/env python3
"""Count the sequences of each sample at every step of the pipeline.

The counts come from files the pipeline already writes:

* ``<sample>.cutadapt.log``: read pairs given to and written by Cutadapt;
* ``<sample>_AA_translated.fasta.log``: reads merged by FLASH (the reads the
  translation step receives), translated VHH, rescued reads and the reason of
  every failure;
* the cluster sizes of ``clusterbig.csv``: sequences in clusters and number of
  clusters, i.e. of cluster representatives;
* ``<sample>_cdr3_summary.tsv``: representatives excluded before nanoCDR-X, sent
  to it, with and without a CDR3, and with a distinct CDR3.

Every step is followed by the losses that separate it from the next one, so
that each row of the table can be checked against its neighbours. A count that
cannot be read is ``None`` (``NA`` in the outputs), never 0.
"""

import csv
import os
import re

CUTADAPT_SUFFIX = ".cutadapt.log"
TRANSLATION_SUFFIX = "_AA_translated.fasta.log"

# Label prefixes of the translation log (bin/translate.py): each label line is
# followed by a line holding the count.
TRANSLATION_LABELS = [
    ("Completed reading reads", "merged"),
    ("Of which passing", "translated"),
    ("of which rescued", "rescued"),
    ("failed because missing start", "no_motif"),
    ("failed because not in frame", "out_of_frame"),
    ("failed because there is a stop codon", "stop_codon"),
    ("failed because the translation contains an undetermined", "undetermined"),
]
TRANSLATION_FAILURES = ["no_motif", "out_of_frame", "stop_codon", "undetermined"]

# key, description, unit, kind ("step", "loss" or "note")
ROWS = [
    ("read_pairs", "Read pairs sequenced", "read pairs", "step"),
    ("removed_cutadapt", "removed by Cutadapt", "read pairs", "loss"),
    ("read_pairs_trimmed", "Read pairs after Cutadapt", "read pairs", "step"),
    ("not_merged", "not merged by FLASH", "read pairs", "loss"),
    ("merged", "Merged reads (FLASH)", "reads", "step"),
    ("no_motif", "no VHH start or end motif", "reads", "loss"),
    ("out_of_frame", "VHH region not a multiple of three nucleotides", "reads", "loss"),
    ("stop_codon", "stop codon in the VHH region", "reads", "loss"),
    ("undetermined", "undetermined residue (X) in the translation", "reads", "loss"),
    ("translated", "Translated VHH sequences", "sequences", "step"),
    ("rescued", "of which rescued with tolerant anchors", "sequences", "note"),
    ("not_clustered", "translated but missing from the clusters", "sequences", "loss"),
    ("clustered", "Sequences in clusters (CD-HIT)", "sequences", "step"),
    ("clusters", "Clusters, one representative each (CD-HIT)", "representatives", "step"),
    ("excluded_length", "outside the 70-150 residues accepted by nanoCDR-X", "representatives", "loss"),
    ("excluded_nonstandard", "non-standard residue, not accepted by nanoCDR-X", "representatives", "loss"),
    ("sent_to_model", "Representatives annotated by nanoCDR-X", "representatives", "step"),
    ("no_cdr3", "no CDR3 called", "representatives", "loss"),
    ("with_cdr3", "Representatives with a CDR3", "representatives", "step"),
    ("non_unique", "same CDR3 as another representative of the sample", "representatives", "loss"),
    ("cdr3_too_long", "CDR3 longer than 50 residues", "representatives", "loss"),
    ("unique", "Distinct CDR3 sequences", "representatives", "step"),
]
# First row of each block: percentages of the steps are relative to it.
BLOCK_START = {"read_pairs": "read_pairs", "clusters": "clusters"}


def _int(text):
    return int(text.replace(",", "").strip())


def sample_name(path, suffix):
    """Sample name of a per-sample file, i.e. its file name without the suffix."""
    name = os.path.basename(path)
    return name[: -len(suffix)] if name.endswith(suffix) else None


def parse_cutadapt_log(path):
    """Return {"read_pairs", "read_pairs_trimmed"} from a Cutadapt report.

    Paired-end reports say "read pairs" and "Pairs written", single-end ones
    "reads" and "Reads written"; both are accepted.
    """
    with open(path) as handle:
        text = handle.read()
    counts = {"read_pairs": None, "read_pairs_trimmed": None}
    found = re.search(r"^Total (?:read pairs|reads) processed:\s+([\d,]+)", text, re.M)
    if found:
        counts["read_pairs"] = _int(found.group(1))
    found = re.search(r"^(?:Pairs|Reads) written \(passing filters\):\s+([\d,]+)", text, re.M)
    if found:
        counts["read_pairs_trimmed"] = _int(found.group(1))
    return counts


def parse_translation_log(path):
    """Return the counts of the translation log; a label missing from the log is None."""
    counts = {key: None for _, key in TRANSLATION_LABELS}
    with open(path) as handle:
        lines = [line.strip() for line in handle]
    for label_line, value_line in zip(lines, lines[1:]):
        for prefix, key in TRANSLATION_LABELS:
            if label_line.startswith(prefix) and value_line.isdigit():
                counts[key] = int(value_line)
    return counts


def parse_cdr3_summary(path):
    """Return {sample: counts} from a <sample>_cdr3_summary.tsv file."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["sample"]: {k: _int(v) for k, v in row.items() if k != "sample"} for row in rows}


def _difference(a, b):
    return None if a is None or b is None else a - b


def account(cutadapt=None, translation=None, cluster_sizes=None, cdr3=None):
    """Return ({row key: count}, [warnings]) for one sample."""
    cutadapt = cutadapt or {}
    translation = translation or {}
    cdr3 = cdr3 or {}
    sizes = cluster_sizes
    counts = {key: None for key, _, _, _ in ROWS}

    counts["read_pairs"] = cutadapt.get("read_pairs")
    counts["read_pairs_trimmed"] = cutadapt.get("read_pairs_trimmed")
    for key in ["merged", "translated", "rescued"] + TRANSLATION_FAILURES:
        counts[key] = translation.get(key)
    if sizes is not None:
        counts["clustered"] = sum(sizes)
        counts["clusters"] = len(sizes)
    for key in ["excluded_length", "excluded_nonstandard", "sent_to_model", "no_cdr3",
                "with_cdr3", "non_unique", "cdr3_too_long", "unique"]:
        counts[key] = cdr3.get(key)

    counts["removed_cutadapt"] = _difference(counts["read_pairs"], counts["read_pairs_trimmed"])
    counts["not_merged"] = _difference(counts["read_pairs_trimmed"], counts["merged"])
    counts["not_clustered"] = _difference(counts["translated"], counts["clustered"])

    warnings = []
    failures = [counts[k] for k in TRANSLATION_FAILURES if counts[k] is not None]
    if counts["merged"] is not None and counts["translated"] is not None and \
            counts["merged"] != counts["translated"] + sum(failures):
        warnings.append(
            f"the translation log does not add up: {counts['merged']:,} reads read, "
            f"{counts['translated']:,} translated and {sum(failures):,} failed")
    for key in ("removed_cutadapt", "not_merged", "not_clustered"):
        if counts[key] is not None and counts[key] < 0:
            description = next(d for k, d, _, _ in ROWS if k == key)
            warnings.append(f"negative count for \"{description}\" ({counts[key]:,})")
    if cdr3:
        if counts["clusters"] is not None and cdr3.get("input_sequences") != counts["clusters"]:
            warnings.append(
                f"nanoCDR-X received {cdr3.get('input_sequences'):,} representatives, "
                f"but CD-HIT produced {counts['clusters']:,} clusters")
        checks = [
            ("input_sequences", ["excluded_length", "excluded_nonstandard", "sent_to_model"]),
            ("sent_to_model", ["no_cdr3", "with_cdr3"]),
            ("with_cdr3", ["unique", "non_unique", "cdr3_too_long"]),
        ]
        for total, parts in checks:
            if cdr3.get(total) != sum(cdr3.get(p, 0) for p in parts):
                warnings.append(f"the CDR3 summary does not add up: {total} is not "
                                f"the sum of {', '.join(parts)}")
    return counts, warnings


def percentage(counts, key):
    """Percentage of a step relative to the first row of its block, or None."""
    start = None
    for row_key, _, _, kind in ROWS:
        if row_key in BLOCK_START:
            start = row_key
        if row_key == key:
            break
    base, value = counts.get(start), counts.get(key)
    if not base or value is None:
        return None
    return 100.0 * value / base


def write_tsv(table, path):
    """Write one line per sample and row: Sample, Step, Description, Unit, Count."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample", "Step", "Description", "Unit", "Count"])
        for sample in sorted(table):
            counts = table[sample]
            for key, description, unit, _ in ROWS:
                value = counts.get(key)
                writer.writerow([sample, key, description, unit, "NA" if value is None else value])

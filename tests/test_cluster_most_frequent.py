"""Unit tests for ``bin/cluster_most_frequent.py``.

Run with ``pytest tests/test_cluster_most_frequent.py`` from the repository root.
"""

import csv
import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA = REPO_ROOT / "modules" / "local" / "clustermostfrequent" / "tests" / "data"


def _load():
    spec = importlib.util.spec_from_file_location(
        "cluster_most_frequent", REPO_ROOT / "bin" / "cluster_most_frequent.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cmf = _load()


@pytest.fixture(scope="module")
def outputs(tmp_path_factory):
    prefix = tmp_path_factory.mktemp("cmf") / "test"
    cmf.main([str(DATA / "test_sequences.fasta"), str(DATA / "test_clusters.clstr"), str(prefix)])
    with open(str(prefix) + ".tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with open(str(prefix) + ".fasta") as handle:
        fasta = cmf.read_fasta(handle)
    with open(DATA / "test_sequences.fasta") as handle:
        seqs = cmf.read_fasta(handle)
    return rows, fasta, seqs


def test_one_row_and_one_sequence_per_cluster(outputs):
    rows, fasta, _ = outputs
    assert [r["cluster_representative"] for r in rows] == ["1", "6", "7"]
    assert list(fasta) == ["1", "6", "7"]


def test_most_frequent_replaces_longest(outputs):
    rows, fasta, seqs = outputs
    assert rows[0]["chosen_member"] == "2"
    assert rows[0]["chosen_copies"] == "3"
    assert rows[0]["same_as_representative"] == "no"
    assert fasta["1"] == seqs["2"]


def test_singleton_keeps_its_sequence(outputs):
    rows, fasta, seqs = outputs
    assert rows[1]["chosen_member"] == "6"
    assert fasta["6"] == seqs["6"]


def test_tie_goes_to_first_member(outputs):
    rows, fasta, seqs = outputs
    assert rows[2]["chosen_member"] == "7"
    assert rows[2]["chosen_copies"] == "2"
    assert rows[2]["same_as_representative"] == "yes"


def test_missing_member_is_an_error(tmp_path):
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">1\nMAQ\n")
    with pytest.raises(ValueError):
        cmf.main([str(seqs), str(DATA / "test_clusters.clstr"), str(tmp_path / "out")])

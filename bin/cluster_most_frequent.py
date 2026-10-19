#!/usr/bin/env python
"""Choose, for every CD-HIT cluster, the sequence on which its CDR3 is annotated.

CD-HIT makes the longest sequence of a cluster its representative. In large clusters the longest
sequence is often an atypical one (an extra stretch after framework 4, a repeated start), and the
CDR3 read on it is wrong or missing. This script picks instead the most frequent sequence of the
cluster: the members are grouped by identical sequence and the largest group wins; ties go to the
member CD-HIT lists first in the cluster (CD-HIT lists members by decreasing length, then in input
order), so the choice is deterministic.

Outputs, for ``PREFIX``:

``PREFIX.fasta``  one sequence per cluster. The header is the identifier of the CD-HIT
                  representative, which is the key the rest of the pipeline uses for the cluster;
                  the sequence is the most frequent sequence of the cluster.
``PREFIX.tsv``    one row per cluster, tab-separated::

    cluster_representative  chosen_member  cluster_size  chosen_copies  chosen_share
    representative_length   chosen_length  same_as_representative

``chosen_copies`` is the number of members identical to the chosen sequence and ``chosen_share``
that number divided by ``cluster_size``. ``same_as_representative`` is ``yes`` when the chosen
sequence is identical to the CD-HIT representative.

Usage: cluster_most_frequent.py SEQUENCES.fasta CLUSTERS.clstr PREFIX
"""

import argparse
import sys
from collections import Counter


def read_fasta(handle):
    """Return {identifier: sequence}. Headers may start with '>' or '@', as CD-HIT accepts both."""
    seqs, name, buf = {}, None, []
    for line in handle:
        line = line.rstrip("\n")
        if line[:1] in (">", "@"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def read_clstr(handle):
    """Yield (representative, [members in CD-HIT order]) for every cluster."""
    members, rep = [], None
    for lineno, line in enumerate(handle, start=1):
        if not line.strip():
            continue
        if line.startswith(">Cluster"):
            if members:
                if rep is None:
                    raise ValueError("cluster ending at line {0} has no representative ('*')".format(lineno - 1))
                yield rep, members
            members, rep = [], None
            continue
        if ">" not in line:
            raise ValueError("line {0} is neither a cluster header nor a member line: {1!r}".format(lineno, line))
        sid = line.split(">", 1)[1].split("...", 1)[0]
        members.append(sid)
        if line.rstrip().endswith("*"):
            rep = sid
    if members:
        if rep is None:
            raise ValueError("last cluster has no representative ('*')")
        yield rep, members


def choose(members, seqs):
    """Return (chosen member, copies of its sequence) for one cluster."""
    counts = Counter(seqs[m] for m in members)
    top = max(counts.values())
    chosen = next(m for m in members if counts[seqs[m]] == top)
    return chosen, top


def main(argv=None):
    parser = argparse.ArgumentParser(description="Choose the most frequent sequence of every CD-HIT cluster.")
    parser.add_argument("sequences", help="FASTA file clustered by CD-HIT")
    parser.add_argument("clusters", help="CD-HIT cluster file (*.clstr)")
    parser.add_argument("prefix", help="prefix of the output files")
    args = parser.parse_args(argv)

    with open(args.sequences) as handle:
        seqs = read_fasta(handle)
    with open(args.clusters) as cl, open(args.prefix + ".fasta", "w") as fa, open(args.prefix + ".tsv", "w") as tsv:
        tsv.write("cluster_representative\tchosen_member\tcluster_size\tchosen_copies\tchosen_share\t"
                  "representative_length\tchosen_length\tsame_as_representative\n")
        for rep, members in read_clstr(cl):
            missing = [m for m in members if m not in seqs]
            if missing:
                raise ValueError("cluster of {0}: {1} member(s) not in {2}, e.g. {3}".format(
                    rep, len(missing), args.sequences, missing[0]))
            chosen, copies = choose(members, seqs)
            fa.write(">{0}\n{1}\n".format(rep, seqs[chosen]))
            tsv.write("{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5}\t{6}\t{7}\n".format(
                rep, chosen, len(members), copies, copies / len(members), len(seqs[rep]),
                len(seqs[chosen]), "yes" if seqs[chosen] == seqs[rep] else "no"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
# Translates the VHH region of merged reads: the region runs from a start motif to the next end
# motif, searched on both strands, and is kept only if in frame, without stop codons, without
# undetermined residues and at least MIN_VHH_LENGTH residues long. When the read holds more than
# one TVSS, the end is the one 7 residues after the tryptophan of WG.G. Reads with no exact motif
# on either strand are searched again with anchors that tolerate one amino-acid substitution (see
# rescue()). Reads are processed in parallel, in input order.
#
# Usage: translate.py <reads.fastq> <output.fasta> [cpus]

import sys
import re
import multiprocessing
from itertools import islice
from Bio.Seq import translate
from multiprocessing import Pool

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# DNA codons of the amino acids used in the motifs (standard genetic code)
CODON = {
    'A': 'GC[ACGT]', 'E': 'GA[AG]', 'L': '(?:CT[ACGT]|TT[AG])', 'M': 'ATG',
    'Q': 'CA[AG]', 'S': '(?:TC[ACGT]|AG[CT])', 'T': 'AC[ACGT]', 'V': 'GT[ACGT]',
}


def motif(peptide):
    """Regular expression matching every DNA encoding of peptide."""
    return ''.join(CODON[aa] for aa in peptide)


# VHH start: every DNA encoding of MAQ, QVQ, EVQ, VQL or LQL (9 nt)
START = re.compile('|'.join(motif(peptide) for peptide in ('MAQ', 'QVQ', 'EVQ', 'VQL', 'LQL')))
# VHH end: every DNA encoding of TVSS (12 nt)
END = re.compile(motif('TVSS'))

# Shortest translation kept as a VHH. A complete VHH, from the framework-1 start to the framework-4
# TVSS, is longer than 100 residues; shorter translations are primer dimers (MAQ followed directly by
# TVSS) or fragments cut by a TVSS inside the VHH, or starting at an LQL of framework 3. 70 is also
# the shortest sequence accepted by nanoCDR-X.
MIN_VHH_LENGTH = 70

# Rescue of reads in which no exact start or end motif is found on either strand. The VHH is
# looked for in the six reading frames with anchors that tolerate one amino-acid substitution,
# each confirmed by residues conserved in the IMGT germline of the alpaca:
# - start: a tripeptide at most one residue away from a germline FR1 start (QVQ, EVQ, QLQ),
#   with the FR1 cysteine (1st-CYS, IMGT position 23) 21 residues later;
# - second cysteine (2nd-CYS, IMGT position 104): a C 88 to 100 residues after the start (94 to
#   96 in the germline; the window allows for CDR1 and CDR2 length differences);
# - end: the first WG.G (J-TRP, IMGT position 118) at least 80 residues after the start,
#   followed 7 residues later by a tetrapeptide at most one residue away from TVSS.
RESCUE_STARTS = ('QVQ', 'EVQ', 'QLQ')
FR4_ANCHOR = re.compile(r'WG.G')


def near(peptide, targets):
    """True when peptide differs in at most one residue from one of targets of the same length."""
    return any(len(peptide) == len(target) and sum(a != b for a, b in zip(peptide, target)) <= 1
               for target in targets)


def rescue(dna):
    """Return (vhh, strand) of a VHH found with the tolerant anchors, or None."""
    for strand, seq in (('forward', dna), ('reverse', dna.translate(COMPLEMENT)[::-1])):
        for frame in range(3):
            coding = seq[frame:]
            protein = translate(coding[:len(coding) // 3 * 3])
            for start in range(len(protein) - 21):
                if protein[start + 21] != 'C' or not near(protein[start:start + 3], RESCUE_STARTS):
                    continue
                if 'C' not in protein[start + 88:start + 101]:
                    continue
                trp = FR4_ANCHOR.search(protein, start + 80)
                if not trp or not near(protein[trp.start() + 7:trp.start() + 11], ('TVSS',)):
                    continue
                vhh = protein[start:trp.start() + 11]
                if '*' not in vhh and 'X' not in vhh:
                    return vhh, strand
    return None


END_ANY = re.compile('(?=' + motif('TVSS') + ')')


def framework4_end(dna, startbase, endbase):
    """Position of the first TVSS in frame with the start whose residue 7 before it starts WG.G.

    In the IMGT germline TVSS starts 7 residues after the tryptophan of WG.G (J-TRP, IMGT 118);
    an earlier TVSS inside the VHH would cut the protein. When no TVSS is in that position,
    the first end motif of the read is kept (endbase).
    """
    for match in END_ANY.finditer(dna, startbase):
        end = match.start()
        if (end - startbase) % 3 == 0 and end - 21 >= startbase and FR4_ANCHOR.fullmatch(translate(dna[end - 21:end - 9])):
            return end
    return endbase


def read_fastq(file_handle):
    """Generator to read FASTQ file."""
    while True:
        headerline = file_handle.readline()
        if not headerline:
            break
        sequenceline = file_handle.readline()
        plusline = file_handle.readline()
        qualityline = file_handle.readline()
        yield (headerline.rstrip(), sequenceline.rstrip(), plusline.rstrip(), qualityline.rstrip())

def find_motifs(dna):
    """Return the positions of the first start motif and of the first end motif after it, -1 when absent.

    An end motif before the start motif delimits no region: the read is treated as having no end, so
    that the other strand and the tolerant anchors are tried. Without a start, the first end of the
    read is returned, for the count of reads with an end motif.
    """
    start = START.search(dna)
    end = END.search(dna, start.end()) if start else END.search(dna)
    return (start.start() if start else -1), (end.start() if end else -1)

def process_read(read):
    """Processes a single read."""
    headerline, sequenceline, _, _ = read
    result = {'reversed': 0}
    dna = sequenceline.rstrip()

    startbase, endbase = find_motifs(dna)
    if startbase == -1 or endbase == -1:
        # In non-directional libraries half of the merged reads hold the VHH on the other strand
        reverse = dna.translate(COMPLEMENT)[::-1]
        reverse_start, reverse_end = find_motifs(reverse)
        if reverse_start != -1 and reverse_end != -1:
            dna, startbase, endbase = reverse, reverse_start, reverse_end
            result['reversed'] = 1
    result['foundstart'] = int(startbase != -1)
    result['foundend'] = int(endbase != -1)

    if startbase == -1 or endbase == -1:
        rescued = rescue(dna)
        if rescued:
            vhh, strand = rescued
            result['status'] = 'rescued'
            result['header'] = headerline
            result['translated'] = vhh + "\n"
            result['rescued'] = f"{headerline}\t{strand}\t{vhh[:3]}\t{vhh[-4:]}\t{len(vhh)}\t{vhh}\n"
            return result
        result['status'] = 'nostartnoend'
        result['notfound'] = f"sequenza non trovata\t{dna}\n"
        return result

    endbase = framework4_end(dna, startbase, endbase)
    targetregion = dna[startbase:endbase + 12]  # 12 is the length of the end sequence
    if len(targetregion) % 3 != 0:
        result['status'] = 'notinframe'
        result['nomultiple'] = f"{headerline}\t{dna}\t{startbase}\t{endbase}\t{targetregion}\n"
        return result

    # Translate the target region
    translated = translate(targetregion)
    if '*' in translated:
        result['status'] = 'withstop'
        return result
    if 'X' in translated:
        result['status'] = 'withundetermined'
        return result
    if len(translated) < MIN_VHH_LENGTH:
        result['status'] = 'tooshort'
        return result

    result['status'] = 'passed'
    result['header'] = headerline
    result['translated'] = translated + "\n"
    return result

def process_reads(reads):
    """Processes a batch of reads."""
    return [process_read(read) for read in reads]

def batches(reads, size):
    """Group reads into lists of size reads."""
    while True:
        batch = list(islice(reads, size))
        if not batch:
            return
        yield batch

def main():
    cpus = int(sys.argv[3]) if len(sys.argv) > 3 else multiprocessing.cpu_count()
    with open(sys.argv[1], 'r') as filein, \
         open(sys.argv[2], 'w') as fileout, \
         open(sys.argv[2] + '.notfound', 'w') as notfound, \
         open(sys.argv[2] + '.log', 'w') as log, \
         open(sys.argv[2] + '.outframe', 'w') as nomultiplefile, \
         open(sys.argv[2] + '.rescued', 'w') as rescuedfile:

        readsread = 0
        readspassing = 0
        foundstart = 0
        foundend = 0
        notinframe = 0
        withstop = 0
        withundetermined = 0
        tooshort = 0
        reversed_reads = 0
        nostartnoend = 0
        rescued_reads = 0

        # Batches of 1000 reads are processed in parallel, a few per CPU at a time, so that
        # only part of the file is held in memory; results are written in input order
        with Pool(processes=cpus) as pool:
            for group in batches(batches(read_fastq(filein), 1000), cpus * 4):
                for results in pool.map(process_reads, group):
                    for result in results:
                        readsread += 1
                        if result['status'] == 'nostartnoend':
                            nostartnoend += 1
                            notfound.write(result['notfound'])
                        elif result['status'] == 'notinframe':
                            notinframe += 1
                            nomultiplefile.write(result['nomultiple'])
                        elif result['status'] == 'withstop':
                            withstop += 1
                        elif result['status'] == 'withundetermined':
                            withundetermined += 1
                        elif result['status'] == 'tooshort':
                            tooshort += 1
                        elif result['status'] == 'rescued':
                            readspassing += 1
                            rescued_reads += 1
                            rescuedfile.write(result['rescued'])
                            fileout.write(result['header'] + '\n')
                            fileout.write(result['translated'])
                        elif result['status'] == 'passed':
                            readspassing += 1
                            fileout.write(result['header'] + '\n')
                            fileout.write(result['translated'])
                        foundstart += result['foundstart']
                        foundend += result['foundend']
                        reversed_reads += result['reversed']

        # Write log file
        log.write("Completed reading reads:\n")
        log.write(f"{readsread}\n")
        log.write("Of which passing and written into FASTA:\n")
        log.write(f"{readspassing}\n")
        log.write("of which rescued with tolerant anchors (no exact start or end motif)\n")
        log.write(f"{rescued_reads}\n")
        log.write("start matched are\n")
        log.write(f"{foundstart}\n")
        log.write("end matched are\n")
        log.write(f"{foundend}\n")
        log.write("failed because missing start or missing end, also with tolerant anchors\n")
        log.write(f"{nostartnoend}\n")
        log.write("failed because not in frame, i.e. not multiple of 3\n")
        log.write(f"{notinframe}\n")
        log.write("failed because there is a stop codon in the sequence\n")
        log.write(f"{withstop}\n")
        log.write("failed because the translation contains an undetermined residue (X)\n")
        log.write(f"{withundetermined}\n")
        log.write(f"failed because the translation is shorter than {MIN_VHH_LENGTH} residues\n")
        log.write(f"{tooshort}\n")
        log.write("located on the reverse strand (reverse-complemented before translation)\n")
        log.write(f"{reversed_reads}\n")

if __name__ == "__main__":
    main()

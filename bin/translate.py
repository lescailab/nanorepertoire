#!/usr/bin/env python
# Translates the VHH region of merged reads: the region runs from a start motif to an end
# motif, searched on both strands, and is kept only if in frame, without stop codons and
# without undetermined residues. Reads are processed in parallel, in input order.
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
# VHH end: every DNA encoding of TVSS (12 nt), except that ACTGTTAGCTCA is matched as
# CTGTTAGCTCA (11 nt), as in the earlier list of motifs, so that the output is unchanged
END = re.compile('(?!ACTGTTAGCTCA)' + motif('TVSS') + '|CTGTTAGCTCA')

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
    """Return the positions of the first start motif and of the first end motif in dna, -1 when absent."""
    start = START.search(dna)
    end = END.search(dna)
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
        result['status'] = 'nostartnoend'
        result['notfound'] = f"sequenza non trovata\t{dna}\n"
        return result

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
         open(sys.argv[2] + '.outframe', 'w') as nomultiplefile:

        readsread = 0
        readspassing = 0
        foundstart = 0
        foundend = 0
        notinframe = 0
        withstop = 0
        withundetermined = 0
        reversed_reads = 0
        nostartnoend = 0

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
        log.write("start matched are\n")
        log.write(f"{foundstart}\n")
        log.write("end matched are\n")
        log.write(f"{foundend}\n")
        log.write("failed because missing start or missing end\n")
        log.write(f"{nostartnoend}\n")
        log.write("failed because not in frame, i.e. not multiple of 3\n")
        log.write(f"{notinframe}\n")
        log.write("failed because there is a stop codon in the sequence\n")
        log.write(f"{withstop}\n")
        log.write("failed because the translation contains an undetermined residue (X)\n")
        log.write(f"{withundetermined}\n")
        log.write("located on the reverse strand (reverse-complemented before translation)\n")
        log.write(f"{reversed_reads}\n")

if __name__ == "__main__":
    main()

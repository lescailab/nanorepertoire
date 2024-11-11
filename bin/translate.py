#!/usr/bin/env python
#This is a script to take input fastq data from FLASH and translate it with biopython
#Input data is in files like:
#/usr/share/sequencing/projects/317/alignments/317-D1-Bst3_S1_L001.relabel.fastq
#Then the translation needs to run from ATGGCACAG (MAQ) up to ACCGTCTCCTCA (TVSS)
#Note that most merged reads start as: CCGGCCATGGCACAG

#### NB THIS SCRIPT ONLY WORKS WITH FASTQ HEADER SHORT

import Bio
from Bio.Seq import Seq
import sys
import gzip
import re

filein = open(sys.argv[1])
fileout = open(sys.argv[2], "w")
log = open(str(sys.argv[2] + ".log"), "w")
nomultiplefile = open(str(sys.argv[2] + ".outframe"), "w")

# Regex pattern for matching sequences
startseq = r'(ATGGCTCAA|ATGGCTCAG|ATGGCCCAA|ATGGCCCAG|ATGGCACAA|ATGGCACAG|ATGGCGCAA|ATGGCGCAG|' \
    r'GTTCAATTA|GTTCAATTG|GTTCAACTT|GTTCAACTC|GTTCAACTA|GTTCAACTG|GTTCAGTTA|GTTCAGTTG|' \
    r'GTTCAGCTT|GTTCAGCTC|GTTCAGCTA|GTTCAGCTG|GTCCAATTA|GTCCAATTG|GTCCAACTT|GTCCAACTC|' \
    r'GTCCAACTA|GTCCAACTG|GTCCAGTTA|GTCCAGTTG|GTCCAGCTT|GTCCAGCTC|GTCCAGCTA|GTCCAGCTG|' \
    r'GTACAATTA|GTACAATTG|GTACAACTT|GTACAACTC|GTACAACTA|GTACAACTG|GTACAGTTA|GTACAGTTG|' \
    r'GTACAGCTT|GTACAGCTC|GTACAGCTA|GTACAGCTG|GTGCAATTA|GTGCAATTG|GTGCAACTT|GTGCAACTC|' \
    r'GTGCAACTA|GTGCAACTG|GTGCAGTTA|GTGCAGTTG|GTGCAGCTT|GTGCAGCTC|GTGCAGCTA|GTGCAGCTG|' \
    r'GAAGTTCAA|GAAGTTCAG|GAAGTCCAA|GAAGTCCAG|GAAGTACAA|GAAGTACAG|GAAGTGCAA|GAAGTGCAG|' \
    r'GAGGTTCAA|GAGGTTCAG|GAGGTCCAA|GAGGTCCAG|GAGGTACAA|GAGGTACAG|GAGGTGCAA|GAGGTGCAG|' \
    r'TTACAATTA|TTACAATTG|TTACAACTT|TTACAACTC|TTACAACTA|TTACAACTG|TTACAGTTA|TTACAGTTG|' \
    r'TTACAGCTT|TTACAGCTC|TTACAGCTA|TTACAGCTG|TTGCAATTA|TTGCAATTG|TTGCAACTT|TTGCAACTC|' \
    r'TTGCAACTA|TTGCAACTG|TTGCAGTTA|TTGCAGTTG|TTGCAGCTT|TTGCAGCTC|TTGCAGCTA|TTGCAGCTG|' \
    r'CTTCAATTA|CTTCAATTG|CTTCAACTT|CTTCAACTC|CTTCAACTA|CTTCAACTG|CTTCAGTTA|CTTCAGTTG|' \
    r'CTTCAGCTT|CTTCAGCTC|CTTCAGCTA|CTTCAGCTG|CTCCAATTA|CTCCAATTG|CTCCAACTT|CTCCAACTC|' \
    r'CTCCAACTA|CTCCAACTG|CTCCAGTTA|CTCCAGTTG|CTCCAGCTT|CTCCAGCTC|CTCCAGCTA|CTCCAGCTG|' \
    r'CTACAATTA|CTACAATTG|CTACAACTT|CTACAACTC|CTACAACTA|CTACAACTG|CTACAGTTA|CTACAGTTG|' \
    r'CTACAGCTT|CTACAGCTC|CTACAGCTA|CTACAGCTG|CTGCAATTA|CTGCAATTG|CTGCAACTT|CTGCAACTC|' \
    r'CTGCAACTA|CTGCAACTG|CTGCAGTTA|CTGCAGTTG|CTGCAGCTT|CTGCAGCTC|CTGCAGCTA|CTGCAGCTG|' \
    r'CAAGTTCAA|CAAGTTCAG|CAAGTCCAA|CAAGTCCAG|CAAGTACAA|CAAGTACAG|CAAGTGCAA|CAAGTGCAG|' \
    r'CAGGTTCAA|CAGGTTCAG|CAGGTCCAA|CAGGTCCAG|CAGGTACAA|CAGGTACAG|CAGGTGCAA|CAGGTGCAG)'
endseq = "ACCGTCTCCTCA"

readsread = 0
readspassing = 0
foundstart = 0
foundend = 0
notinframe = 0
withstop = 0
nostartnoend = 0

linecounter = 0
for line in filein:
    linecounter += 1
    if linecounter % 4 == 1:
        readsread += 1
        headerline = line
    if linecounter % 4 == 2:
        dna = Seq(line.rstrip())
        if len(re.findall(startseq, str(dna))) > 0:
            startbase = dna.find(re.findall(startseq, str(dna))[0])
            foundstart += 1
        else:
            startbase = -1
        endbase = dna.find(endseq)
        if endbase > -1:
            foundend += 1
        if startbase == -1 or endbase == -1: #read fails if missing the start or end sequence
            nostartnoend += 1
            log.write("sequenza non trovata" + "\t" + str(dna) + "\n")
            continue
        targetregion = dna[startbase:endbase+len(endseq)]
        if len(str(targetregion)) % 3 != 0: #read fails if the grabbed region is not a multiple of three
            notinframe += 1
            nomultiplefile.write(headerline + "\t" + str(dna) + "\t" + str(startbase) + "\t" + str(endbase) + "\t" + str(targetregion) + "\n")
            continue
        translated = targetregion.translate()
        if "*" in str(translated): #read fails if it has a stop codon
            withstop += 1
            continue
        readspassing += 1
        fileout.write(headerline)
        fileout.write(str(translated)+"\n")

print("Completed reading reads:")
print(readsread)
print("Of which passing and written into FASTA:")
print(readspassing)
print("start matched are")
print(foundstart)
print("end matched are")
print(foundend)
print("failed because missing start or missing end")
print(nostartnoend)
print("failed because not in frame, i.e. not multiple of 3")
print(notinframe)
print("failed because there is a stop codon in the sequence")
print(withstop)

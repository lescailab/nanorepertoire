#!/usr/bin/env python3
import sys
import csv
import subprocess
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--prefix", required=True)
    args, unknown = parser.parse_known_args()

    input_file = args.input
    prefix = args.prefix

    # 1. Convert Input to CSV for nanocdrx
    sequences = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            if i+1 >= len(lines): break
            header = lines[i].strip()
            seq = lines[i+1].strip()
            # Clean header as getcdr3 does
            identifier = header.replace("@", "").replace(">", "")
            sequences.append({'identifier': identifier, 'input': seq})

    temp_csv_in = "nanocdrx_input.csv"
    temp_csv_out = "nanocdrx_output.csv"

    with open(temp_csv_in, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['identifier', 'input'])
        writer.writeheader()
        for item in sequences:
            writer.writerow({'identifier': item['identifier'], 'input': item['input']})

    # 2. Run nanocdrx
    # Pass any extra arguments to predict_cdrs
    cmd = ["predict_cdrs", "-i", temp_csv_in, "-o", temp_csv_out] + unknown
    subprocess.run(cmd, check=True)

    # 3. Process Output
    fasta_out_name = f"{prefix}_cdr3.fasta"
    tsv_out_name = f"{prefix}_cdr3.tsv"
    hist_out_name = f"{prefix}_cdr3.hist"

    unique_cdr3s = set()
    cdr3_lengths = defaultdict(int)

    # Initialize histogram for 0-50 like getcdr3
    for i in range(51):
        cdr3_lengths[i] = 0

    with open(temp_csv_out, 'r') as f_in, \
         open(fasta_out_name, 'w') as f_fasta, \
         open(tsv_out_name, 'w') as f_tsv:

        reader = csv.DictReader(f_in)
        # TSV Header from getcdr3: ID\tCDR3\tsequence\tunique
        f_tsv.write("ID\tCDR3\tsequence\tunique\n")

        for row in reader:
            identifier = row['identifier']
            full_sequence = row['input']
            cdr3 = row.get('predicted_cdr3', '')

            # Logic from getcdr3.py
            status = ""

            if not cdr3 or cdr3 == "nan":
                status = "no-cdr3"
                f_tsv.write(f"{identifier}\tNA\t{full_sequence}\t{status}\n")
                continue

            if len(cdr3) < 1 or len(cdr3) > 50:
                status = "non-unique"
                f_tsv.write(f"{identifier}\t{cdr3}\t{full_sequence}\t{status}\n")
                continue

            if cdr3 in unique_cdr3s:
                status = "non-unique"
                f_tsv.write(f"{identifier}\t{cdr3}\t{full_sequence}\t{status}\n")
                continue

            # It is unique and valid
            unique_cdr3s.add(cdr3)
            cdr3_lengths[len(cdr3)] += 1

            # Write to FASTA
            # getcdr3 header: fastaheader (which is >ID)
            f_fasta.write(f">{identifier}\n{cdr3}\n")

            # Write to TSV
            f_tsv.write(f"{identifier}\t{cdr3}\t{full_sequence}\tunique\n")

    # 4. Write Histogram
    with open(hist_out_name, 'w') as f_hist:
        f_hist.write("Size,Count\n")
        for i in range(51):
            f_hist.write(f"{i},{cdr3_lengths[i]}\n")

if __name__ == "__main__":
    main()

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
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">") or line.startswith("@"):
                # Save the pending sequence if it exists
                if current_id:
                    full_seq = "".join(current_seq).upper()
                    # Strict cleaning: keep only standard amino acids
                    full_seq = "".join(c for c in full_seq if c in "ACDEFGHIKLMNPQRSTVWY")
                    # Biological VHH length filter (70-160 AA)
                    if 70 <= len(full_seq) <= 160:
                        sequences.append({'identifier': current_id, 'input': full_seq})
                # Start a new sequence
                # Replace comma in identifier to avoid CSV issues, use only the first part as ID
                current_id = str(line.lstrip(">@").split()[0]).replace(",", "_")
                current_seq = []
            else:
                current_seq.append(str(line))
        # Don't forget the last sequence in the file
        if current_id:
            full_seq = "".join(current_seq).upper()
            full_seq = "".join(c for c in full_seq if c in "ACDEFGHIKLMNPQRSTVWY")
            if 70 <= len(full_seq) <= 160:
                sequences.append({'identifier': current_id, 'input': full_seq})

    # 2. Run nanocdrx in Batches
    # Processing in chunks (e.g., 500) reduces memory usage and prevents inhomogeneous array shape errors on large datasets.
    batch_size = 500
    all_results = []
    
    import os
    seq_list = list(sequences)
    for i in range(0, len(seq_list), batch_size):
        batch = seq_list[i:i + batch_size]
        batch_num = (i // batch_size) + 1
        # print(f"Processing batch {batch_num} ({len(batch)} sequences)...")
        
        temp_csv_in = f"nanocdrx_input_b{batch_num}.csv"
        temp_csv_out = f"nanocdrx_output_b{batch_num}.csv"
        
        with open(temp_csv_in, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['identifier', 'input'])
            writer.writeheader()
            for item in batch:
                writer.writerow(item)
        
        # Run predict_cdrs on this batch
        cmd = ["predict_cdrs", "-i", temp_csv_in, "-o", temp_csv_out] + unknown
        try:
            subprocess.run(cmd, check=True)
            # Read results back
            with open(temp_csv_out, 'r') as f_out:
                reader = csv.DictReader(f_out)
                for row in reader:
                    all_results.append(row)
        except subprocess.CalledProcessError as e:
            # print(f"Error in batch {batch_num}: {e}")
            # If a batch fails, we skip it but continue with others to maximize output
            continue
        finally:
            # Cleanup temp files
            if os.path.exists(temp_csv_in): os.remove(temp_csv_in)
            if os.path.exists(temp_csv_out): os.remove(temp_csv_out)

    # 3. Process Final Output
    fasta_out_name = f"{prefix}_cdr3.fasta"
    tsv_out_name = f"{prefix}_cdr3.tsv"
    hist_out_name = f"{prefix}_cdr3.hist"

    unique_cdr3s = set()
    cdr3_lengths = {}

    # Initialize histogram for 0-50 like getcdr3
    for i in range(51):
        cdr3_lengths[i] = 0

    with open(fasta_out_name, 'w') as f_fasta, \
         open(tsv_out_name, 'w') as f_tsv:

        # TSV Header
        f_tsv.write("ID\tCDR3\tsequence\tunique\n")

        for row in all_results:
            identifier = row['identifier']
            full_sequence = row['input']
            cdr3 = row.get('predicted_cdr3', '')

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

            # Valid unique CDR3
            unique_cdr3s.add(cdr3)
            length = len(cdr3)
            cdr3_lengths[length] = cdr3_lengths.get(length, 0) + 1

            f_fasta.write(f">{identifier}\n{cdr3}\n")
            f_tsv.write(f"{identifier}\t{cdr3}\t{full_sequence}\tunique\n")

    # 4. Write Histogram
    with open(hist_out_name, 'w') as f_hist:
        f_hist.write("Size,Count\n")
        for i in range(51):
            count = cdr3_lengths.get(i, 0)
            f_hist.write(f"{i},{count}\n")

if __name__ == "__main__":
    main()

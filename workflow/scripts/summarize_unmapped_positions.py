#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--unmapped", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    record = next(SeqIO.parse(args.ref, "fasta"))
    chrom = record.id
    seq = str(record.seq)

    with open(args.out, "w") as out:
        out.write("sample\tchrom\tpos\tref_base\n")
        with open(args.unmapped) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                pos = int(line)
                ref_base = seq[pos-1] if 1 <= pos <= len(seq) else "N"
                out.write(f"{args.sample}\t{chrom}\t{pos}\t{ref_base}\n")

if __name__ == "__main__":
    main()
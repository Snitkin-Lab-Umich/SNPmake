#!/usr/bin/env python3
import argparse
import csv
import subprocess

def load_depth_table(depth_path):
    depth = {}
    with open(depth_path, "r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return depth

        # expected: Locus,Total_Depth,Average_Depth_sample,Depth_for_SAMPLE...
        for row in reader:
            if not row:
                continue
            locus = row[0].strip()
            if ":" not in locus:
                continue
            chrom, pos = locus.split(":")
            try:
                pos = int(pos)
                d = int(float(row[-1]))
            except ValueError:
                continue
            depth[(chrom, pos)] = d
    return depth

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--vcf", required=True, help="indel_prox_mask.vcf.gz")
    p.add_argument("--depth", required=True, help="*_depth_of_coverage (CSV)")
    p.add_argument("--out", required=True, help="output TSV")
    args = p.parse_args()

    depth = load_depth_table(args.depth)

    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n", args.vcf]
    res = subprocess.run(cmd, check=True, text=True, capture_output=True)

    with open(args.out, "w", newline="") as out:
        out.write("sample\tchrom\tpos\tref\talt\tfilter\tdoc_depth\tdoc_depth_missing\n")
        for line in res.stdout.splitlines():
            if not line:
                continue
            chrom, pos_s, ref, alt, flt = line.split("\t")
            pos = int(pos_s)
            d = depth.get((chrom, pos))
            out.write(
                f"{args.sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t{flt}\t"
                f"{'' if d is None else d}\t{1 if d is None else 0}\n"
            )

if __name__ == "__main__":
    main()
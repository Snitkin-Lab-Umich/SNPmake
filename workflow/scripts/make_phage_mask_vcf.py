#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys

def run(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True)

def fetch_ref_seq(ref_fasta, chrom, start, end):
    res = run(["samtools", "faidx", ref_fasta, f"{chrom}:{start}-{end}"])
    lines = res.stdout.splitlines()
    seq = "".join([l.strip() for l in lines if l and not l.startswith(">")]).upper()
    return seq

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True, help="Reference FASTA")
    p.add_argument("--regions_json", required=True, help="PHASTEST predicted_phage_regions.json")
    p.add_argument("--chrom", required=True, help="Chrom/contig name in reference/VCF (e.g., contig_1)")
    p.add_argument("--vcf_sample", required=True, help="Sample name to use in VCF header (must match PASS VCF sample)")
    p.add_argument("--out", required=True, help="Output VCF (uncompressed .vcf)")
    p.add_argument(
        "--mask_completeness",
        default="intact,questionable,incomplete",
        help="Comma-separated completeness classes to mask (default: all)"
    )
    args = p.parse_args()

    allowed = {x.strip().lower() for x in args.mask_completeness.split(",") if x.strip()}

    with open(args.regions_json) as f:
        regions = json.load(f)

    kept = []
    for r in regions:
        comp = str(r.get("completeness", "")).strip().lower()
        if not comp or comp in allowed:
            kept.append(r)

    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FILTER=<ID=MASKED_PHAGE,Description="PHASTEST-predicted prophage region masked to N in consensus">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.vcf_sample}\n")
        
        for r in kept:
            start = int(r["start"])
            stop = int(r["stop"])
            if stop < start:
                start, stop = stop, start

            ref_seq = fetch_ref_seq(args.ref, args.chrom, start, stop)
            if len(ref_seq) != (stop - start + 1):
                raise RuntimeError(f"Reference fetch length mismatch for {args.chrom}:{start}-{stop}")

            for i, pos in enumerate(range(start, stop + 1)):
                ref_base = ref_seq[i] if ref_seq else "N"
                out.write(f"{args.chrom}\t{pos}\t.\t{ref_base}\tN\t.\tMASKED_PHAGE\t.\tGT\t1/1\n")

if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr)
        raise
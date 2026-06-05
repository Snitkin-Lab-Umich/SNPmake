#!/usr/bin/env python3
import argparse
import csv
from cyvcf2 import VCF

def load_depth_table(depth_path):
    """
    Reads DepthOfCoverage file:
      Locus,Total_Depth,Average_Depth_sample,Depth_for_<sample>
      contig_1:42553,63,63.00,63
    Returns dict[(chrom, pos)] = depth (int)
    """
    depth = {}
    with open(depth_path, "r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return depth

        # Handle leading '*' in file
        header0 = header[0].lstrip("* ").strip()
        if header0.lower() != "locus":
            # If file has comment lines rewind style parsing
            # fall back to manual line parsing
            fh.seek(0)
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # skip header lines
                if "Locus" in line:
                    continue
                parts = [p.strip() for p in line.lstrip("* ").split(",")]
                if len(parts) < 2:
                    continue
                locus = parts[0]
                d = parts[-1]
                if ":" not in locus:
                    continue
                chrom, pos = locus.split(":")
                try:
                    depth[(chrom, int(pos))] = int(float(d))
                except ValueError:
                    continue
            return depth

        # Normal CSV path
        for row in reader:
            if not row:
                continue
            row0 = row[0].lstrip("* ").strip()
            if ":" not in row0:
                continue
            chrom, pos = row0.split(":")
            try:
                pos = int(pos)
            except ValueError:
                continue
            # the last column is sample depth
            try:
                d = int(float(row[-1]))
            except ValueError:
                d = None
            if d is not None:
                depth[(chrom, pos)] = d
    return depth

def info_get(rec, key):
    """cyvcf2: rec.INFO.get(key) returns None if absent"""
    try:
        return rec.INFO.get(key)
    except Exception:
        return None

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True, help="Final PASS SNP VCF (.vcf.gz)")
    p.add_argument("--depth", required=True, help="DepthOfCoverage file")
    p.add_argument("--sample", required=True, help="Sample name")
    p.add_argument("--out", required=True, help="Output TSV")
    args = p.parse_args()

    depth_map = load_depth_table(args.depth)
    vcf = VCF(args.vcf)

    # If VCF has one sample, we can attempt GT
    has_sample = (vcf.samples is not None and len(vcf.samples) >= 1)

    cols = [
        "sample","chrom","pos","ref","alt",
        "qual","filter",
        "vcf_DP","vcf_MQ","vcf_AF1","vcf_FQ","vcf_DP4",
        "gt",
        "doc_depth"
    ]

    with open(args.out, "w", newline="") as out:
        out.write("\t".join(cols) + "\n")

        for rec in vcf:
            # Feeding final VCF which already contains only passing records,
            # but keep FILTER anyway for transparency
            chrom = rec.CHROM
            pos = rec.POS
            ref = rec.REF
            alt = rec.ALT[0] if rec.ALT else "."

            # Skip non-SNP/weird records if any
            if len(ref) != 1 or len(alt) != 1:
                continue
            if alt in (".",):
                continue

            vcf_dp = info_get(rec, "DP")
            vcf_mq = info_get(rec, "MQ")
            vcf_af1 = info_get(rec, "AF1")
            vcf_fq = info_get(rec, "FQ")
            vcf_dp4 = info_get(rec, "DP4")

            gt = "."
            if has_sample:
                try:
                    # cyvcf2 genotype: [a1, a2, phased]
                    g = rec.genotypes[0]
                    if g is not None and len(g) >= 2:
                        gt = f"{g[0]}/{g[1]}"
                except Exception:
                    gt = "."

            doc_depth = depth_map.get((chrom, pos), "")

            row = [
                args.sample, chrom, str(pos), ref, alt,
                str(rec.QUAL) if rec.QUAL is not None else "",
                str(rec.FILTER) if rec.FILTER is not None else "",
                "" if vcf_dp is None else str(vcf_dp),
                "" if vcf_mq is None else str(vcf_mq),
                "" if vcf_af1 is None else str(vcf_af1),
                "" if vcf_fq is None else str(vcf_fq),
                "" if vcf_dp4 is None else str(vcf_dp4),
                gt,
                str(doc_depth)
            ]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
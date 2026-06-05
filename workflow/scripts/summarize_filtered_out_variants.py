import argparse
import csv
from cyvcf2 import VCF

def load_depth_table(depth_path):
    depth = {}
    with open(depth_path, "r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return depth

        header0 = header[0].lstrip("* ").strip()
        if header0.lower() != "locus":
            fh.seek(0)
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
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
            try:
                d = int(float(row[-1]))
            except ValueError:
                d = None
            if d is not None:
                depth[(chrom, pos)] = d
    return depth

def info_get(rec, key):
    try:
        return rec.INFO.get(key)
    except Exception:
        return None

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw", required=True, help="Raw SNP VCF after 5bp-indel removal (.vcf.gz)")
    p.add_argument("--final", required=True, help="Final PASS SNP VCF (.vcf.gz)")
    p.add_argument("--depth", required=True, help="DepthOfCoverage file")
    p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    depth_map = load_depth_table(args.depth)

    # Load FINAL positions into a set (chrom,pos)
    final_pos = set()
    vcf_final = VCF(args.final)
    for rec in vcf_final:
        final_pos.add((rec.CHROM, rec.POS))

    vcf_raw = VCF(args.raw)

    cols = [
        "sample","chrom","pos","ref","alt",
        "qual",
        "DP","MQ","AF1","FQ","DP4",
        "doc_depth","doc_depth_missing"
    ]

    with open(args.out, "w", newline="") as out:
        out.write("\t".join(cols) + "\n")

        for rec in vcf_raw:
            key = (rec.CHROM, rec.POS)
            if key in final_pos:
                continue  # not filtered out

            ref = rec.REF
            alt = rec.ALT[0] if rec.ALT else "."

            # Keep SNPs only 
            if len(ref) != 1 or len(alt) != 1:
                continue
            if alt in (".",):
                continue

            DP = info_get(rec, "DP")
            MQ = info_get(rec, "MQ")
            AF1 = info_get(rec, "AF1")
            FQ = info_get(rec, "FQ")
            DP4 = info_get(rec, "DP4")

            doc_depth = depth_map.get(key, "")
            doc_missing = "1" if doc_depth == "" else "0"

            row = [
                args.sample, rec.CHROM, str(rec.POS), ref, alt,
                str(rec.QUAL) if rec.QUAL is not None else "",
                "" if DP is None else str(DP),
                "" if MQ is None else str(MQ),
                "" if AF1 is None else str(AF1),
                "" if FQ is None else str(FQ),
                "" if DP4 is None else str(DP4),
                str(doc_depth),
                doc_missing
            ]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
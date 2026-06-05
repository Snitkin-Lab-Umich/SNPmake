import argparse
import subprocess
from Bio import SeqIO

def load_reference(ref_fasta):
    ref = {}
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        ref[rec.id] = str(rec.seq).upper()
    return ref

def load_positions_txt(path):
    s = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            s.add(int(line))
    return s

def bcftools_query_lines(cmd):
    res = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return res.stdout.splitlines()

def snp_positions_from_vcf(vcf_gz, chrom=None):
    """
    Return set of (chrom,pos) for SNP-only records:
      len(REF)==1 and all ALT alleles len==1 (multi-allelic SNPs allowed).
    """
    cmd = ["bcftools", "query"]
    if chrom:
        cmd += ["-r", chrom]
    cmd += ["-f", "%CHROM\t%POS\t%REF\t%ALT\n", vcf_gz]

    s = set()
    for line in bcftools_query_lines(cmd):
        if not line:
            continue
        c, pos_s, ref, alt = line.split("\t")
        pos = int(pos_s)

        if alt in (".", "") or alt.startswith("<"):
            continue

        alts = alt.split(",")
        if len(ref) != 1:
            continue
        if any(len(a) != 1 for a in alts):
            continue

        s.add((c, pos))
    return s

def positions_from_vcf(vcf_gz, chrom=None):
    """Return set of (chrom,pos) for all records in VCF (no TYPE filtering)."""
    cmd = ["bcftools", "query"]
    if chrom:
        cmd += ["-r", chrom]
    cmd += ["-f", "%CHROM\t%POS\n", vcf_gz]

    s = set()
    for line in bcftools_query_lines(cmd):
        if not line:
            continue
        c, pos_s = line.split("\t")
        s.add((c, int(pos_s)))
    return s

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--raw", required=True, help="*_aln_mpileup_raw.vcf.gz")
    p.add_argument("--indel_removed", required=True, help="*_aln_mpileup_raw.vcf_5bp_indel_removed.vcf.gz")
    p.add_argument("--chrom", required=False, default=None,
                   help="Contig name (recommended; required if --unmapped_positions is used)")
    p.add_argument("--unmapped_positions", required=False, default=None,
                   help="*_unmapped.bed_positions (1-based positions, one per line)")
    p.add_argument("--phage_mask_vcf", required=False, default=None,
                   help="optional phage mask VCF (.vcf.gz) to exclude from INDEL_PROX")
    p.add_argument("--vcf_sample", required=True)
    p.add_argument("--out", required=True, help="output .vcf (uncompressed)")
    args = p.parse_args()

    if args.unmapped_positions and not args.chrom:
        raise SystemExit("--chrom is required when using --unmapped_positions (unmapped file has positions only).")

    ref = load_reference(args.ref)

    raw_snps = snp_positions_from_vcf(args.raw, chrom=args.chrom)
    kept_snps = snp_positions_from_vcf(args.indel_removed, chrom=args.chrom)

    # INDEL_PROX = SNPs removed by 5bp indel filter
    indel_prox = set(raw_snps - kept_snps)

    # Exclude UNMAPPED positions (dash wins; keep categories exclusive)
    if args.unmapped_positions:
        unmapped_pos = load_positions_txt(args.unmapped_positions)
        indel_prox -= {(args.chrom, p) for p in unmapped_pos}

    # Exclude PHAGE positions (if enabled; keep categories exclusive)
    if args.phage_mask_vcf:
        phage_pos = positions_from_vcf(args.phage_mask_vcf, chrom=args.chrom)
        indel_prox -= phage_pos

    indel_prox = sorted(indel_prox, key=lambda x: (x[0], x[1]))

    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FILTER=<ID=MASKED_INDEL_PROX,Description="SNP removed by 5bp-indel filter (excluding unmapped/phage); masked to N in consensus">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.vcf_sample}\n")

        for chrom, pos in indel_prox:
            seq = ref.get(chrom)
            ref_base = seq[pos-1] if (seq and 1 <= pos <= len(seq)) else "N"
            out.write(f"{chrom}\t{pos}\t.\t{ref_base}\tN\t.\tMASKED_INDEL_PROX\t.\tGT\t1/1\n")

if __name__ == "__main__":
    main()
import argparse
import subprocess
from Bio import SeqIO

def load_reference(ref_fasta):
    ref = {}
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        ref[rec.id] = str(rec.seq)
    return ref

def bcftools_positions(vcf_gz):
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\n", vcf_gz]
    res = subprocess.run(cmd, check=True, text=True, capture_output=True)
    pos_set = set()
    for line in res.stdout.splitlines():
        if not line:
            continue
        chrom, pos = line.split("\t")
        pos_set.add((chrom, int(pos)))
    return pos_set

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--raw", required=True)
    p.add_argument("--final", required=True)
    p.add_argument("--vcf_sample", required=True)   
    p.add_argument("--out", required=True)
    args = p.parse_args()

    ref = load_reference(args.ref)
    raw_pos = bcftools_positions(args.raw)
    final_pos = bcftools_positions(args.final)
    filtered = sorted(raw_pos - final_pos)

    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FILTER=<ID=MASKED_FILTERED_OUT,Description="Site present in raw SNP VCF but absent from final SNP VCF; masked to N in consensus">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.vcf_sample}\n")
        for chrom, pos in filtered:
            seq = ref.get(chrom)
            ref_base = seq[pos-1] if (seq and 1 <= pos <= len(seq)) else "N"
            out.write(f"{chrom}\t{pos}\t.\t{ref_base}\tN\t.\tMASKED_FILTERED_OUT\t.\tGT\t1/1\n")

if __name__ == "__main__":
    main()
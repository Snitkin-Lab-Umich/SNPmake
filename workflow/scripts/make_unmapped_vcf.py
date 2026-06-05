import argparse
import subprocess
import sys
import tempfile

def run(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True)

def fetch_ref_bases(ref_fasta, chrom, positions):
    """
    Fetch reference base for each position using samtools faidx in batch via -r regions_file.
    Returns dict[pos] = base.
    """
    regions = [f"{chrom}:{p}-{p}" for p in positions]
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tf:
        for r in regions:
            tf.write(r + "\n")
        region_file = tf.name

    try:
        res = run(["samtools", "faidx", "-r", region_file, ref_fasta])
    finally:
        try:
            import os
            os.unlink(region_file)
        except Exception:
            pass

    bases = {}
    current_pos = None
    for line in res.stdout.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            hdr = line[1:].strip()              
            chrom_part, rng = hdr.split(":")
            start, end = rng.split("-")
            current_pos = int(start)
        else:
            base = line.strip().upper()
            if current_pos is not None:
                bases[current_pos] = base[0] if base else "N"
                current_pos = None
    return bases

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True)
    p.add_argument("--unmapped", required=True)
    p.add_argument("--chrom", required=True, help="Reference contig name (e.g., contig_1)")
    p.add_argument("--vcf_sample", required=True, help="Sample name to use in VCF header (must match PASS VCF sample)")
    p.add_argument("--out", required=True, help="Output VCF (uncompressed .vcf)")
    args = p.parse_args()

    # read positions
    positions = []
    with open(args.unmapped) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            positions.append(int(line))

    positions = sorted(set(positions))
    bases = fetch_ref_bases(args.ref, args.chrom, positions) if positions else {}

    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FILTER=<ID=UNMAPPED,Description="Position flagged as unmapped/low coverage; masked in consensus">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.vcf_sample}\n")

        for pos in positions:
            ref_base = bases.get(pos, "N")
            # Mask as N and set GT to 1/1 (homozygous alt) to indicate masked position in consensus
            out.write(f"{args.chrom}\t{pos}\t.\t{ref_base}\tN\t.\tUNMAPPED\t.\tGT\t1/1\n")

if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr)
        raise
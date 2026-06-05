__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
import gzip
from snakemake.shell import shell


def open_text(path):
    """
    Open plain or gzipped text.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


snp_vcf_file = snakemake.input.snp_vcf
raw_vcf_file = snakemake.input.raw_vcf

output_vcf = snakemake.output.vcf
excluded_positions_file = snakemake.output.excluded_positions
output_vcf_gz = snakemake.output.vcf_gz


def collect_og_style_indel_positions(raw_vcf_file):
    """
    Mimic OG SNPkit indel detection.

    So we collect POS from raw VCFs whose INFO field starts with 'INDEL;'.
    """
    indel_positions = set()

    with open_text(raw_vcf_file) as raw:
        for line in raw:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            pos = int(fields[1])
            info = fields[7]

            if info.startswith("INDEL;"):
                indel_positions.add(pos)

    return indel_positions


def remove_5_bp_snp_indel_og_style(snp_vcf_file, raw_vcf_file, output_vcf, excluded_positions_file):
    """
    Remove SNP records whose POS falls within +/-5 bp of an OG style raw INDEL position.

    Important:
    - This uses raw VCF INFO.startswith('INDEL;') to find indels matching OG SNPkit behavior.
    - It writes unique excluded positions that are sorted.
    """
    indel_positions = collect_og_style_indel_positions(raw_vcf_file)

    exclude_positions = set()
    for indel_pos in indel_positions:
        for p in range(indel_pos - 5, indel_pos + 6):
            if p > 0:
                exclude_positions.add(p)

    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    os.makedirs(os.path.dirname(excluded_positions_file), exist_ok=True)

    with open(excluded_positions_file, "w") as out:
        for pos in sorted(exclude_positions):
            out.write(f"{pos}\n")

    with open(output_vcf, "w") as out, open_text(snp_vcf_file) as snps:
        for line in snps:
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue

            pos = int(fields[1])

            if pos not in exclude_positions:
                out.write(line)

    return output_vcf


remove_5_bp_snp_indel_og_style(
   snp_vcf_file,
   raw_vcf_file,
   output_vcf,
   excluded_positions_file
)


shell(
    """
    set -euo pipefail

    bgzip -f -c {output_vcf} > {output_vcf_gz}
    tabix -f -p vcf {output_vcf_gz}
    """
)
__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell

outdir_temp = snakemake.params.get("outdir_temp", "")
prefix = snakemake.params.get("prefix", "") 
ref_genome = snakemake.params.get("ref_genome", "")

#samclip command
shell("samclip --ref {ref_genome} --max 10 < {snakemake.input.aligned_sam_out} > {snakemake.output.clipped_sam_out}")

# samtools command
shell("samtools view -Sb {snakemake.output.clipped_sam_out} > {snakemake.output.bam_out} && samtools sort {snakemake.output.bam_out} -m 500M -@ 0 -o {snakemake.output.sorted_bam_out} -T {outdir_temp} &> {snakemake.log.post_align_sam_to_bam_log}")
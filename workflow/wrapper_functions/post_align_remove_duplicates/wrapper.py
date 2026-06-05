__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell

outdir = snakemake.params.get("outdir", "")
outdir_dups_removed = snakemake.params.get("outdir_dups_removed", "")
prefix = snakemake.params.get("prefix", "") 

#picard command 
shell("picard MarkDuplicates -REMOVE_DUPLICATES true -INPUT {snakemake.input.sorted_bam_out} -OUTPUT {snakemake.output.bam_duplicates_removed_out} -METRICS_FILE {outdir_dups_removed}/{prefix}_markduplicates_metrics -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT &> {snakemake.log.picard_dups_log}")

# samtools command
shell("samtools sort {snakemake.output.bam_duplicates_removed_out} -m 500M -@ 0 -o {snakemake.output.dups_rmvd_sorted_bam_out} -T {outdir}/{prefix}_aln_sort_temp"
      " && samtools index {snakemake.output.dups_rmvd_sorted_bam_out}  &> {snakemake.log.post_align_remove_duplicates_log}")

__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell

# params
ref_genome = snakemake.params.get("ref_genome", "")

# snps
dp_snp_filter = snakemake.params.get("dp_snp_filter", "")
fq_snp_filter = snakemake.params.get("fq_snp_filter", "")
mq_snp_filter = snakemake.params.get("mq_snp_filter", "")
qual_snp_filter = snakemake.params.get("qual_snp_filter", "")
af_snp_filter = snakemake.params.get("af_snp_filter", "")

# indels
dp_indel_filter = snakemake.params.get("dp_indel_filter", "")
mq_indel_filter = snakemake.params.get("mq_indel_filter", "")
qual_indel_filter = snakemake.params.get("qual_indel_filter", "")
af_indel_filter = snakemake.params.get("af_indel_filter", "")

# snps
gatk_snp_filter_parameter_expression= "%s && %s && %s && %s && %s" % (dp_snp_filter, fq_snp_filter, mq_snp_filter, qual_snp_filter, af_snp_filter)
shell("gatk VariantFiltration -R {ref_genome} -O {snakemake.output.filter_snp_vcf} --variant {snakemake.input.final_raw_snp_vcf} --filter-expression \"{gatk_snp_filter_parameter_expression}\" --filter-name PASS_filter &> {snakemake.log.gatk_snp}") 
shell("grep '#\|PASS_filter' {snakemake.output.filter_snp_vcf} > {snakemake.output.filter_snp_final}")
shell("""
       bgzip -c {snakemake.output.filter_snp_final} > {snakemake.output.zipped_filtered_snp_vcf} &&
       tabix -p vcf -f {snakemake.output.zipped_filtered_snp_vcf}
    """)

shell("""
       bgzip -c {snakemake.output.filter_snp_vcf} > {snakemake.output.zipped_filter_snp_vcf} &&
       tabix -p vcf -f {snakemake.output.zipped_filter_snp_vcf}
    """)

# indels       
gatk_indel_filter_parameter_expression="%s && %s && %s && %s" % (dp_indel_filter, mq_indel_filter, qual_indel_filter, af_indel_filter)
shell("gatk VariantFiltration -R {ref_genome} -O {snakemake.output.filter_indel_vcf} --variant {snakemake.input.final_raw_indel_vcf} --filter-expression \"{gatk_indel_filter_parameter_expression}\" --filter-name PASS_filter&> {snakemake.log.gatk_indels}")
shell("grep '#\|PASS_filter' {snakemake.output.filter_indel_vcf} > {snakemake.output.filter_indel_final}")

shell("""
       bgzip -c {snakemake.output.filter_indel_final} > {snakemake.output.zipped_filtered_indel_vcf} &&
       tabix -p vcf -f {snakemake.output.zipped_filtered_indel_vcf}
    """)

shell("""
       bgzip -c {snakemake.output.filter_indel_vcf} > {snakemake.output.zipped_filter_indel_vcf} &&
       tabix -p vcf -f {snakemake.output.zipped_filter_indel_vcf}
    """)
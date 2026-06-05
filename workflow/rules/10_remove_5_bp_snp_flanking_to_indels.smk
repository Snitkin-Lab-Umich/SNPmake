# remove snps with 5bp of an indel 
rule remove_5_bp_snp_flanking_to_indels:
    input:
        snp_vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_filter_snp_final.vcf",
        raw_vcf="results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf"
    output:
        vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf",
        excluded_positions="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_excluded_positions.txt",
        vcf_gz="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",
        tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz.tbi" 
    benchmark: 
        "benchmarks/{prefix}/remove_5bp_snp_flanking_to_indels/{sample}.benchmark.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    wrapper:
        "file:workflow/wrapper_functions/remove_5_bp_snp_flanking_to_indels"
    
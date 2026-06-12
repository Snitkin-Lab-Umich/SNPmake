# variant calling 
# gatk
# calling snp/indel and subset of variants using gatk
rule gatk_call_indels:
    input:
        index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam",
    output:
        final_raw_vcf= "results/{prefix}/gatk_varcall/{sample}/{sample}_aln_mpileup_raw.vcf",
        indel_file = temp("results/{prefix}/gatk_varcall/{sample}/{sample}_indel.vcf"),
        zipped_indel_vcf = "results/{prefix}/gatk_varcall/{sample}/{sample}_indel.vcf.gz"
    params:
        ref_genome = config["reference_genome"]
    log:
        "logs/{prefix}/{sample}/gatk/{sample}_gatk.log"
    benchmark: 
        "benchmarks/{prefix}/gatk_call_indels/{sample}.benchmark.tsv"
    threads: 3
    resources:
        mem_mb=4000,
        runtime=60
    wrapper:
        "file:workflow/wrapper_functions/gatk_call_indel"

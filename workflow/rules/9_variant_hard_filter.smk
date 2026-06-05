rule variant_hard_filter:
    input:
        final_raw_snp_vcf = "results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf",
        final_raw_indel_vcf = "results/{prefix}/gatk_varcall/{sample}/{sample}_indel.vcf"
    output:
        filter_snp_vcf = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_snp.vcf",
        filter_snp_final = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_snp_final.vcf",
        filter_indel_vcf = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_indel.vcf",
        filter_indel_final = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_indel_final.vcf",
        zipped_filtered_snp_vcf = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_snp_final.vcf.gz",
        zipped_filtered_indel_vcf = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_indel_final.vcf.gz",
        zipped_filter_snp_vcf = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_snp.vcf.gz",
        zipped_filter_indel_vcf = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_indel.vcf.gz",
    params:
        ref_genome = config["reference_genome"],
        dp_indel_filter = config["dp_indel_filter"],
        mq_indel_filter = config["mq_indel_filter"],
        qual_indel_filter = config["qual_indel_filter"],
        af_indel_filter = config["af_indel_filter"],
        dp_snp_filter = config["dp_snp_filter"],
        fq_snp_filter = config["fq_snp_filter"],
        mq_snp_filter = config["mq_snp_filter"],
        qual_snp_filter = config["qual_snp_filter"],
        af_snp_filter = config["af_snp_filter"]
    log:
        gatk_snp = "logs/{prefix}/gatk/{sample}/{sample}_filtered_snps.log",
        gatk_indels = "logs/{prefix}/gatk/{sample}/{sample}_filtered_indels.log"
    benchmark: 
        "benchmarks/{prefix}/variant_hard_filter/{sample}.benchmark.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=10
    wrapper:
        "file:workflow/wrapper_functions/hard_filter"
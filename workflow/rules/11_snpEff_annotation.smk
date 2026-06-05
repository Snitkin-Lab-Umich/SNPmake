rule snpEff_annotation:
    input:
        remove_snps_5_bp_snp_indel_file = "results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf",
        filter_snp_final = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_snp_final.vcf",
        filter_indel_final = "results/{prefix}/filtered_vcf/{sample}/{sample}_filter_indel_final.vcf"
    output:
        csv_summary_file = "results/{prefix}/annotated_files/{sample}/{sample}_ANN.csv",
        annotated_vcf = "results/{prefix}/annotated_files/{sample}/{sample}_ANN.vcf",
        annotated_filter_snp_final = "results/{prefix}/annotated_files/{sample}/{sample}_filter_snp_final_ANN.vcf",
        annotated_filter_indel_final = "results/{prefix}/annotated_files/{sample}/{sample}_filter_indel_final_ANN.vcf",
        summary_file_filter_snp_final = "results/{prefix}/annotated_files/{sample}/{sample}_filter_snp_final_ANN.csv",
        summary_file_filter_indel_final = "results/{prefix}/annotated_files/{sample}/{sample}_filter_indel_final_ANN.csv"
    params:
        snpeff_parameters = config["snpeff_parameters"],
        snpEff_db = REF_NAME, 
        base_dir = my_basedir
    log:
        "logs/{prefix}/snpEff_annotation/{sample}/{sample}_snpEff.log"
    benchmark: 
        "benchmarks/{prefix}/snpEff_annotation/{sample}.benchmark.tsv"
    threads: 3
    resources:
        mem_mb=1000,
        runtime=15
    wrapper:
        "file:workflow/wrapper_functions/snpEff_annotation"

def parse_bed_file(final_bed_unmapped_file):
    unmapped_positions_array = []
    with open(final_bed_unmapped_file, 'r') as fp:
        for line in fp:
            line_array = line.split('\t')
            lower_index = int(line_array[1]) + 1
            upper_index = int(line_array[2]) + 1
            for positions in range(lower_index,upper_index):
                unmapped_positions_array.append(positions)
    only_unmapped_positions_file = final_bed_unmapped_file + "_positions"
    f1=open(only_unmapped_positions_file, 'w+')
    for i in unmapped_positions_array:
        p_string = str(i) + "\n"
        f1.write(p_string)
    return only_unmapped_positions_file

# determine statistics of file
rule alignment_stats:
    input:
        index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam"
    output:
        alignment_stats = "results/{prefix}/stats/{sample}/{sample}_alignment_stats.tsv" 
    singularity:
        "docker://staphb/samtools:1.23.1"
    log:
        "logs/{prefix}/post_align/{sample}/{sample}_stats.log"
    benchmark: 
        "benchmarks/{prefix}/alignment_stats/{sample}.benchmark.tsv"
    shell:
        "samtools flagstat {input.index_sorted_dups_rmvd_bam} > {output.alignment_stats} &> {log}" 

# determine coverage of bam file       
rule gatk_coverage_depth_statistics:
    input:
        index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam",
        intervals = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME)
    output:
        gatk_depthCoverage_summary = "results/{prefix}/stats/{sample}/{sample}_depth_of_coverage.sample_summary"
    params:
        outdir = "results/{prefix}/stats/{sample}",
        ref_genome = config["reference_genome"], 
        prefix = "{sample}",
    log:
        "logs/{prefix}/post_align/{sample}/{sample}_coverage_depth.log"
    singularity:
        "docker://broadinstitute/gatk:4.6.2.0" 
    benchmark: 
        "benchmarks/{prefix}/gatk_coverage_depth_statistics/{sample}.benchmark.tsv"
    shell:
        "gatk DepthOfCoverage -R {params.ref_genome} -O {params.outdir}/{params.prefix}_depth_of_coverage -I {input.index_sorted_dups_rmvd_bam} --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {input.intervals} &> {log}"

# bedtools
rule bedtools_extract_coverage:
    input:
        index_sorted_dups_rmvd_bam ="results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam"
    output:
        unmapped_bed = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed"
    log:
       "logs/{prefix}/bedtools_unmapped/{sample}/{sample}.log"
    singularity:
        "docker://staphb/bedtools:2.31.1"
    benchmark: 
        "benchmarks/{prefix}/bedtools_extract_coverage/{sample}.benchmark.tsv"
    shell:
        "bedtools genomecov -ibam {input.index_sorted_dups_rmvd_bam} -bga | awk '$4==0' > {output.unmapped_bed}"

# returns unmapped positions file    
rule parse_bed_file_find_unmapped_regions:
    input:
        unmapped_bed = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed"
    output:
        unmapped_bam_positions = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed_positions"
    benchmark: 
        "benchmarks/{prefix}/parse_bed_file_find_unmapped_regions/{sample}.benchmark.tsv"  
    run:
        parse_bed_file(input.unmapped_bed)

# prepared reference size and window files
rule prepare_reference_windows:
    output:
        reference_size_file="results/{prefix}/ref_genome_files/{ref_name}.size",
        reference_window_file = "results/{prefix}/ref_genome_files/{ref_name}.bed"
    params:
        ref_genome = config["reference_genome"]
    # benchmark: 
    #     "benchmarks/{prefix}/prepare_reference_windows/benchmark.tsv"
    wrapper:
        "file:workflow/wrapper_functions/prepare_reference_files"

rule bedcoverage:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/post_align/{wildcards.sample}/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
        reference_window_file = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME)
    output:
        bedgraph_cov = f"results/{{prefix}}/bedtools/{{sample}}/bedgraph_coverage/{{sample}}.bedcov"
    singularity:
        "docker://staphb/bedtools:2.31.1"
    log:
        "logs/{prefix}/{sample}/bedgraph_coverage/{sample}_bedcov.log"
    benchmark: 
        "benchmarks/{prefix}/bedcoverage/{sample}.benchmark.tsv"
    shell:
        "bedtools coverage -abam {input.index_sorted_dups_rmvd_bam} -b {input.reference_window_file} > {output.bedgraph_cov} 2>&1 | tee {log}"
        
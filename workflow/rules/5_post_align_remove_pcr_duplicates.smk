# remove duplicates and sort and index bam file with duplicates removed 
rule post_align_remove_pcr_duplicates:
    input:
        sorted_bam_out = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort.bam"
    output:
        bam_duplicates_removed_out = temp("results/{prefix}/post_align/{sample}/remove_duplicates/{sample}_aln_marked.bam"),
        dups_rmvd_sorted_bam_out = temp("results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam")
    params:
        outdir_dups_removed = "results/{prefix}/post_align/{sample}/remove_duplicates",
        outdir = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/",
        prefix = "{sample}"
    log:
        picard_dups_log = "logs/{prefix}/post_align/{sample}/{sample}_picard.log",
        post_align_remove_duplicates_log= "logs/{prefix}/post_align/{sample}/{sample}_samtools.log"
    benchmark: 
        "benchmarks/{prefix}/post_align_remove_pcr_duplicates/{sample}.benchmark.tsv"
    threads: 2
    resources:
        mem_mb=5000,
        runtime=15
    wrapper:
        "file:workflow/wrapper_functions/post_align_remove_duplicates"
    
# samclip and sort bam file
rule removed_clipped_reads:  
    input:
        aligned_sam_out = "results/{prefix}/align_reads/{sample}/{sample}_aln.sam"
    output:
        clipped_sam_out = "results/{prefix}/post_align/{sample}/samclip/{sample}_clipped.sam",
        bam_out = "results/{prefix}/post_align/{sample}/aligned_bam/{sample}_aln.bam",
        sorted_bam_out = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort.bam"
    params:
        outdir_temp = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort_temp",
        prefix = "{sample}",
        ref_genome= config["reference_genome"]
    log:
        samclip_log = "logs/{prefix}/post_align/{sample}/{sample}_samclip.log",
        post_align_sam_to_bam_log= "logs/{prefix}/post_align/{sample}/{sample}.log"
    benchmark: 
        "benchmarks/{prefix}/remove_clipped_reads/{sample}.benchmark.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15
    wrapper:
        "file:workflow/wrapper_functions/sam_to_bam"

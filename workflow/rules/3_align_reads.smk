rule align_reads:
    input:
        r1=align_r1,
        r2=align_r2
    output:
        aligned_sam_out="results/{prefix}/align_reads/{sample}/{sample}_aln.sam"
    params:
        ref_genome=config["reference_genome"]
    log:
        "logs/{prefix}/align_reads/{sample}/{sample}.log"
    singularity:
        "docker://staphb/bwa:0.7.19"
    benchmark: 
        "benchmarks/{prefix}/align_reads/{sample}.benchmark.tsv"
    threads: 8
    resources:
        mem_mb=3000,
        runtime=15
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})

        bash workflow/scripts/align_reads.sh \
            {input.r1} \
            {threads} \
            {params.ref_genome} \
            {input.r1} \
            {input.r2} \
            {output.aligned_sam_out} \
            > {log} 2>&1
        """

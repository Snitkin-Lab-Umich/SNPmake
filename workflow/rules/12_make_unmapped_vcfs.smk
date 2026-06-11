# Create unmapped vcfs
rule make_unmapped_vcfs:
    input:
        unmapped = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed_positions",
        ref=REF_GENOME,
        remove_snps_5_bp_snp_indel_file = "results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",   
    output:
        vcf="results/{prefix}/unmapped_vcf/{sample}/{sample}_unmapped.vcf.gz",
        tbi="results/{prefix}/unmapped_vcf/{sample}/{sample}_unmapped.vcf.gz.tbi"
    params:
        outdir="results/{prefix}/unmapped_vcf/{sample}",
        chrom=config["reference_contig"]
    # singularity:
    #     "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/make_unmapped_vcfs/{sample}.benchmark.tsv"
    envmodules:
        "Bioinformatics",
        "htslib",
        "bcftools",
        "samtools", 
        "python/3.13.2"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        r"""
        set -euo pipefail
        
        VCF_SAMPLE=$(bcftools query -l {input.remove_snps_5_bp_snp_indel_file} | head -n 1)

        python workflow/scripts/make_unmapped_vcf.py \
          --ref {input.ref} \
          --unmapped {input.unmapped} \
          --chrom {params.chrom} \
          --vcf_sample "$VCF_SAMPLE" \
          --out {params.outdir}/{wildcards.sample}_unmapped.vcf

        bgzip -f {params.outdir}/{wildcards.sample}_unmapped.vcf
        tabix -f -p vcf {output.vcf}
        """

rule summarize_unmapped_positions:
    input:
        unmapped = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed_positions",
        ref=REF_GENOME
    output:
        tsv="results/{prefix}/unmapped_positions/{sample}_unmapped_positions.tsv"
    params:
        sample="{sample}"
    benchmark: 
        "benchmarks/{prefix}/summarize_unmapped_positions/{sample}.benchmark.tsv"
    conda:
        "envs/make_filtered_mask_vcf.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        r"""
        set -euo pipefail # three error handlong behaviors: exit on error, treat unset variables as error, and fail if any command in a pipeline fails
        python workflow/scripts/summarize_unmapped_positions.py \
          --ref {input.ref} \
          --unmapped {input.unmapped} \
          --sample {params.sample} \
          --out {output.tsv}
        """

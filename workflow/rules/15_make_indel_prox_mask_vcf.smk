rule make_indel_prox_mask_vcf:
    input:
        ref=REF_GENOME,
        raw="results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf.gz",
        indel_removed="results/{prefix}/gatk_varcall/{sample}/{sample}_indel.vcf.gz",
        unmapped="results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed_positions",
        pass_vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",
        phage_mask=lambda wc: f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz" if PHAGE_ENABLED else [],
        phage_tbi=lambda wc: f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz.tbi" if PHAGE_ENABLED else [],
    output:
        vcf_gz="results/{prefix}/indel_prox_mask_vcf/{sample}/{sample}_indel_prox_mask.vcf.gz",
        tbi="results/{prefix}/indel_prox_mask_vcf/{sample}/{sample}_indel_prox_mask.vcf.gz.tbi"
    params:
        outdir="results/{prefix}/indel_prox_mask_vcf/{sample}",
        chrom=config["reference_contig"],
        phage_arg=phage_arg
    # singularity:
    #     "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/make_indel_prox_mask_vcf/{sample}.benchmark.tsv"
    # envmodules:
    #     "Bioinformatics",
    #     "bcftools",
    #     "htslib",
    #     "python/3.13.2"
    conda:
        "envs/make_filtered_mask_vcf.yaml"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        VCF_SAMPLE=$(bcftools query -l {input.pass_vcf} | head -n 1)

        python workflow/scripts/make_indel_prox_mask_vcf.py \
          --ref {input.ref} \
          --raw {input.raw} \
          --indel_removed {input.indel_removed} \
          --chrom {params.chrom} \
          --unmapped_positions {input.unmapped} \
          {params.phage_arg} \
          --vcf_sample "$VCF_SAMPLE" \
          --out {params.outdir}/{wildcards.sample}_indel_prox_mask.vcf

        bgzip -f {params.outdir}/{wildcards.sample}_indel_prox_mask.vcf
        tabix -f -p vcf {output.vcf_gz}
        """

rule summarize_indel_prox_positions:
    input:
        vcf="results/{prefix}/indel_prox_mask_vcf/{sample}/{sample}_indel_prox_mask.vcf.gz",
        depth="results/{prefix}/stats/{sample}/{sample}_depth_of_coverage.sample_summary"
    output:
        tsv="results/{prefix}/indel_prox_positions/{sample}_indel_prox_positions.tsv"
    # singularity:
    #     "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/summarize_indel_prox_positions/{sample}.benchmark.tsv"
    envmodules:
        "Bioinformatics",
        "bcftools" ,
        "python/3.13.2",
    # conda:
    #     "envs/cyvcf2.yaml"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        python workflow/scripts/summarize_indel_prox_positions.py \
          --sample {wildcards.sample} \
          --vcf {input.vcf} \
          --depth {input.depth} \
          --out {output.tsv}
        """

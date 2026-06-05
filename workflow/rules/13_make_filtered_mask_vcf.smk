rule make_filtered_mask_vcf:
    input:
        raw="results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf",
        final="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf",
        pass_vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf",
        ref=REF_GENOME
    output:
        vcf_gz="results/{prefix}/filtered_mask_vcf/{sample}/{sample}_filtered_mask.vcf.gz",
        tbi="results/{prefix}/filtered_mask_vcf/{sample}/{sample}_filtered_mask.vcf.gz.tbi"
    params:
        outdir="results/{prefix}/filtered_mask_vcf/{sample}"
    # singularity:
    #     "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/make_filtered_mask_vcf/{sample}.benchmark.tsv"
    conda:
        "envs/make_filtered_mask_vcf.yaml"
    # envmodules:
    #     "Bioinformatics",
    #     "htslib",
    #     "bcftools",
    #     "python/3.13.2"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        r"""
        set -euo pipefail

        VCF_SAMPLE=$(bcftools query -l {input.pass_vcf} | head -n 1)

        python workflow/scripts/make_filtered_mask_vcf.py \
          --ref {input.ref} \
          --raw {input.raw} \
          --final {input.final} \
          --vcf_sample "$VCF_SAMPLE" \
          --out {params.outdir}/{wildcards.sample}_filtered_mask.vcf
        
        bgzip -f {params.outdir}/{wildcards.sample}_filtered_mask.vcf
        tabix -f -p vcf {output.vcf_gz}
        """

rule summarize_filtered_out_variants:
    input:
        raw="results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf",
        final="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",
        depth="results/{prefix}/stats/{sample}/{sample}_depth_of_coverage.sample_summary"
    output:
        tsv="results/{prefix}/filtered_out_variants/{sample}_filtered_out_variants.tsv"
    benchmark: 
        "benchmarks/{prefix}/summarize_filtered_out_variants/{sample}.benchmark.tsv"
    conda:
        "envs/cyvcf2.yaml"
    # envmodules:
    #     "Bioinformatics",
    #     "python/3.13.2",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        """
        set -euo pipefail
        python workflow/scripts/summarize_filtered_out_variants.py \
          --raw {input.raw} \
          --final {input.final} \
          --depth {input.depth} \
          --sample {wildcards.sample} \
          --out {output.tsv}
        """

rule summarize_pass_alt_variants:
    input:
        vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",
        depth="results/{prefix}/stats/{sample}/{sample}_depth_of_coverage.sample_summary"
    output:
        tsv="results/{prefix}/pass_alt_variants/{sample}_pass_alt_variants.tsv"
    benchmark: 
        "benchmarks/{prefix}/summarize_pass_alt_variants/{sample}.benchmark.tsv"
    conda:
        "envs/cyvcf2.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        r"""
        set -euo pipefail
        python workflow/scripts/summarize_pass_variants_with_depth.py \
          --vcf {input.vcf} \
          --depth {input.depth} \
          --sample {wildcards.sample} \
          --out {output.tsv}
        """
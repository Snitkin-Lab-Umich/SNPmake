def consensus_vcf_list(wc):
    lst = [
        f"results/{wc.prefix}/unmapped_vcf/{wc.sample}/{wc.sample}_unmapped.vcf.gz",
        f"results/{wc.prefix}/filtered_mask_vcf/{wc.sample}/{wc.sample}_filtered_mask.vcf.gz",
        f"results/{wc.prefix}/indel_prox_mask_vcf/{wc.sample}/{wc.sample}_indel_prox_mask.vcf.gz",
    ]
    if PHAGE_ENABLED:
        lst.append(f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz")
    lst.append(f"results/{wc.prefix}/filtered_vcf/{wc.sample}/{wc.sample}_5bp_indel_removed.vcf.gz")
    return " ".join(lst)

# Make reference allele consensus fasta for each sample using the unmapped variants, filtered mask variants, and pass variants
rule consensus_ref_allele_unmapped_variants:
    input:
        ref=REF_GENOME,
        pass_vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",
        unmapped_vcf="results/{prefix}/unmapped_vcf/{sample}/{sample}_unmapped.vcf.gz",
        filtered_vcf="results/{prefix}/filtered_mask_vcf/{sample}/{sample}_filtered_mask.vcf.gz",                                 
        indel_prox_vcf="results/{prefix}/indel_prox_mask_vcf/{sample}/{sample}_indel_prox_mask.vcf.gz",
        phage_vcf_gz=lambda wc: f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz" if PHAGE_ENABLED else [],
        phage_tbi=lambda wc: f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz.tbi" if PHAGE_ENABLED else [],
    output:
        fasta="results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}/{sample}_ref_allele_unmapped_variants.fa"
    params:
        outdir="results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}",
        vcf_list=consensus_vcf_list
    singularity:
        "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/consensus_ref_allele_unmapped_variants/{sample}.benchmark.tsv"
    # envmodules:
    #     "Bioinformatics",
    #     "htslib",
    #     "bcftools"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        VCF_SAMPLE=$(bcftools query -l {input.pass_vcf} | head -n 1)

        # Build a single sample overlay VCF by concatenation (not merge)
        bcftools concat -a -D {params.vcf_list} -O z -o {params.outdir}/{wildcards.sample}_tmp_concat.vcf.gz
        bcftools sort -O z -o {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz {params.outdir}/{wildcards.sample}_tmp_concat.vcf.gz
        tabix -f -p vcf {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz

        # Drop het genotypes by setting them to missing
        bcftools +setGT {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz -O z \
        -o {params.outdir}/{wildcards.sample}_consensus_sources_homonly.vcf.gz \
        -- -t q -n . -i 'GT="het"'
        tabix -f -p vcf {params.outdir}/{wildcards.sample}_consensus_sources_homonly.vcf.gz

        # Now missing GT (formerly het) becomes N in the consensus
        bcftools consensus -s "$VCF_SAMPLE" -M N -f {input.ref} {params.outdir}/{wildcards.sample}_consensus_sources_homonly.vcf.gz > {output.fasta}

        # bcftools consensus -s "$VCF_SAMPLE" -f {input.ref} {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz > {output.fasta} 
        sed -i 's/>.*/>{wildcards.sample}/g' {output.fasta}
        """
rule make_msa:
    input:
        ref=REF_GENOME,
        cons=expand(
            "results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}/{sample}_ref_allele_unmapped_variants.dash.fa",
            sample=SAMPLES, prefix=PREFIX
        )
    output:
        msa="results/{prefix}/alignment/{prefix}_genome_aln_w_alt_allele_unmapped.fa"
    # benchmark:
    #     "benchmarks/{prefix}/make_msa/benchmark.tsv"
    threads: 2
    resources:
        mem_mb=1000,
        runtime=15
    shell:
        r"""
        set -euo pipefail

        # concatenate reference first, then all samples
        cat {input.ref} {input.cons} > {output.msa}
        """

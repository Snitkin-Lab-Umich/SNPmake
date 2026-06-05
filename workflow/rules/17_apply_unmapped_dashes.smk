rule apply_unmapped_dashes:
    input:
        fasta="results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}/{sample}_ref_allele_unmapped_variants.fa",
        unmapped="results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed_positions"
    output:
        fasta="results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}/{sample}_ref_allele_unmapped_variants.dash.fa"
    benchmark:
        "benchmarks/{prefix}/apply_unmapped_dashes/{sample}.benchmark.tsv"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        r"""
        set -euo pipefail
        python workflow/scripts/apply_unmapped_dashes.py \
          --fasta {input.fasta} \
          --unmapped {input.unmapped} \
          --out {output.fasta}
        """
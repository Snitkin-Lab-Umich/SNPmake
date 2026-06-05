def ref_id():
    return config.get("existing_msa_ref_id", config["reference_contig"])

rule merge_with_existing_msa:
    input:
        existing=lambda wc: config["existing_msa"],
        ref=REF_GENOME,
        add=expand(
            "results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}/{sample}_ref_allele_unmapped_variants.dash.fa",
            sample=SAMPLES, prefix=PREFIX
        )
    output:
        msa="results/{prefix}/alignment/{prefix}_merged_with_existing.fa"
    params:
        ref_id=ref_id(),
        dup_policy=config.get("existing_msa_duplicate_policy", "error")
    threads: 2
    resources:
        mem_mb=1000,
        runtime=15
    shell:
        r"""
        set -euo pipefail
        python workflow/scripts/merge_reference_msa.py \
          --existing_msa {input.existing} \
          --ref_fasta {input.ref} \
          --ref_id {params.ref_id} \
          --duplicate_policy {params.dup_policy} \
          --add {input.add} \
          --out {output.msa}
        """


def phage_arg(wc):
    if PHAGE_ENABLED:
        return f"--phage_mask_vcf results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz"
    return ""
##### PHAGE MASKING RULES #####
if PHAGE_ENABLED:
    rule make_phage_mask_vcf:
        input:
            ref=REF_GENOME,
            regions_json=lambda wc: os.path.join(PHASTEST_RESULTS, "predicted_phage_regions.json"),
            pass_vcf="results/{prefix}/{sample}/filtered_vcf/{sample}_5bp_indel_removed.vcf.gz", 
        output:
            vcf="results/{prefix}/phage_mask_vcf/{sample}/{sample}_phage_mask.vcf"
        params:
            chrom=config["reference_contig"]
        singularity:
            "docker://staphb/bcftools:1.23.1"
        benchmark: 
            "benchmarks/{prefix}/make_phage_mask_vcf/{sample}.benchmark.tsv"
        # envmodules:
        #     "Bioinformatics",
        #     "samtools",
        #     "bcftools"
        threads: 1
        resources:
            mem_mb=1000,
            runtime=10
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.vcf})

            VCF_SAMPLE=$(bcftools query -l {input.pass_vcf} | head -n 1)

            python workflow/scripts/make_phage_mask_vcf.py \
            --ref {input.ref} \
            --regions_json {input.regions_json} \
            --chrom {params.chrom} \
            --vcf_sample "$VCF_SAMPLE" \
            --out {output.vcf}
            """

    rule bgzip_tabix_phage_files:
        input:
            vcf="results/{prefix}/phage_mask_vcf/{sample}/{sample}_phage_mask.vcf"
        output:
            vcf_gz="results/{prefix}/phage_mask_vcf/{sample}/{sample}_phage_mask.vcf.gz",
            tbi="results/{prefix}/phage_mask_vcf/{sample}/{sample}_phage_mask.vcf.gz.tbi"
        singularity:
            "docker://staphb/htslib:1.23"
        benchmark: 
            "benchmarks/{prefix}/bgzip_tabix_phage_files/{sample}.benchmark.tsv"
        threads: 1
        resources:
            mem_mb=1000,
            runtime=10
        shell:
            r"""
            set -euo pipefail
            bgzip -f {input.vcf}
            tabix -f -p vcf {output.vcf_gz}
            """
    
    rule phastest_regions_tsv:
        input:
            regions_json=lambda wc: os.path.join(PHASTEST_RESULTS, "predicted_phage_regions.json")
        output:
            tsv="results/{prefix}/qc/phastest_phage_regions.tsv"
        # benchmark: 
        #     "benchmarks/{prefix}/phastest_regions_tsv/benchmark.tsv"
        threads: 1
        resources:
            mem_mb=1000,
            runtime=10
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv})

            python workflow/scripts/phage_regions.py \
            --regions_json {input.regions_json} \
            --out {output.tsv}
            """
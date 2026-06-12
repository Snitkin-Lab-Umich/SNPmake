# calling snps with samtools
# variant_calling
rule bcftools_call_snps:
    input:
        index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam",
    output:
        final_raw_vcf = temp("results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf"),
        raw_vcf_gz="results/{prefix}/samtools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf.gz",
    log:
        "logs/{prefix}/samtools_varcall/{sample}/{sample}_bcftools_call_snps.log"
    params:
        ref_genome = config["reference_genome"]
    singularity:
        "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/bcftools_call_snps/{sample}.benchmark.tsv"
    shell:
        """
        bcftools mpileup -f {params.ref_genome} {input.index_sorted_dups_rmvd_bam} | bcftools call -Ov -v -c -o {output.final_raw_vcf} &> {log}
        bgzip -c {output.final_raw_vcf} > {output.raw_vcf_gz}
        tabix -f -p vcf {output.raw_vcf_gz}
        """
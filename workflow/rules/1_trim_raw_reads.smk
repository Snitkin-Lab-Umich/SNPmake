rule trim_raw_reads:
    input:
        r1=input_r1,
        r2=input_r2
    output:
        r1="results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2="results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz",
        r1_unpaired="results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_unpaired.fastq.gz",
        r2_unpaired="results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_unpaired.fastq.gz"
    params:
        adapter_filepath=config["adaptor_filepath"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"]
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    benchmark: 
        "benchmarks/{prefix}/trim_raw_reads/{sample}.benchmark.tsv"
    threads: 6
    resources:
        mem_mb=6000,
        runtime=15
    shell:
        r"""
        set -euo pipefail
        trimmomatic PE \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          -threads {threads} \
          ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} \
          SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} \
          MINLEN:{params.minlength} \
          HEADCROP:{params.headcrop_length} \
          &> {log}
        """
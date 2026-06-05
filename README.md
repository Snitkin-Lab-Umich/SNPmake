# SNPmake - Microbial variant calling and reference-based alignment workflow

SNPmake is a modular Snakemake workflow for microbial variant calling, variant filtering, masking, and reference-based consensus alignment generation.

### Summary

As part of the SOP in the [Snitkin lab](https://github.com/Snitkin-Lab-Umich/Data-Flow-SOP), this pipeline can be used to perform reference based microbial SNP/indel calling from Illumina paired-end reads and generate a whole genome alignment suitable for downstream comparative genomic analyses.


In short, it performs the following steps:

* Accepts either raw or trimmed paired-end Illumina FASTQ reads.
* If `read_input_type: raw`, trims and filters low-quality bases and adapter sequences using [Trimmomatic](https://github.com/usadellab/Trimmomatic)(v0.39).
* If starting from raw reads, downsampling is performed after trimming using [seqtk](https://github.com/lh3/seqtk) to reduce very high coverage samples.
* Aligns reads to a user provided reference genome using [BWA-MEM](https://github.com/lh3/bwa).
* Removes clipped reads, converts SAM to BAM, sorts alignments, and indexes BAM files using [samtools](https://github.com/samtools/samtools).
* Removes PCR duplicates using Picard/samtools based post-alignment processing.
* Generates alignment statistics using [samtools flagstat](https://www.htslib.org/doc/samtools-flagstat.html).
* Calculates depth-of-coverage statistics using [GATK DepthOfCoverage](https://gatk.broadinstitute.org/).
* Identifies unmapped/zero coverage regions using [bedtools genomecov](https://bedtools.readthedocs.io/).
* Calls SNPs using [bcftools](https://samtools.github.io/bcftools/).
* Calls indels using [GATK](https://gatk.broadinstitute.org/).
* Applies hard filters to SNPs and indels using user-defined thresholds in `config/config.yaml`.
* Removes SNPs within 5 bp of indels.
* Annotates variants using [snpEff](https://pcingola.github.io/SnpEff/).
* Creates a reference based multi-FASTA alignment from all sample consensus sequences.
* Alternatively, the pipeline merges newly generated consensus sequences with an existing alignment, with strict reference sequence compatibility checks.

<!-- 
The workflow generates all output in the output prefix folder set in the config file. Instructions on setup are found [below](#setup-config-samples-and-profile-files). Each workflow step gets its own individual folder as shown. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and some of their contents.**

```
INSERT tree dir here
```
--> 

## Installation

> If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname.

```
cd /scratch/esnitkin_root/esnitkin1/your_uniqname/
```

> Clone the GitHub repository onto your system.
```
git clone https://github.com/Snitkin-Lab-Umich/SNPmake.git
```

> Ensure you have successfully cloned SNPmake. Type `ls` and you should see the newly created directory **_SNPmake_**. Move to the newly created directory.
```
cd SNPmake
```
> Load Bioinformatics, Snakemake, Singularity, and mamba modules from Great Lakes modules.
```
module load Bioinformatics snakemake singularity mamba
```

This workflow makes use of Singularity containers, conda environments, and Great Lakes environment modules. Several containers are available through the [StaPH-B Docker Builds](https://github.com/StaPH-B/docker-builds) project. If you are running SNPmake on Great Lakes, load the modules as shown above. If you are running it on another computing platform, ensure that Snakemake, Singularity/Apptainer, and conda/mamba are installed.

## Setup config and samples files

**_If you are just testing this pipeline, the config and sample files may already be loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix, reference genome, and input read paths in the config file to understand which variables need to be modified when running your own samples on SNPmake._**

### Config

SNPmake supports either raw reads or already trimmed reads as input. If raw reads are provided, the workflow performs read trimming and downsampling before alignment. If trimmed reads are provided (assuming it has been downsampled), reads are passed directly to the alignment step.

As input, the Snakemake workflow takes a config file where you can set the path to `config/samples.csv`, input FASTQ directory, input read type, reference genome, filtering thresholds, workflow mode, and optional existing MSA settings. Instructions on how to modify `config/config.yaml` are included in the config file.


### Samples

Add samples to `config/samples.csv`. The file should be a comma separated file containing at least one column:
sample_id

Example:
```
sample_id
IMPALA_001
IMPALA_002
IMPALA_003
```

`sample_id` should be the prefix extracted from your FASTQ reads.

For example, if your input read directory contains:
```
IMPALA_001_R1_trim_paired.fastq.gz
IMPALA_001_R2_trim_paired.fastq.gz
```

then your sample ID should be:
```IMPALA_001```

If using raw reads and your input read directory contains:
```
IMPALA_001_R1.fastq.gz
IMPALA_001_R2.fastq.gz
```

then your sample ID should still be:
```IMPALA_001```

You can create a `samples.csv` file for trimmed reads using the following loop. Replace `path_to_your_trimmed_reads` with the actual path to your trimmed FASTQ files.

```
echo "sample_id" > config/samples.csv

for read1 in path_to_your_trimmed_reads/*_R1_trim_paired.fastq.gz; do
    sample_id=$(basename "$read1" | sed 's/_R1_trim_paired.fastq.gz//g')
    echo "$sample_id"
done >> config/samples.csv

```

You can create a `samples.csv` file for raw reads using:

```
echo "sample_id" > config/samples.csv

for read1 in path_to_your_raw_reads/*_R1.fastq.gz; do
    sample_id=$(basename "$read1" | sed 's/_R1.fastq.gz//g')
    echo "$sample_id"
done >> config/samples.csv
```


## Quick start

### Preview SNPmake steps with a dry run

```
snakemake -s workflow/Snakefile --dryrun
```

### Run SNPmake locally

For small test datasets, SNPmake can be run locally:

```
snakemake -s workflow/Snakefile  --profile profile  --configfile config/config.yaml --cores all 
```

### Run SNPmake on Great Lakes using the Slurm profile as a batch job

Submitting the workflow as a batch job is recommended for full runs on Great Lakes.

Change these `SBATCH` commands before running:

* `--job-name` to a more descriptive name.
* `--mail-user` to your email address.
* `--time` depending on the number of samples you have. This is the runtime for the Snakemake controller job and should be longer than the total expected workflow runtime.
* `--account` to the appropriate Slurm account.
* `--partition` if needed.

Save the below script as something like `run_SNPmake.sbat`.
```
#!/bin/bash
#SBATCH --job-name=run_SNPmake
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard
set -euo pipefail

#Load necessary modules
module load Bioinformatics snakemake singularity mamba

#Run Snakemake pipeline
snakemake -s workflow/Snakefile  --profile profile --configfile config/config.yaml
```

Submit the batch job:
```sbatch run_SNPmake.sbat```

<!-- 
### Run only VCF generation

To run only the core VCF generation steps, set the following in `config/config.yaml`:
workflow_mode: vcf_only

Then run:
snakemake -s workflow/Snakefile --profile profile --configfile config/config.yaml

### Run full workflow

To generate VCFs, masks, consensus sequences, and the final alignment, set:
workflow_mode: full

Then run:
snakemake -s workflow/Snakefile --profile profile --configfile config/config.yaml

![Alt text](images/SNPkit2_dag.svg)
-->

## Dependencies

### Near Essential

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Conda](https://docs.conda.io/en/latest/)
* [mamba](https://mamba.readthedocs.io/en/latest/)

### Tool stack used in workflow

* [Trimmomatic](https://github.com/usadellab/Trimmomatic)
* [seqtk](https://github.com/lh3/seqtk)
* [BWA](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools/samtools)
* [bcftools](https://samtools.github.io/bcftools/)
* [GATK](https://gatk.broadinstitute.org/)
* [bedtools](https://bedtools.readthedocs.io/)
* [Picard](https://broadinstitute.github.io/picard/)
* [snpEff](https://pcingola.github.io/SnpEff/)
* [htslib](https://www.htslib.org/)
* [pigz](https://zlib.net/pigz/)
* [Pandas](https://pandas.pydata.org/)
* [Biopython](https://biopython.org/)



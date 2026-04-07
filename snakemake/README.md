# Snakemake Instructions

## About

This `snakemake` workflow automates the `long-read-rna-seq-pipeline` from start to finish. It has been tested and run with `slurm` as a job scheduler.

It takes `.bam` files from one or many samples as input. It returns an aggregate `.gtf` file and sample-specific transcript quantifications.

## Table of Contents

- [Snakemake Instructions](#snakemake-instructions)
  - [About](#about)
  - [Table of Contents](#table-of-contents)
  - [Environment Set Up](#environment-set-up)
    - [Additional Dependencies](#additional-dependencies)
  - [Configuration](#configuration)
    - [Output Directories](#output-directories)
    - [Asset Paths](#asset-paths)
    - [Resource Parameters](#resource-parameters)
    - [Samples](#samples)
    - [Other](#other)
  - [Example Usage](#example-usage)
  - [References](#references)

## Environment Set Up

To get started, you will need to have `conda` installed. See [install instructions here](<https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>) for more details.

Once `conda` is installed, you can make set up the custom environment for running the pipeline with [./install](./install). This will generate a new environment called `conda_env`.

### Additional Dependencies

The pipeline additionally requires a genome file in `fasta` format and a transcript annotation file in `gtf` format. Example human data can be downloaded from [GENCODE](<https://www.gencodegenes.org/>) as follows:

```{bash}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
```

## Configuration

The `snakemake` relies parameters specified in [snakemake_config.yaml](./snakemake_config.yaml).

### Output Directories

- `annotation_dir: '/path/to/final/gtf/output'`
  - Expected contents:
    - `gffcompare` output directory
    - `output.transcripts.gtf` (merged, sample-dependent transcript annotation)
    - `output.transcripts.updated.gtf.gz` (merged, sample-dependent transcript annotation with ORF calls)
    - `output.transcripts.updated.gtf.gz.tbi` (index for merged, sample-dependent transcript annotation with ORF calls)
- `ckpt_dir: '/path/to/checkpoint/output'`
  - *For now, this is used to facilitate parallel ORF calls by chromosome. This could probably be a temporary directory?*
- `log_dir: '/path/to/log/output'`
- `sample_dir: '/path/to/sample/output'`
  - Expected contents:
    - `output.gtf` (sample-specific transcript annotation)
    - `transcript_counts.tsv` (sample-specific transcript quantifications from the merged `output.transcripts.updated.gtf.gz`)

### Asset Paths

- `genome_path: '/path/to/genome/fasta'`
- `gencode_annotation_path: '/path/to/gencode/annotation'`
- `refTSS: '/path/to/refTSS/bed'`
- `polyA: '/path/to/polyA/bed'`

### Resource Parameters

- `{job}_threads: number of threads`
- `{job_mem_gb}: memory (GB)`
- `{job_time_hr}: time (hrs)`

### Samples

The pipeline takes `.bam` files as input. Samples are specified as follows:

```{markdown}
samples:
    sample_A:
        - bam: '/path/to/sample_A/bam'
    sample_B:
        - bam: '/path/to/sample_B/bam'
```

### Other

- `conda_wrapper: /path/to/conda_wrapper`
  - This field is added automatically during `./install`

## Example Usage

Once you've updated the [snakemake_config](./snakemake_config.yaml), you can run the `snakemake` with [./run](./run).

**Add toy-data and DAG example later*

## References

Scripts associated with the `snakemake_profile`, `conda` installation, and job submission were lifted from [ESPRESSO](<`https://github.com/Xinglab/espresso/tree/v1.6.0/snakemake`>).

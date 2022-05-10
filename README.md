# TBSeqPipe

## Introduction

TBSeqPipe is a flexible and user-friendly pipeline based on snakemake workflow for analyzing WGS data of *Mycobacterium tuberculosis* complex isolates. Taking illumina WGS data as input, this workflow preforms some basic analysis tasks as well as some downstream high-level analysis steps. TBSeqPipe generates a final summary report to better integrate and present results from all analysis modules.

## Workflow

![Workflow](/flowchart/flowchart.png)

## Installation

### Environment

#### Conda

Conda can function as a package manager and is available [here](https://docs.conda.io/en/latest/miniconda.html). If you have conda make sure the bioconda and conda-forge channels are added:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Snakemake

The Snakemake workflow management system is a tool to create reproducible and scalable data analyses. Detailed intsruction could be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Quick installation:

* Install mamba first (mamba provides a faster and more roboust way for conda packages installation):
```bash
conda install -n base -c conda-forge mamba
```
* Install snakemake using mamba:
```bash
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### Clone the repository

```bash
git clone git@github.com:KevinLYW366/TBSeqPipe.git
```

### Activate the environment

```bash
conda activate snakemake
```

### Kraken database

A pre-built 8 GB database [MiniKraken DB_8GB](https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz) is the suggested reference database for TBSeqPipe. It is constructed from complete bacterial, archaeal, and viral genomes in RefSeq.

## Set up configuration

To run the complete workflow do the following:

* Create an sample list file for all the samples you want to analyze with one ID per line.
* Copy all FASTQ files of your samples into one directory.
* Customize the workflow based on your need in: `config/configfile.yaml`. Parameters in "Required Parameters" section must be entered manually:
  * `sample_list`: `/path/to/sample_list_file`
  * `data_dir`: `/path/to/fastq_files`
  * `fastq_read_id_format`, `fastq_suffix_format` and `data_dir_format`: give values based on the FASTQ file directory structure and the format of FASTQ file names
  * `kraken_db`: `/path/to/minikraken_20171019_8GB`

## Usage

1. Move to the directory of TBSeqPipe.

```bash
cd /path/to/TBSeqPipe
```

2. A dry-run is recommended at first to check if everything is okay.

```bash
snakemake -r -p -n
```

3. If no error message shows up, let's do a formal run (feel free to modify "-j 40" which controls the CPU cores used in parallel).

```bash
snakemake --use-conda -r -p -j 40
```

### Note

#### Crashed and burned (Unlocking)

After the workflow was killed (Snakemake didn’t shutdown), the workflow directory will be still locked. If you are sure, that snakemake is no longer running `(ps aux | grep snake)`.

Unlock the working directory:

```bash
snakemake *.snakemake --unlock
```

#### Rerun incomplete

If Snakemake marked a file as incomplete after a crash, delete and produce it again.

```bash
snakemake *.snakemake --ri
```
 
## License

The code is available under the [GNU GPLv3 license](https://choosealicense.com/licenses/gpl-3.0/). The text and data are availabe under the [CC-BY license](https://choosealicense.com/licenses/cc-by-4.0/).

## Questions and Issues

For contacting the developer and issue reports please go to [Issues](https://github.com/KevinLYW366/TBSeqPipe/issues).

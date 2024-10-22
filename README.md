
# Nextflow example workflow

This repository contains a completed example of converting a simple [bash script](assets/workflow_bash.sh) into a nextflow workflow, made for the [2024 Nextflow vs Snakemake workshop](https://sagc-bioinformatics.github.io/nextflow-vs-snakemake-2024/).

The version of the workflow that was written on the day can be found in [main.nf.WORKSHOP_VERSION](main.nf.WORKSHOP_VERSION).

## Input data

Example data has been sourced from the [Snakemake tutorial](https://github.com/snakemake/snakemake-tutorial-data.git).

Data includes some short reads in fasta format and a portion of the *Saccharomyces cerevisiae* genome (if you're wondering why mapping statistics are so poor, it's incomplete).

## Configuration

The `nextflow.config` file is currently set up to use conda.

There is an `environment.yml` file under `assets/`. 

If you want, you can create a conda environment manually with `conda env create -f environment.yml`, and delete the config in `nextflow.config`.
This will actually make the pipeline execute slightly faster than allowing nextflow to manage the conda env for you.

Otherwise, no configuration is provided. The `local` executor will be used.

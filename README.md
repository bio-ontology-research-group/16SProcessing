# 16SProcessing
A Nextflow workflow to run vsearch for short read 16S data

## Before running:
- Install Nextflow
- Install Singularity
- Change path of data folder in nextflow.config runOptions

## To run:
nextflow run maintest.nf --in_dir test_16S_reads -with-singularity BORG_16S_3.sif

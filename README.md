# 16SProcessing
A Nextflow workflow to run vsearch for short read 16S data

## Before running:
- Install Nextflow
- Install Singularity
- Create 'data' folder in main directory and download and unzip:
    * https://www.drive5.com/sintax/rdp_16s_v18.fa.gz
    * https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip
- Change path of data folder in nextflow.config runOptions

## To run:
nextflow run maintest.nf --in_dir test_16S_reads -with-singularity BORG_16S.sif

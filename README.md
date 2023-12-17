# 16SProcessing
A Nextflow workflow to run vsearch for short read 16S data. The workflow includes the following steps:
1. Trimming (cutadapt) and QC (fastQC)
2. Merging and filtering reads (vsearch)
3. Clustering and chimera removal (vsearch and swarm)
4. OTU classification (vsearch sintax)

#### Notes:
- The workflow is configured to work with docker or singularity. The singularity profile works with SLURM by default, an sbatch job can be submitted with the available example script.
- Fastq files must be named *_L001_R{1,2}_001.fastq.gz

## Before running:
- Install Nextflow
- Install Singularity or Docker
- Move to 'data' folder in main directory and download and unzip:
    * https://www.drive5.com/sintax/rdp_16s_v18.fa.gz
    * https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip
    * Run unzip -p silva.gold.bacteria.zip | sed -e "s/[.-]//g" > gold.fasta
    * Run gunzip rdp_16s_v18.fa.gz
- Change path of data folder in nextflow.config runOptions

## To run:
nextflow run 16SProcessing.nf --in_dir directory/with/fastq/files -profile (docker OR singularity)

##### Workflow can be tested with:
nextflow run 16SProcessing.nf --in_dir test_16S_reads -profile (docker OR singularity)

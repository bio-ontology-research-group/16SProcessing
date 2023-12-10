#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BORG/vsearch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : 
----------------------------------------------------------------------------------------
*/

params.in_dir = './' // Default input directory
params.fwd_primer = 'CCTACGGGNGGCWGCAG' // Default forward primer sequence
params.rev_primer = 'GACTACHVGGGTATCTAATCC' // Default reverse primer sequence
params.db_path = "${baseDir}/data/gold.fasta"
params.rdp_path = "${baseDir}/data/rdp_16s_v18.fa"
/*
    ================================================================================
                                Preprocessing and QC
    ================================================================================
*/


process CutAdapt {
    tag "${pair_id}"
    publishDir "results/trimmed_reads", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    def (read1, read2) = reads
    """
    cutadapt -g ${params.fwd_primer} -G ${params.rev_primer} \
        -o ${pair_id}_L001_R1_001_trimmed.fastq.gz \
        -p ${pair_id}_L001_R2_001_trimmed.fastq.gz \
        ${read1} ${read2}
    """
}

process FastQC {
    tag "${pair_id}"
    publishDir "results/fastQC", mode: 'copy'

    input:
    tuple val(pair_id), path(read_files)

    output:
    path "*.html", emit: fastqc_html
    path "*.zip", emit: fastqc_zip

    script:
    def (read1, read2) = read_files
    """
    fastqc ${read1} ${read2}
    """
}

/*
    ================================================================================
                                        Merging
    ================================================================================
*/

process VsearchMerge {
    tag "${pair_id}"
    publishDir "results/intermediate_files", mode: 'copy'

    input:
    tuple val(pair_id), path(trimmed_reads)

    output:
    tuple val(pair_id), path("*.merged.fastq"), emit: merged_fastq

    script:
    def (read1, read2) = trimmed_reads
    """
    vsearch --fastq_mergepairs ${read1} --reverse ${read2} --fastqout ${pair_id}.merged.fastq
    """
}

process VsearchStats {
    tag "${pair_id}"
    publishDir "results/intermediate_files", mode: 'copy'

    input:
    tuple val(pair_id), path(merged_fastq)

    output:
    path "${pair_id}.stats"

    script:
    """
    vsearch --fastq_eestats ${merged_fastq} --output ${pair_id}.stats
    """
}

process VsearchFilter {
    tag "${pair_id}"
    publishDir "results/intermediate_files", mode: 'copy'

    input:
    tuple val (pair_id), path(merged_fastq)

    output:
    tuple val(pair_id), path("*.filtered.fasta"), emit: filtered_fasta

    script:
    """
    vsearch --fastq_filter ${merged_fastq} --fastq_maxee 1.0 --fastq_minlen 10 --fastq_maxlen 500 --fastq_maxns 0 --fastaout ${pair_id}.filtered.fasta --fasta_width 0
    """
}

process VsearchDereplicate {
    tag "${pair_id}"
    publishDir "results/intermediate_files", mode: 'copy'

    input:
    tuple val(pair_id), path(filtered_fasta)

    output:
    path "${pair_id}.derep.fasta"

    script:
    """
    vsearch --derep_fulllength ${filtered_fasta} --strand plus --output ${pair_id}.derep.fasta --sizeout --relabel ${pair_id} --fasta_width 0
    """
}

process MergeAll{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path derep_fasta_files

    output:
    path "all.fasta"

    script:
    """
    cat ${derep_fasta_files.join(' ')} > all.fasta
    """
}

process VsearchDerepAll{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path all_merged_fasta

    output:
    path "derep.fasta"

    script:
    """
    vsearch --derep_fulllength ${all_merged_fasta} --sizein --sizeout --fasta_width 0 --uc all.derep.uc --output derep.fasta
    """
}

/*
    ================================================================================
                                        Clustering
    ================================================================================
*/

process VsearchCluster{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path all_derep_fasta

    output:
    path "centroids.fasta"

    script:
    """
    vsearch --cluster_size derep.fasta --id 0.98 --strand plus --sizein --sizeout --fasta_width 0 --centroids centroids.fasta
    """
}

process SwarmCluster{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path all_derep_fasta
    path centroids

    output:
    path "centroids.fasta"

    script:
    """
    swarm ${all_derep_fasta} --differences 1 --fastidious --seeds ${centroids} --usearch-abundance --output /dev/null
    """
}

process SingletonRemoval{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path swarm_centroids

    output:
    path "sorted.fasta"

    script:
    """
    vsearch --sortbysize ${swarm_centroids} --sizein --sizeout --fasta_width 0 --minsize 2 --output sorted.fasta
    """
}

process DenovoChimeraRemoval{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path nonsingleton_clusters

    output:
    path "denovo.nonchimeras.fasta"

    script:
    """
    vsearch --uchime_denovo ${nonsingleton_clusters} --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras denovo.nonchimeras.fasta
    """
}

process ReferenceChimeraRemoval{
    publishDir "results/merged_intermediates", mode: 'copy'

    input:
    path denovo_nonchimeras

    output:
    path "nonchimeras.fasta"

    script:
    """
    vsearch --uchime_ref ${denovo_nonchimeras} --db ${params.db_path} --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --nonchimeras nonchimeras.fasta
    """
}

/*
    ================================================================================
                                        OTUs
    ================================================================================
*/

process RelabelOTUs{
    publishDir "results/OTU_tables", mode: 'copy'

    input:
    path ref_nonchimeras

    output:
    path "otus.fasta"

    script:
    """
    vsearch --fastx_filter ${ref_nonchimeras} --sizein --sizeout --fasta_width 0 --relabel OTU_ --fastaout otus.fasta
    """
}

process MapOTUs{
    publishDir "results/OTU_tables", mode: 'copy'

    input:
    path all_merged_fasta
    path relabeled_otus

    output:
    path "otutab.txt"

    script:
    """
    vsearch --usearch_global ${all_merged_fasta} --db ${relabeled_otus} --id 0.98 \
        --strand plus --sizein --sizeout --fasta_width 0 --qmask none --dbmask none \
        --otutabout otutab.txt 
    """
}

process ClassifyOTUs{
    publishDir "results/OTU_tables", mode: 'copy'

    input:
    path relabeled_otus

    output:
    path "reads.sintax"

    script:
    """
    vsearch -sintax ${relabeled_otus} -db ${params.rdp_path} -tabbedout reads.sintax \
        -strand both --sintax_cutoff 0.8
    """
}

// process MapOTUClassification{
//     publishDir "results/OTU_tables", mode: 'copy'

//     input:
//     path sintax_classification

//     output:

//     script:
//     """

//     """
// }


workflow {
    read_pairs = Channel.fromFilePairs("${params.in_dir}/*_L001_R{1,2}_001.fastq.gz", size: 2)
    CutAdapt(read_pairs)
    FastQC(CutAdapt.out.trimmed_reads)

    //Merge
    vsearch_merge_out = VsearchMerge(CutAdapt.out)
    VsearchStats(vsearch_merge_out)
    vsearch_filter_out = VsearchFilter(vsearch_merge_out)
    derep_out = VsearchDereplicate(vsearch_filter_out)
    derep_fasta_files = derep_out.collect()
    all_merged_fasta = MergeAll(derep_fasta_files)
    all_derep_fasta = VsearchDerepAll(all_merged_fasta)

    //Cluster
    centroids = VsearchCluster(all_derep_fasta)
    swarm_centroids = SwarmCluster(all_derep_fasta, centroids)
    nonsingleton_clusters = SingletonRemoval(swarm_centroids)
    denovo_nonchimeras = DenovoChimeraRemoval(nonsingleton_clusters)
    ref_nonchimeras = ReferenceChimeraRemoval(denovo_nonchimeras)

    //OTUs
    relabeled_otus = RelabelOTUs(ref_nonchimeras)
    otu_table = MapOTUs(all_merged_fasta, relabeled_otus)
    sintax_classification = ClassifyOTUs(relabeled_otus) 
}
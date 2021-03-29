#!/usr/bin/env nextflow

// Nextflow.configuration
aux=config.aux_location
data=config.data_location // data_location or btdata_location
output=config.output_location
genome_dir=config.genome_location

large_core=config.large_core
small_core=config.small_core

// Parameters

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"


////////////////////////////////////////////////
// ** - Pull in fq files (paired)
////////////////////////////////////////////////

fqs = Channel.fromPath(data + "${params.dir}/*.f[a-z]*q.gz")
                        .map { n -> [ n.getName(), n ] }


////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trim_reads {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.json'

   input:
       tuple val(id), file(reads) from fqs

   output:
       tuple id_out, file("${id_out}.fq.gz") into trimmed_fqs
       tuple file("*.html"), file("*.json")  into trim_log

  script:
      id_out = id.replace('.fastq.gz', '')

   """
       fastp -i $reads -o ${id_out}.fq.gz -y -l 15 -h ${id_out}.html -j ${id_out}.json
   """
}
trimmed_fqs.set { trimmed_reads_bwa }


////////////////////////////////////////////////
// ** - Fetch genome
////////////////////////////////////////////////

hosturl="https://www.vectorbase.org/download/aedes-albopictus-foshanscaffoldsaalof1fagz"

process fetch_ref {

    publishDir "${output}/reference/", mode: 'copy'

    output:
        file("reference.fa.gz") into reference_fa

    """
        echo '${hosturl}'
        wget ${hosturl} -O reference.fa.gz
    """
}
reference_fa.into { bwa_index }


////////////////////////////////////////////////
// ** - Index Genome (bwa)
////////////////////////////////////////////////

process build_bwa_index {

    publishDir "${output}/reference/", mode: 'copy'

    cpus large_core

    input:
        file("reference.fa.gz") from bwa_index

    output:
        file "reference.*" into bwa_indices

    """
        zcat reference.fa.gz > reference.fa
        bwa index reference.fa
    """
}


////////////////////////////////////////////////
// ** - bwa mapping
////////////////////////////////////////////////

process bwa_align {
    publishDir "${output}/bwa_stats/", mode: 'copy', pattern: '*align_.txt'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(reads) from trimmed_reads_bwa
        file bwa_indices from bwa_indices.first()

    output:
        file("${id}_align.txt") into bwa_stats
        file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indexes

    script:
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
        index_base = bwa_indices[0].toString() - ~/.fa[.a-z]*/

        """
        bwa aln -o 0 -n 0 -t ${large_core} ${index_base}.fa ${reads} > ${id}.sai
        bwa samse ${index_base}.fa ${id}.sai ${reads} > ${id}.sam
        samtools view -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -b ${id}.bam
        samtools flagstat ${id}.bam > ${id}_align.txt
        """
}

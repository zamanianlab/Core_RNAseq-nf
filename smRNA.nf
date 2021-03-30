#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
aux=params.aux
genome=params.genome

big=params.big
small=params.small

// Global Params
params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"

// flag for fastqc and multiqc (--qc)
params.qc = false


////////////////////////////////////////////////
// ** - Pull in fq files
////////////////////////////////////////////////

fqs = Channel.fromPath(input + "${params.dir}/*.f[a-z]*q.gz")
                        .map { n -> [ n.getName(), n ] }


////////////////////////////////////////////////
// ** - Trim reads
////////////////////////////////////////////////

process trim_reads {

  publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.{json,html}'

  cpus small
  tag { id }

  input:
    tuple val(id), file(reads) from fqs

  output:
    tuple id_out, file("${id_out}.fq.gz") into trimmed_fqs
    tuple file("*.html"), file("*.json")  into trim_log

  script:
      id_out = id.replace('.fastq.gz', '')


  """
     fastp -i $reads -a AACTGTAGGCACCATCAAT -o ${id_out}.fq.gz -y -l 15 -h ${id_out}.html -j ${id_out}.json
  """
}
trimmed_fqs.into { trimmed_reads_bwa; trimmed_reads_qc }




////////////////////////////////////////////////
// ** - Fetch viral genomes
////////////////////////////////////////////////

Ae_bangkok_url="https://vectorbase.org/common/downloads/Current_Release/AalbopictusFPA/fasta/data/VectorBase-51_AalbopictusFPA_Genome.fasta"
Ae_mo_url="https://vectorbase.org/common/downloads/Current_Release/AalbopictusFPA/gff/data/VectorBase-51_AalbopictusFPA.gff"

process fetch_ref {

    publishDir "${genome}/reference/", mode: 'copy'

    output:
        file("reference.fa") into reference_fa

    """
        echo '${genome_url}'
        wget ${genome_url} -O reference.fa
        echo '${annot_url}'
        wget ${annot_url} -O geneset.gff
    """
}
reference_fa.into { bwa_index }


////////////////////////////////////////////////
// ** - Index Genome (bwa)
////////////////////////////////////////////////

process build_bwa_index {

    cpus huge

    input:
        file("reference.fa") from bwa_index

    output:
        file "reference.*" into bwa_indices

    """
        bwa index reference.fa
    """
}


////////////////////////////////////////////////
// ** - bwa mapping
////////////////////////////////////////////////

process bwa_align {
    publishDir "${output}/${params.dir}/bwa_stats/", mode: 'copy', pattern: '*align_.txt'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus big
    tag { id }

    input:
        tuple val(id), file(reads) from trimmed_reads_bwa
        file bwa_indices from bwa_indices.first()

    output:
        file("${id}_align.txt") into bwa_stats
        file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indexes

    script:
        fa_prefix = reads[0].toString() - ~/(\.fq\.gz)$/
        index_base = bwa_indices[0].toString() - ~/.fa[.a-z]*/

        """
        bwa aln -o 0 -n 0 -t ${task.cpus} ${index_base}.fa ${reads} > ${id}.sai
        bwa samse ${index_base}.fa ${id}.sai ${reads} > ${id}.sam
        samtools view -@ ${task.cpus} -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${task.cpus} -m 8G -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -@ ${task.cpus} -b ${id}.bam
        samtools flagstat ${id}.bam > ${id}_align.txt
        """
}

////////////////////////////////////////////////
// ** - Fetch genome
////////////////////////////////////////////////

genome_url="https://vectorbase.org/common/downloads/Current_Release/AalbopictusFPA/fasta/data/VectorBase-51_AalbopictusFPA_Genome.fasta"
annot_url="https://vectorbase.org/common/downloads/Current_Release/AalbopictusFPA/gff/data/VectorBase-51_AalbopictusFPA.gff"

process fetch_ref {

    publishDir "${genome}/reference/", mode: 'copy'

    output:
        file("reference.fa") into reference_fa

    """
        echo '${genome_url}'
        wget ${genome_url} -O reference.fa
        echo '${annot_url}'
        wget ${annot_url} -O geneset.gff
    """
}
reference_fa.into { bwa_index }


////////////////////////////////////////////////
// ** - Index Genome (bwa)
////////////////////////////////////////////////

process build_bwa_index {

    cpus huge

    input:
        file("reference.fa") from bwa_index

    output:
        file "reference.*" into bwa_indices

    """
        bwa index reference.fa
    """
}

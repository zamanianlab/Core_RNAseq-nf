#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
aux=params.aux
genome=params.genome

huge=params.huge
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

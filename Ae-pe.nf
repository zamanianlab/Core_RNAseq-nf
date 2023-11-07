#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
aux=params.aux

big=params.big
small=params.small

// Global Parameters

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "prjn: $params.rlen"


////////////////////////////////////////////////
// ** - Pull in fq files (paired) and indexed genome
////////////////////////////////////////////////

Channel.fromFilePairs(input + "/${params.dir}/*_{1,2}.fq.gz", flat: true)
        .set { fq_pairs }

star_indices = Channel.fromPath(output + "/Aeaeg_index/*" )

////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trim_reads {

   cpus  big
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.json'

   input:
       tuple val(id), file(forward), file(reverse) from fq_pairs

  output:
//    tuple id, file("${id}_R1.fq.gz"), file("${id}_R2.fq.gz") into trimmed_fqs
    tuple id, file("${id}_1.fq.gz"), file("${id}_2.fq.gz") into trimmed_fqs	
    tuple file("*.html"), file("*.json")  into trim_log

  """
	fastp -i $forward -I $reverse -w ${task.cpus} -o ${id}_trimmed_1.fq.gz -O ${id}_trimmed_2.fq.gz -y -l 150 -h ${id}.html -j ${id}.json	
  """
}
trimmed_fq_pairs.set { trimmed_reads_star; trimmed_reads_qc }

////////////////////////////////////////////////
// ** Align reads using STAR
////////////////////////////////////////////////

process star_align {

    publishDir "${output}/${params.dir}/star", mode: 'copy', pattern: '*.Log.final.out'
    publishDir "${output}/${params.dir}/star", mode: 'copy', pattern: '*.flagstat.txt'
    publishDir "${output}/${params.dir}/counts", mode: 'copy', pattern: '*.ReadsPerGene.tab'
    //publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    //publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus big
    tag { id }
    maxForks 6


    input:
        file("STAR_index/*") from star_indices
        tuple val(id), file(forward), file(reverse) from trimmed_reads_star

    output:
        tuple file("${id}.Log.final.out"), file("${id}.flagstat.txt") into alignment_logs_star
        tuple id, file("${id}.bam"), file("${id}.bam.bai") into bam_files_star
        file("${id}.ReadsPerGene.tab") into star_counts

    script:

        """
          STAR --runThreadN ${task.cpus} --runMode alignReads --genomeDir STAR_index\
            --outSAMtype BAM Unsorted --readFilesCommand zcat \
            --outFileNamePrefix ${id}. --readFilesIn ${forward} ${reverse}\
            --peOverlapNbasesMin 10 \
            --quantMode GeneCounts --outSAMattrRGline ID:${id}
          samtools sort -@ ${task.cpus} -m 24G -o ${id}.bam ${id}.Aligned.out.bam
          rm *.Aligned.out.bam
          samtools index -@ ${task.cpus} -b ${id}.bam
          samtools flagstat ${id}.bam > ${id}.flagstat.txt
          cat ${id}.ReadsPerGene.out.tab | cut -f 1,2 > ${id}.ReadsPerGene.tab
        """
// remove -m 12G
}
bam_files_star.into {bam_files_qc}

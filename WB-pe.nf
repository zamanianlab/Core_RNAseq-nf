#!/usr/bin/env nextflow

// Params from config files (system-dependent)

data=params.data
output=params.output
aux=params.aux

println "data: $data"

large_core=params.large_core
small_core=params.small_core

// Global Params

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

params.release = null
if( !params.release ) error "Missing release parameter"
println "release: $params.release"

params.species = null
if( !params.species ) error "Missing species parameter"
println "species: $params.species"

params.prjn = null
if( !params.species ) error "Missing prjn parameter"
println "prjn: $params.prjn"

// flag for final stringtie_table_counts process (--stc)
params.stc = false

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"

//params.rtype = null
//if( !params.dir ) error "Missing read type parameter (single = SE; paired = PE)"
//println "rtype: $params.rtype"


////////////////////////////////////////////////
// ** - Pull in fq files (paired)
////////////////////////////////////////////////

Channel.fromFilePairs(data + "${params.dir}/*_R{1,2}_001.f[a-z]*q.gz", flat: true)
        .set { fq_pairs }
        .println()


////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trim_reads {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.json'

   input:
       tuple val(id), file(forward), file(reverse) from fq_pairs

   output:
       tuple id, file("${id}_R1.fq.gz"), file("${id}_R2.fq.gz") into trimmed_fq_pairs
       tuple file("*.html"), file("*.json")  into trim_log

   """
       fastp -i $forward -I $reverse -o ${id}_R1.fq.gz -O ${id}_R2.fq.gz -y -l 50 -h ${id}.html -j ${id}.json
   """
}
trimmed_fq_pairs.set { trimmed_reads_hisat }


////////////////////////////////////////////////
// ** - Fetch genome (fa.gz) and gene annotation file (gtf.gz)
////////////////////////////////////////////////

process fetch_genome {

    output:
        file("geneset.gtf.gz") into geneset_gtf
        file("reference.fa.gz") into reference_fa

    script:

        prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${params.release}/species/${params.species}/${params.prjn}"

    """
        echo '${prefix}'
        curl ${prefix}/${params.species}.${params.prjn}.${params.release}.canonical_geneset.gtf.gz > geneset.gtf.gz
        curl ${prefix}/${params.species}.${params.prjn}.${params.release}.genomic.fa.gz > reference.fa.gz
    """
}
geneset_gtf.into { geneset_hisat; geneset_stringtie }
reference_fa.into { reference_hisat }


////////////////////////////////////////////////
// ** - HiSat2/Stringtie pipeline
////////////////////////////////////////////////

// Create HiSat2 Index using reference genome and annotation file
extract_exons = file("${aux}/scripts/hisat2_extract_exons.py")
extract_splice = file("${aux}/scripts/hisat2_extract_splice_sites.py")

process build_hisat_index {

    cpus large_core

    input:
        file("geneset.gtf.gz") from geneset_hisat
        file("reference.fa.gz") from reference_hisat

    output:
        file "*.ht2" into hs2_indices

    """
        zcat geneset.gtf.gz | python ${extract_splice} - > splice.ss
        zcat geneset.gtf.gz | python ${extract_exons} - > exon.exon
        zcat reference.fa.gz > reference.fa
        hisat2-build -p ${large_core} --ss splice.ss --exon exon.exon reference.fa reference.hisat2_index
    """

}

// Alignment and stringtie combined
process hisat2_stringtie {

    publishDir "${output}/expression", mode: 'copy', pattern: '**/*'
    publishDir "${output}/expression", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(forward), file(reverse) from trimmed_reads_hisat
        file("geneset.gtf.gz") from geneset_stringtie
        file hs2_indices from hs2_indices.first()

    output:
        file "${id}.hisat2_log.txt" into alignment_logs
        file("${id}/*") into stringtie_exp
        file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indexes

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
        hisat2 -p ${large_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
        samtools view -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -b ${id}.bam
        zcat geneset.gtf.gz > geneset.gtf
        stringtie ${id}.bam -p ${large_core} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf
        rm *.gtf
    """
}

////////////////////////////////////////////////
// ** - STRINGTIE table counts
////////////////////////////////////////////////

prepDE = file("${aux}/scripts/prepDE.py")
process stringtie_counts_final {

    echo true

    publishDir "${output}/counts", mode: 'copy', pattern: '*.csv'

    cpus small_core

    when:
      params.stc

    output:
        file ("gene_count_matrix.csv") into gene_count_matrix
        file ("transcript_count_matrix.csv") into transcript_count_matrix

    """
        python ${prepDE} -i ${output}/expression -l ${rlen} -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}

#!/usr/bin/env nextflow

// Params from config files (system-dependent)

data=params.data
output=params.output
aux=params.aux

large_core=params.large_core
small_core=params.small_core

// Global Params

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

params.release = null
if( !params.release ) error "Missing WB release parameter"
println "release: $params.release"

params.species = null
if( !params.species ) error "Missing WB species parameter"
println "species: $params.species"

params.prjn = null
if( !params.prjn ) error "Missing WB prjn parameter"
println "prjn: $params.prjn"

params.se = false
println "se: $params.se"

// flags for final stringtie_table_counts process (--stc)
params.stc = false

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"


////////////////////////////////////////////////
// ** - Pull in fq files (paired vs unpaired)
////////////////////////////////////////////////

if ( !params.se ) {
  Channel.fromFilePairs(data + "${params.dir}/*_R{1,2}_001.f[a-z]*q.gz", flat: true)
          .set { fqs }
} else if ( params.se ) {
  fqs = Channel.fromPath(data + "${params.dir}/*.f[a-z]*q.gz")
                          .map { n -> [ n.getName(), n ] }
} else exit 1, 'error loading fqs'


////////////////////////////////////////////////
// ** TRIM READS (SE or PE)
////////////////////////////////////////////////

process trim_reads_pe {

   cpus small_core
   tag { id }
   publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.json'

   when:
       !params.se

   input:
       tuple val(id), file(forward), file(reverse) from fqs

   output:
       tuple id, file("${id}_R1.fq.gz"), file("${id}_R2.fq.gz") into trimmed_fq_pairs
       tuple file("*.html"), file("*.json")  into trim_log

   """
       fastp -i $forward -I $reverse -w ${small_core} -o ${id}_R1.fq.gz -O ${id}_R2.fq.gz -y -l 50 -h ${id}.html -j ${id}.json
   """
}
trimmed_fqs.set { trimmed_reads_hisat }

process trim_reads_se {

   cpus small_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.json'

   when:
       params.se

   input:
       tuple val(id), file(reads) from fqs

   output:
       tuple id_out, file("${id_out}.fq.gz") into trimmed_fqs
       tuple file("*.html"), file("*.json")  into trim_log

  script:
      id_out = id.replace('.fastq.gz', '')

   """
       fastp -i $reads -o ${id_out}.fq.gz -y -l 50 -h ${id_out}.html -j ${id_out}.json
   """
}
trimmed_fqs.into { trimmed_reads_hisat }


////////////////////////////////////////////////
// ** - Fetch genome (fa.gz) and gene annotation file (gtf.gz)
////////////////////////////////////////////////

process fetch_genome {

    cpus small_core

    output:
        file("geneset.gtf.gz") into geneset_gtf
        file("reference.fa.gz") into reference_fa

    script:

        prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${params.release}/species/${params.species}/${params.prjn}"

    """
        echo '${prefix}'
        wget -c ${prefix}/${params.species}.${params.prjn}.${params.release}.canonical_geneset.gtf.gz -O geneset.gtf.gz
        wget -c ${prefix}/${params.species}.${params.prjn}.${params.release}.genomic.fa.gz -O reference.fa.gz
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

// Alignment and stringtie combined (handles both SE and PE)
process hisat2_stringtie {

    publishDir "${output}/${params.dir}/expression", mode: 'copy', pattern: '**/*'
    publishDir "${output}/${params.dir}/expression", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

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

        if (rtype == "PE")
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
        else if (rtype == "SE")
            """
            hisat2 -p ${large_core} -x $index_base -U ${reads} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
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
        else
            """
            """
}

////////////////////////////////////////////////
// ** - STRINGTIE table counts
////////////////////////////////////////////////

prepDE = file("${aux}/scripts/prepDE.py")
process stringtie_counts_final {

    echo true

    publishDir "${output}/${params.dir}/counts", mode: 'copy', pattern: '*.csv'

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

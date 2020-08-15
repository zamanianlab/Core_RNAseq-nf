#!/usr/bin/env nextflow

// Params from config files (system-dependent)

data=params.data // data = btc seq, data2 = uploaded seq
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

// flags for final stringtie_table_counts process (--stc)
params.stc = false

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"


////////////////////////////////////////////////
// ** - Pull in fq files
////////////////////////////////////////////////

fqs = Channel.fromPath(data + "/${params.dir}/*.f[a-z]*q.gz")
                        .map { n -> [ n.getName(), n ] }


////////////////////////////////////////////////
// ** - Trim reads
////////////////////////////////////////////////

process trim_reads_se {

  publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.{json,html}'

  cpus small_core
  tag { id }

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
trimmed_fqs.into { trimmed_reads_hisat; trimmed_reads_qc; geneset_qc }


////////////////////////////////////////////////
// ** - multiQC of trimmed fqs
////////////////////////////////////////////////

process fastqc {

    publishDir "${output}/${params.dir}/fastqc", mode: 'copy', pattern: '*_fastqc.{zip,html}'

    cpus small_core
    tag { id }

    input:
    tuple val(id), file(reads) from trimmed_reads_qc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:

    """
      fastqc -q $reads
    """
}

process multiqc {
  publishDir "${output}/${params.dir}/fastqc", mode: 'copy', pattern: 'multiqc_report.html'

  cpus small_core
  tag { id }

    input:
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report

    script:

    """
      multiqc .
    """
}

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

// Alignment and stringtie
process hisat2_stringtie {

    publishDir "${output}/${params.dir}/expression", mode: 'copy', pattern: '**/*'
    publishDir "${output}/${params.dir}/expression", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(reads) from trimmed_reads_hisat
        file("geneset.gtf.gz") from geneset_stringtie
        file hs2_indices from hs2_indices.first()

    output:
        file "${id}.hisat2_log.txt" into alignment_logs
        file("${id}/*") into stringtie_exp
        tuple id, file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indexes

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

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

}


////////////////////////////////////////////////
// ** - Stringtie table counts
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
        python2 ${prepDE} -i ${output}/${params.dir}/expression -l ${params.rlen} -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}


////////////////////////////////////////////////
// ** - Post-alignment QC
////////////////////////////////////////////////

process align_analysis {

    cpus small_core

}

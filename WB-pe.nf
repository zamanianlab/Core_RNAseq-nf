#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
aux=params.aux

big=params.big
small=params.small

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

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"

// flag for fastqc and multiqc (--qc)
params.qc = false


////////////////////////////////////////////////
// ** - Pull in fq files
////////////////////////////////////////////////

Channel.fromFilePairs(input + "/${params.dir}/*_R{1,2}_001.f[a-z]*q.gz", flat: true)
          .set { fqs }


////////////////////////////////////////////////
// ** - Trim reads
////////////////////////////////////////////////

process trim_reads {

  publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.{json,html}'

  cpus small
  tag { id }

  input:
    tuple val(id), file(forward), file(reverse) from fqs

  output:
    tuple id, file("${id}_R1.fq.gz"), file("${id}_R2.fq.gz") into trimmed_fqs
    tuple file("*.html"), file("*.json")  into trim_log

  """
    fastp -i $forward -I $reverse -w ${task.cpus} -o ${id}_R1.fq.gz -O ${id}_R2.fq.gz -y -l 50 -h ${id}.html -j ${id}.json
  """
}
trimmed_fqs.into { trimmed_reads_hisat; trimmed_reads_star; trimmed_reads_qc }


////////////////////////////////////////////////
// ** - multiQC of trimmed fqs
////////////////////////////////////////////////

process fastqc {

    publishDir "${output}/${params.dir}/fastqc", mode: 'copy', pattern: '*_fastqc.{zip,html}'

    cpus small
    tag { id }

    when:
      params.qc

    input:
    tuple val(id), file(forward), file(reverse) from trimmed_reads_qc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:

    """
      fastqc -q $forward $reverse -t ${task.cpus}
    """
}

process multiqc {

    publishDir "${output}/${params.dir}/fastqc", mode: 'copy', pattern: 'multiqc_report.html'

    cpus small

    when:
      params.qc

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

    cpus small

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
geneset_gtf.into { geneset_hisat; geneset_stringtie; geneset_star; geneset_qc }
reference_fa.into { reference_hisat; reference_star }


////////////////////////////////////////////////
// ** - STAR pipeline
////////////////////////////////////////////////

// Build STAR Index using reference genome and annotation file
process star_index {

    cpus big

    input:
        file("geneset.gtf.gz") from geneset_star
        file("reference.fa.gz") from reference_star

    output:
        file("STAR_index/*") into star_indices

    script:
        overhang = params.rlen - 1

    """
        zcat reference.fa.gz > reference.fa
        zcat geneset.gtf.gz > geneset.gtf
        mkdir STAR_index

        STAR --runThreadN ${task.cpus} --runMode genomeGenerate  --genomeDir STAR_index \
          --genomeFastaFiles reference.fa \
          --sjdbGTFfile geneset.gtf \
          --sjdbOverhang ${overhang}
    """

}

process star_align {

    publishDir "${output}/${params.dir}/star", mode: 'copy', pattern: '*.Log.final.out'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus big
    tag { id }
    maxForks 4

    input:
        file("STAR_index/*") from star_indices
        tuple val(id), file(forward), file(reverse) from trimmed_reads_star

    output:
        file("${id}.Log.final.out") into alignment_logs_star
        tuple id, file("${id}.bam"), file("${id}.bam.bai") into bam_files_star

    script:

        """
          STAR --runThreadN ${task.cpus} --runMode alignReads --outSAMtype BAM Unsorted \
            --readFilesCommand zcat --genomeDir STAR_index \
            --outFileNamePrefix ${id}. --readFilesIn  ${forward} ${reverse}
          samtools flagstat ${id}.Aligned.out.bam
          samtools sort -@ ${task.cpus} -m 16G -o ${id}.bam ${id}.Aligned.out.bam
          rm *.Aligned.out.bam
          samtools index -@ ${task.cpus} -b ${id}.bam
        """
// remove -m16G
}


// STAR --runMode alignReads --outSAMtype BAM Unsorted --readFilesCommand zcat --genomeDir /path/to/STAR/genome/folder --outFileNamePrefix {sample name}  --readFilesIn  /path/to/R1 /path/to/R2

////////////////////////////////////////////////
// ** - HiSat2/Stringtie pipeline
////////////////////////////////////////////////

// Build HiSat2 Index using reference genome and annotation file
extract_exons = file("${aux}/scripts/hisat2_extract_exons.py")
extract_splice = file("${aux}/scripts/hisat2_extract_splice_sites.py")

process hisat_index {

    cpus big

    input:
        file("geneset.gtf.gz") from geneset_hisat
        file("reference.fa.gz") from reference_hisat

    output:
        file "*.ht2" into hs2_indices

    """
        zcat geneset.gtf.gz | python ${extract_splice} - > splice.ss
        zcat geneset.gtf.gz | python ${extract_exons} - > exon.exon
        zcat reference.fa.gz > reference.fa
        hisat2-build -p ${task.cpus} --ss splice.ss --exon exon.exon reference.fa reference.hisat2_index
    """

}

// HiSat2 alignment
process hisat_align {

    publishDir "${output}/${params.dir}/hisat", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/${params.dir}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus big
    tag { id }
    maxForks 4

    input:
        tuple val(id), file(forward), file(reverse) from trimmed_reads_hisat
        file hs2_indices from hs2_indices.first()

    output:
        file "${id}.hisat2_log.txt" into alignment_logs_hisat
        tuple id, file("${id}.bam"), file("${id}.bam.bai") into bam_files_hisat

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

        """
          hisat2 -p ${task.cpus} -x $index_base -1 ${forward} -2 ${reverse} --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" --summary-file ${id}.hisat2_log.txt | \
             samtools view -@ ${task.cpus} -bS > ${id}.unsorted.bam
          samtools flagstat ${id}.unsorted.bam
          samtools sort -@ ${task.cpus} -m 16G -o ${id}.bam ${id}.unsorted.bam
          rm *.unsorted.bam
          samtools index -@ ${task.cpus} -b ${id}.bam
        """

}
bam_files_hisat.into {bam_files_stringtie; bam_files_htseq; bam_files_qc}

// Stringtie
process stringtie {

    publishDir "${output}/${params.dir}/hisat", mode: 'copy', pattern: '**/*'

    cpus small
    tag { id }
    maxForks 4

    input:
        file("geneset.gtf.gz") from geneset_stringtie
        tuple val(id), ("${id}.bam"), file("${id}.bam.bai") from bam_files_stringtie

    output:
        file("${id}/*") into stringtie_exp

    script:

        """
          zcat geneset.gtf.gz > geneset.gtf
          stringtie ${id}.bam -p ${task.cpus} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf
          rm *.gtf
        """

}


// Stringtie table counts [collect stringtie outputs to confirm stringtie is complete before generating counts]
prepDE = file("${aux}/scripts/prepDE.py")
process stringtie_counts {

    publishDir "${output}/${params.dir}/counts", mode: 'copy', pattern: '*.csv'

    cpus small

    input:
      val(id) from stringtie_exp.toSortedList()

    output:
      file ("gene_count_matrix.csv") into gene_count_matrix
      file ("transcript_count_matrix.csv") into transcript_count_matrix

    """
      python2 ${prepDE} -i ${output}/${params.dir}/hisat -l ${params.rlen} -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}

////////////////////////////////////////////////
// ** - Post-alignment QC
////////////////////////////////////////////////

process align_analysis {

    publishDir "${output}/${params.dir}/align_qc", mode: 'copy', pattern: '*_QC.txt'

    cpus small

    when:
      params.qc

    input:
        file("geneset.gtf.gz") from geneset_qc
        tuple val(id), file(bam), file(bai) from bam_files_qc

    output:
        file("*_QC.txt") into align_qc

    script:

    """
      zcat geneset.gtf.gz > geneset.gtf
      awk '{ if (\$0 ~ "transcript_id") print \$0; else print \$0" transcript_id "";"; }' geneset.gtf | gtf2bed - > geneset.gtf.bed
      cat geneset.gtf.bed | sed '/\tgene\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print \$1, \$2, \$3, \$4, \$6, \$8}' > geneset.gene.bed
      cat geneset.gtf.bed | sed '/\texon\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print \$1, \$2, \$3, \$4, \$6, \$8}' > geneset.exon.bed
      cat geneset.gtf.bed | sed '/\tfive_prime_utr\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print \$1, \$2, \$3, \$4, \$6, \$8}' > geneset.5utr.bed
      cat geneset.gtf.bed | sed '/\tthree_prime_utr\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print \$1, \$2, \$3, \$4, \$6, \$8}' > geneset.3utr.bed

      echo -n "total," >> ${bam}_QC.txt
      samtools view -c ${bam} >> ${bam}_QC.txt
      echo -n "mapped," >> ${bam}_QC.txt
      samtools view -F 0x4 -c ${bam} >> ${bam}_QC.txt
      echo -n "unique," >> ${bam}_QC.txt
      samtools view -F 0x4 -q 60 -c ${bam} >> ${bam}_QC.txt

      samtools view -L geneset.gene.bed -h ${bam} > tmp.sam
      echo -n "gene," >> ${bam}_QC.txt
      samtools view -F 0x4 -q 60 -c tmp.sam >> ${bam}_QC.txt
      rm tmp.sam
      samtools view -L geneset.exon.bed -h ${bam} > tmp.sam
      echo -n "exon," >> ${bam}_QC.txt
      samtools view -F 0x4 -q 60 -c tmp.sam >> ${bam}_QC.txt
      rm tmp.sam
      samtools view -L geneset.5utr.bed -h ${bam} > tmp.sam
      echo -n "5utr," >> ${bam}_QC.txt
      samtools view -F 0x4 -q 60 -c tmp.sam >> ${bam}_QC.txt
      rm tmp.sam
      samtools view -L geneset.3utr.bed -h ${bam} > tmp.sam
      echo -n "3utr," >> ${bam}_QC.txt
      samtools view -F 0x4 -q 60 -c tmp.sam >> ${bam}_QC.txt
      rm tmp.sam
    """
}

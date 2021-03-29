#!/usr/bin/env nextflow

// Nextflow.configuration
aux=config.aux_location
data=config.data_location
output=config.output_location
aedesgenome=config.aedesgenome_location

large_core=config.large_core
small_core=config.small_core

// Parameters

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

////////////////////////////////////////////////
// ** - Pull in fq files (paired)
////////////////////////////////////////////////

fqs = Channel.fromPath(data + "${params.dir}/*.fastq.gz")
                        .map { n -> [ n.getName(), n ] }


// rRNAs = file(GHdata + "smRNA/rRNA/ascaris_suum_rRNA.fasta")
// tRNAs = file(GHdata + "smRNA/tRNA/ascaris_suum_tRNA.fasta")
bm_miRNAs_mature = file(GHdata + "smRNA/miRNA/brugia_malayi_mature_b.fasta")
bm_miRNAs_prec = file(GHdata + "smRNA/miRNA/brugia_malayi_stemloop_b.fasta")
ce_miRNAs_mature = file(GHdata + "smRNA/miRNA/caenorhabditis_elegans_mature_b.fasta")
ce_miRNAs_prec = file(GHdata + "smRNA/miRNA/caenorhabditis_elegans_stemloop_b.fasta")
ae_miRNAs_mature = file(GHdata + "smRNA/miRNA/aedes_aegypti_mature_b.fasta")
ae_miRNAs_prec = file(GHdata + "smRNA/miRNA/aedes_aegypti_stemloop_b.fasta")

// ** - Fetch reference genome (fa.gz)
release="WBPS11"
species="brugia_malayi"
prjn="PRJNA10729"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

process fetch_parasite_ref {

    publishDir "${output}/reference/", mode: 'copy'

    output:
        file("parasite.fa.gz") into parasite_ref

    """
        echo '${prefix}'
        curl ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz > parasite.fa.gz
    """
}
parasite_ref.into { parasite_bwa; parasite_bowtie; parasite_mirdeep }

// ** - Fetch host genome (fa.gz)
hosturl="https://www.vectorbase.org/download/aedes-aegypti-lvpagwgchromosomesaaegl5fagz"

process fetch_host_ref {

    publishDir "${output}/reference/", mode: 'copy'

    output:
        file("host.fa.gz") into host_ref

    """
        echo '${hosturl}'
        wget ${hosturl} -O host.fa.gz
    """
}
host_ref.into { host_bwa; host_bowtie; host_mirdeep }


//** TRIM READS
process trimmomatic {
    cpus small_core
    tag { name }

    publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*_trimout.txt'

    input:
        set val(name), file(reads) from fq_set

    output:
        set val(name_out), file(name_out) into fq_trim
        file("*_trimout.txt") into trim_log

    script:
    name_out = name.replace('.fastq.gz', '_trim.fq.gz')

    """
        trimmomatic SE -phred33 -threads ${large_core} ${reads} ${name_out} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 &> ${reads}_trimout.txt
    """
}
fq_trim.into { fq_trim_bwa; fq_trim_bowtie; fq_trim_contam; fq_trim_mirdeepMAP_P; fq_trim_mirdeepMIR_P ; fq_trim_mirdeepMAP_H; fq_trim_mirdeepMIR_H}



//INDEX PARASITE GENOME - BOWTIE
process bowtie_index_parasite {

    publishDir "${output}/reference/", mode: 'copy'

    cpus large_core

    input:
        file("parasite.fa.gz") from parasite_bowtie

    output:
        file "parasite_bowtie*.ebwt" into parasite_bowtie_indices

    script:

    """
        zcat parasite.fa.gz > parasite.fa
        bowtie-build parasite.fa parasite_bowtie
    """
}


// Mirdeep2 mapper.pl (parasite genome)
process mirDeep2_mapper_parasite {
    cpus large_core
    tag { name }

    input:
        set val(name), file(reads) from fq_trim_mirdeepMAP_P
        file bowtieindex from parasite_bowtie_indices.first()

    output:
        file("${fa_prefix}_parasite_map.arf") into reads_vs_parasite_genome_arf
        file("${fa_prefix}_parasite_collapsed.fa") into reads_parasite_collapsed

    script:
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/

        """
        zcat ${reads} > ${fa_prefix}.fa
        mapper.pl ${fa_prefix}.fa -e -h -j -l 18 -m -p parasite_bowtie -s ${fa_prefix}_parasite_collapsed.fa -t ${fa_prefix}_parasite_map.arf -v
        """
}
reads_parasite_collapsed.into {reads_parasite_collapsed_Q; reads_parasite_collapsed_M}


// Mirdeep2 quantifier.pl (map to predefined parasite mature/precursor seqs)
process quantifier_pl_parasite {

    publishDir "${output}/quantifier_parasite/${fa_prefix}/", mode: 'copy'

    cpus large_core
    tag { collapsed_reads }

    input:
        file collapsed_reads from reads_parasite_collapsed_Q

    output:
        file "*" into quantifier_paraste_out

    script:
        fa_prefix = collapsed_reads[0].toString() - ~/(_parasite_collapsed)(\.fa)$/

        """
        quantifier.pl -p ${bm_miRNAs_prec} -m ${bm_miRNAs_mature} -r ${collapsed_reads} -y now
        """
}


//INDEX HOST GENOMES - BOWTIE
process bowtie_index_host {

    publishDir "${output}/reference/", mode: 'copy'

    cpus large_core

    input:
        file("host.fa.gz") from host_bowtie

    output:
        file "host_bowtie*.ebwt" into host_bowtie_indices

    script:

    """
        zcat host.fa.gz > host.fa
        bowtie-build host.fa host_bowtie
    """
}

// Mirdeep2 mapper.pl (host genome)
process mirDeep2_mapper_host {
    cpus large_core
    tag { name }

    input:
        set val(name), file(reads) from fq_trim_mirdeepMAP_H
        file bowtieindex from host_bowtie_indices.first()

    output:
        file("${fa_prefix}_host_map.arf") into reads_vs_host_genome_arf
        file("${fa_prefix}_host_collapsed.fa") into reads_host_collapsed

    script:
        fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/

        """
        zcat ${reads} > ${fa_prefix}.fa
        mapper.pl ${fa_prefix}.fa -e -h -j -l 18 -m -p host_bowtie -s ${fa_prefix}_host_collapsed.fa -t ${fa_prefix}_host_map.arf -v
        """
}
reads_host_collapsed.into {reads_host_collapsed_Q; reads_host_collapsed_M}


// Mirdeep2 quantifier.pl (map to predefined host mature/precursor seqs)
process quantifier_pl_host {

    publishDir "${output}/quantifier_host/${fa_prefix}/", mode: 'copy'

    cpus large_core
    tag { collapsed_reads }

    input:
        file collapsed_reads from reads_host_collapsed_Q

    output:
        file "*" into quantifier_host_out

    script:
        fa_prefix = collapsed_reads[0].toString() - ~/(_host_collapsed)(\.fa)$/

        """
        quantifier.pl -p ${ae_miRNAs_prec} -m ${ae_miRNAs_mature} -r ${collapsed_reads} -y now
        """
}


// // Mirdeep2 quantifier.pl (map to predefined host mature/precursor seqs)
// process quantifier_pl_host {
//
//     publishDir "${output}/quantifier_host/", mode: 'copy'
//     cpus large_core
//     tag { reads }
//
//     input:
//         set val(id), file(reads) from fq_trim_mirdeepQ2
//
//     script:
//         fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
//
//         """
//         zcat ${reads} > ${fa_prefix}.fa
//         quantifier.pl -p ${ae_miRNAs_prec} -m ${ae_miRNAs_mature} -r ${fa_prefix}.fa -y now
//         """
// }



//
//
// // Mirdeep2 mirdeep2.pl
// process mirDeep2_pl {
//     cpus large_core
//     tag { reads }
//
//     input:
//         file("parasite.fa.gz") from parasite_mirdeep
//         file reads_vs_parasite_genome_arf from reads_vs_parasite_genome_arf
//         file reads_parasite_collapsed from reads_parasite_collapsed
//
//         """
//         zcat parasite.fa.gz > parasite.fa
//         cat parasite.fa | awk '{print \$1}' > parasite_temp.fa
//         miRDeep2.pl ${reads_parasite_collapsed} parasite_temp.fa ${reads_vs_parasite_genome_arf} ${bm_miRNAs_mature} ${ce_miRNAs_mature} ${bm_miRNAs_prec} -P
//         """
// }


// //INDEX GENOMES - BWA
// process build_bwa_index {
//
//     publishDir "${output}/reference/", mode: 'copy'
//
//     cpus large_core
//
//     input:
//         file("parasite.fa.gz") from parasite_bwa
//         file("host.fa.gz") from host_bwa
//
//     output:
//         file "parasite.*" into bwa_parasite_indices
//         file "host.*" into bwa_host_indices
//
//     """
//         zcat parasite.fa.gz > parasite.fa
//         bwa index parasite.fa
//         zcat host.fa.gz > host.fa
//         bwa index host.fa
//     """
// }

// ALIGN TRIMMED READS TO PARASITE GENOME (BWA)
// process align {
//     publishDir "${output}/bwa_stats/", mode: 'copy'
//
//     cpus large_core
//     tag { id }
//
//     input:
//         set val(id), file(reads) from fq_trim1
//         file(parasite_bwaindex) from bwa_parasite_indices.first()
//
//     output:
//         file("bwa_parasite_align.txt") into bwa_stats
//
//     script:
//         fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
//
//         """
//         bwa aln -o 0 -n 0 -t ${large_core} parasite.fa ${reads} > ${id}.sai
//         bwa samse parasite.fa ${id}.sai ${reads} > ${id}.sam
//         samtools view -bS ${id}.sam > ${id}.unsorted.bam
//         rm *.sam
//         samtools flagstat ${id}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
//         rm *.unsorted.bam
//         samtools index -b ${id}.bam
//         samtools flagstat ${id}.bam > bwa_parasite_align.txt
//         """
// }

// // Map rRNAs and tRNAs
// process map_rRNAs_tRNAs {
//     publishDir "${output}/stats/", mode: 'copy'
//
//     cpus large_core
//     tag { reads }
//
//     input:
//         file reads from fq_trim2
//         file rRNA_fa from rRNAs
//         file tRNA_fa from tRNAs
//
//     output:
//         file("bwa_rRNA_align.txt") into bwa_rRNA_alignstats
//         file("bwa_tRNA_align.txt") into bwa_tRNA_alignstats
//
//     script:
//         fa_prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
//
//         """
//         bwa index ${rRNA_fa}
//
//         bwa aln -o 0 -n 0 -t ${large_core} ${rRNA_fa} ${reads} > ${fa_prefix}.sai
//         bwa samse ${rRNA_fa} ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
//         samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
//         samtools flagstat ${fa_prefix}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${fa_prefix}_rRNA.bam ${fa_prefix}.unsorted.bam
//         samtools index -b ${fa_prefix}_rRNA.bam
//         samtools flagstat ${fa_prefix}_rRNA.bam > bwa_rRNA_align.txt
//
//         bwa index ${tRNA_fa}
//
//         bwa aln -o 0 -n 0 -t ${large_core} ${tRNA_fa} ${reads} > ${fa_prefix}.sai
//         bwa samse ${tRNA_fa} ${fa_prefix}.sai ${reads} > ${fa_prefix}.sam
//         samtools view -bS ${fa_prefix}.sam > ${fa_prefix}.unsorted.bam
//         samtools flagstat ${fa_prefix}.unsorted.bam
//         samtools sort -@ ${large_core} -o ${fa_prefix}_tRNA.bam ${fa_prefix}.unsorted.bam
//         samtools index -b ${fa_prefix}_tRNA.bam
//         samtools flagstat ${fa_prefix}_tRNA.bam > bwa_tRNA_align.txt
//
//         """
// }

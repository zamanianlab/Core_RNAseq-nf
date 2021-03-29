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
trimmed_fqs.set { trimmed_reads_mirdeep }


////////////////////////////////////////////////
// ** - Load in genome and indices
////////////////////////////////////////////////

geneset_gtf = file("${genome_dir}/Other/Aedes_aegypti/annotation/geneset_h.gtf.gz")
genome_fa = file("${genome_dir}/Other/Aedes_aegypti/genome.fa")
bowtie2_indices = Channel.fromPath("${genome_dir}/Other/Aedes_aegypti/bowtie2Index/*").collect()

aae_mature = file(aux + "mirbase/aae_mature.fa")
aae_prec = file(aux + "mirbase/aae_pre.fa")


////////////////////////////////////////////////
// ** - mirDeep2 pipeline
////////////////////////////////////////////////

// Mirdeep2 mapper.pl (map to genome)
process mirDeep2_mapper {

    cpus small_core
    tag { id }

    input:
        tuple val(id), file(reads) from trimmed_reads_mirdeep
        file bowtie2_indices from bowtie2_indices

    output:
        file("${id}_map.arf") into reads_vs_genome_arf
        tuple val(id), file("${id}_collapsed.fa") into reads_collapsed

    script:
        index_base = bowtie2_indices[0].toString() - ~/./

      """
        zcat ${reads} > ${id}.fa
        mapper.pl ${id}.fa -e -h -j -l 18 -m -p ${index_base} -s ${id}_collapsed.fa -t ${id}_map.arf -v
      """
}
reads_collapsed.into {reads_collapsed_Q; reads_collapsed_M}


// Mirdeep2 quantifier.pl (map to predefined mature/precursor seqs)
process mirDeep2_quantifier {

    publishDir "${output}/quantifier/${id}/", mode: 'copy'

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(collapsed_reads) from reads_collapsed_Q

    output:
        file "*" into quantifier_out

    """
        quantifier.pl -p ${aae_prec} -m ${aae_mature} -r ${collapsed_reads} -y now
    """
}

// Mirdeep2 mirdeep2.pl
process mirDeep2_pl {

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(collapsed_reads) from reads_collapsed_M
        file reads_vs_genome_arf from reads_vs_genome_arf
        file("genome.fa") from genome_fa

        """
        cat genome.fa | awk '{print \$1}' > genome_clean.fa
        miRDeep2.pl ${collapsed_reads} genome_clean.fa ${reads_vs_genome_arf} ${aae_mature} none ${aae_prec} -P
        """
}

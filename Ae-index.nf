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



////////////////////////////////////////////////
// ** - Fetch genome and gene annotation files
////////////////////////////////////////////////

genome_url="https://vectorbase.org/common/downloads/Current_Release/AaegyptiLVP_AGWG/fasta/data/VectorBase-65_AaegyptiLVP_AGWG_Genome.fasta"
annot_url="https://vectorbase.org/common/downloads/Current_Release/AaegyptiLVP_AGWG/gff/data/VectorBase-65_AaegyptiLVP_AGWG.gff"

process fetch_ref {

    publishDir "${output}/${params.dir}/", mode: 'copy'

    output:
        file("reference.fa") into reference_fa
        file("geneset.gff") into geneset_gff

    """
        echo '${genome_url}'
        wget ${genome_url} -O reference.fa
        echo '${annot_url}'
        wget ${annot_url} -O geneset.gff
    """
}
reference_fa.into { star_index }


geneset_gff.into { geneset_star}
reference_fa.into { reference_star }

////////////////////////////////////////////////
// ** - Index Genome (STAR)
////////////////////////////////////////////////
/*
process build_star_index {

    cpus huge

    publishDir "${output}/${params.dir}/", mode: 'copy'

    input:
        file("reference.fa") from star_index

    output:
        file "reference.*" into star_indices

    """
        bwa index reference.fa
    """
}
*/

// Build STAR Index using reference genome and annotation file
process star_index {

    cpus huge

    publishDir "${output}/${params.dir}/", mode: 'copy'

    input:
        file("geneset.gff") from geneset_star
        file("reference.fa") from reference_star

    output:
        file("STAR_index/*") into star_indices

    script:
        overhang = params.rlen - 1

    """
        mkdir STAR_index

        STAR --runThreadN ${task.cpus} --runMode genomeGenerate  --genomeDir STAR_index \
          --genomeFastaFiles reference.fa \
          --sjdbGTFfile geneset.gff \
          --sjdbOverhang ${overhang}
    """

}
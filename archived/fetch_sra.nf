#!/usr/bin/env nextflow

// Edit nextflow.configuration!

aux=config.aux_location
data=config.data_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// ** - Get txt file of SRA accession IDs from 'auxillary' folder
sra_file = Channel.fromPath(aux + "SRR_Acc_List_UGA.txt")
sra_file.into { sra_file_fetch; sra_file_convert }

// ** - Download SRA files based on text file list of SRA accession IDs (goes to ncbi folder)
// process fetch_SRA {
    
//     cpus small_core

//     input:
//         file("SRR_Acc_List_UGA.txt") from sra_file_fetch

//     script:

//     sra_list="SRR_Acc_List_UGA.txt"

//     """ 

//     while read line     
//     do           
//         echo \$line
//         prefetch \$line 
//     done <${sra_list} 

//     """
// }

// ** - Covert SRA files to fastqs (comment out until fetch_SRA{} complete)
process sra_to_fastq {

    cpus small_core

    publishDir "${data}/fq/", mode: 'move'
    
    input:
        file("SRR_Acc_List_UGA.txt") from sra_file_convert

    output:
        file("*")

    script:

    sra_list="SRR_Acc_List_UGA.txt"

    """ 

    while read line     
    do           
        fastq-dump --gzip --split-files ~/ncbi/public/sra/\$line.sra
    done <${sra_list} 

    """

}

// Can remove the sra files from NCBI folder when finished
# WB_RNAseq-nf
nf RNA-seq pipelines for Wormbase Parasite and Vectorbase species

## Contents

### Nextflow scripts
- WB-se.nf (WBP species single-end reads)
- WB-pe.nf (WBP species paired-end reads)
- Ae-pe.nf (Aedes aegypti paired-end reads)

### Nextflow config files
- chtc.config (UW CHTC server in Docker environment)
- chtc-local.config (local CHTC troubleshooting in Docker environment)
- brc.config (UW BRC server)

## Output dirs and files
- /trim_stats/: read trimming log files
- /fastqc/: post-trimming fastqc and multiqc outputs
- /expression/: HISAT2 outputs and log files
- /bams/: sorted .bam files
- /counts/: Stringtie gene and transcript count matrices
- /align_qc/: read distribution stats
- /trace/: Nextflow timeline, report, and trace files

## Example of core command (WBP)
- WBP pipelines download and index required genomes.

`nextflow run -resume [WB-pe.nf|WB-se.nf] [nextflow options] --dir [fastq dir] --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729" --rlen "100"`


## Example (Ae. aegypti)
- VP pipelines use pre-indexed genomes.

`nextflow run Ae-pe.nf [nextflow options] --dir [fastq dir] --rlen "100"`

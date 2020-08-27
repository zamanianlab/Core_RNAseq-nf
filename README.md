# WB_RNAseq-nf
nf RNA-seq pipelines for Wormbase Parasite and Vectorbase species

## Main Nextflow scripts
- WB-se.nf (WBP species single-end reads)
- WB-pe.nf (WBP species paired-end reads)
- Ae-pe.nf (Aedes aegypti paired-end reads)

## Output dirs and files
- /trim_stats/: read trimming log files
- /fastqc/: post-trimming fastqc and multiqc outputs
- /expression/: HISAT2 outputs and log files
- /bams/: sorted .bam files
- /counts/: Stringtie gene and transcript count matrices
- /align_qc/: read distribution stats
- /trace/: Nextflow timeline, report, and trace files

## Example of core command (WBP)
- Trim reads, download and index WBP reference genome, and align cleaned reads to reference (HISAT2)

`nextflow run -resume [WB-pe.nf|WB-se.nf] [nextflow options] --dir [fastq dir] --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729" --rlen "100"`




## Example (Ae. aegypti)

- Step 1: trim reads, download and index reference genome, and align cleaned reads to reference (HISAT2)

`nextflow run Ae-pe.nf [nextflow options] --dir [fastq dir] --rlen "100"`

- Step 2: re-run with `-resume` and append `--stc` flag to get stringtie count tables

`nextflow run -resume Ae-pe.nf [nextflow options] --dir [fastq dir] --rlen "100" --stc`

## Other

- fetch_sra.nf (fetches SRA data using .txt file list of SRA accession ID)

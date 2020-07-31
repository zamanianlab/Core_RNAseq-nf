# WB_RNAseq-nf
nf RNA-seq pipelines for Wormbase Parasite and Vectorbase species

## Main Nextflow scripts
- WB-se.nf (WBP species single-end reads)
- WB-pe.nf (WBP species paired-end reads)
- Ae-pe.nf (Aedes aegypti paired-end reads)

## Output dirs and files
- /trim_stats/:  read trimming log files
- /counts/:  "gene_count_matrix.csv"; "transcript_count_matrix.csv" from Stringtie
- /expression/: HISAT2 outputs and log files
- /bams/: sorted .bam files

## Example (WBP)

- Step 1: trim reads, download and index reference genome, and align cleaned reads to reference (HISAT2)

`nextflow run [WB-pe.nf|WB-se.nf] [nextflow options] --dir [fastq dir] --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729" --rlen "100"`

- Step 2: re-run with `-resume` and append `--stc` flag to get stringtie count tables

`nextflow run -resume [WB-pe.nf|WB-se.nf] [nextflow options] --dir [fastq dir] --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729" --rlen "100" --stc`

## Example (Ae. aegypti)

- Step 1: trim reads, download and index reference genome, and align cleaned reads to reference (HISAT2)

`nextflow run Ae-pe.nf [nextflow options] --dir [fastq dir] --rlen "100"`

- Step 2: re-run with `-resume` and append `--stc` flag to get stringtie count tables

`nextflow run -resume Ae-pe.nf [nextflow options] --dir [fastq dir] --rlen "100" --stc`

## Other

- fetch_sra.nf (fetches SRA data using .txt file list of SRA accession ID)

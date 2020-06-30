# WB_RNAseq-nf
nf RNA-seq pipeline for WB/WBP species

## Main Nextflow scripts
- main-se.nf (single-end reads)
- main-pe.nf (paired-end reads)

## Output
- /trim_stats/ (read trimming log files)
- /counts/ ("gene_count_matrix.csv" & "transcript_count_matrix.csv" from Stringtie)
- /expression/ (HISAT2 outputs and log files)
- /bams/ (sorted .bam files)
  
## Example

- Step 1: trim reads (fastp or Trimmomatic), download and index reference genome, and align cleaned reads to reference (HISAT2)

`nextflow run [main-pe.nf|main-se.nf] [nextflow options] --dir [fastq dir] --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729"`

- Step 2: re-run with `-resume` and append `--stc` flag to get stringtie count tables

`nextflow run -resume main-pe.nf [nextflow options] --dir [fastq dir] --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729" --stc`

## Other

- fetch_sra.nf (fetches SRA data using .txt file list of SRA accession ID)

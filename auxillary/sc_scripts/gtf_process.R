#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

#read in gtf
gtf <- read.table("geneset.gtf", sep = "\t", header = FALSE, quote = "")

#read in contig lengths
contig.len <- read.table("contig_lengths.txt", sep = "\t", header = FALSE, quote = "")
colnames(contig.len) <- c("seqname","len")

#read in 3' UTR extension value
ext <- as.numeric(args[1])

#apply extension to gene, transcript, last exon, and 3' UTR
gtf.3ext <- gtf %>%
  mutate(transcript_id = str_match(V9, "transcript_id \"(.*?)\";")[,2]) %>%
  mutate(exon_number = str_match(V9, "exon_number \"(\\d+)\";")[,2]) %>%
  mutate(exon_number = replace_na(exon_number, as.character(0))) %>% 
  mutate(exon_number = as.numeric(exon_number)) %>% 
  group_by(transcript_id) %>% 
  mutate(group_max = max(exon_number)) %>% 
  mutate(group_max = ifelse(exon_number == group_max & V3 == "exon", "yes", "no")) %>% 
  ungroup() %>% 
  mutate(V5 = ifelse(group_max == "yes" & V7 == "+", as.numeric(V5) + ext, V5)) %>% 
  mutate(V4 = ifelse(group_max == "yes" & V7 == "-", as.numeric(V4) - ext, V4)) %>% 
  mutate(V4 = ifelse(group_max == "yes" & V4 < 1, 1, V4)) %>%
  mutate(V4 = ifelse(V3 == "three_prime_utr" & V7 =="-", as.numeric(V4) - ext, V4)) %>%
  mutate(V4 = ifelse(V3 == "three_prime_utr" & V4 < 1, 1, V4)) %>%
  mutate(V5 = ifelse(V3 == "transcript" & V7 =="+", as.numeric(V5) + ext, V5)) %>%
  mutate(V4 = ifelse(V3 == "transcript" & V7 =="-", as.numeric(V4) - ext, V4)) %>%
  mutate(V4 = ifelse(V3 == "transcript" & V4 < 1, 1, V4)) %>%
  mutate(V5 = ifelse(V3 == "gene" & V7 =="+", as.numeric(V5) + ext, V5)) %>%
  mutate(V4 = ifelse(V3 == "gene" & V7 =="-", as.numeric(V4) - ext, V4)) %>%
  mutate(V4 = ifelse(V3 == "gene" & V4 < 1, 1, V4))
colnames(gtf.3ext) <- c("seqname","source","feature","start","end","score","strand","frame","attribute", "transcript_id", "exon_number", "group_max")

#truncate if exceeds boundaries of contig 
gtf.3ext <- left_join(gtf.3ext,contig.len, by = "seqname") %>%
  mutate(end = ifelse(strand == "+" & end >= len, len, end)) %>%
  select(-len, -group_max, -exon_number, -transcript_id)

#generate new gtf
write.table(gtf.3ext,"geneset.3ext.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

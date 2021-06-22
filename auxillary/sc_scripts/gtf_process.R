#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

gtf <- read.table("geneset.gtf", sep = "\t", header = FALSE, quote = "")
ext <- as.numeric(args[1])

gtf.3ext <- gtf %>%
  mutate(V5 = ifelse(V3 == "three_prime_utr" & V7 =="+", as.numeric(V5) + ext, V5)) %>%
  mutate(V4 = ifelse(V3 == "three_prime_utr" & V7 =="-", as.numeric(V4) - ext, V4)) %>%
  mutate(V4 = ifelse(V3 == "three_prime_utr" & V4 < 1, 1, V4)) %>%
  mutate(V5 = ifelse(V3 == "transcript" & V7 =="+", as.numeric(V5) + ext, V5)) %>%
  mutate(V4 = ifelse(V3 == "transcript" & V7 =="-", as.numeric(V4) - ext, V4)) %>%
  mutate(V4 = ifelse(V3 == "transcript" & V4 < 1, 1, V4)) %>%
  mutate(V5 = ifelse(V3 == "gene" & V7 =="+", as.numeric(V5) + ext, V5)) %>%
  mutate(V4 = ifelse(V3 == "gene" & V7 =="-", as.numeric(V4) - ext, V4)) %>%
  mutate(V4 = ifelse(V3 == "gene" & V4 < 1, 1, V4))

write.table(gtf.3ext,"geneset.3ext.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#!/usr/bin/env Rscript

require(tidyverse)
require(Rgb)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript process_gtf.R input_file.gtf output_file.gtf", call.=FALSE)
}

input_file = args[1]
output_file = args[2]

gencode_gtf = read.gtf(gzfile(input_file), attr = "split")
gencode_genes = gencode_gtf[gencode_gtf$feature=="gene",c("seqname","start","end","strand","gene_id","gene_name")] %>% distinct() %>% dplyr::rename(chr = seqname)
gencode_genes$gene_length = abs(gencode_genes$end-gencode_genes$start)

write.table(gencode_genes, file = output_file, sep = "\t", quote = F, row.names = F, col.names = T)

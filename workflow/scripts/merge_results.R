# merge together all suez chunks into a single file
# match in chr and pos from VCF
# sort, bgzip and tabix for random access

require(doMC)
require(tidyverse)
require(argparser)
library(gtools)

mixedrank = function(x) order(gtools::mixedorder(x))
registerDoMC(1)
## -------------------------------------------------------------------------------------------------------------------------------------------------

ap = arg_parser(description = "Merge eQTLs results from multiple chunks", name = "merge_results.R")
ap = add_argument(ap, "--vcf_pos", help = "Path to file with positions of each variant in the vcf. (bcftools query -f '%CHROM\t%POS\t%ID\n' file.vcf", type = "character")
ap = add_argument(ap, "--gtf_annot", help = "GTF annotation with gene names and gene positions")
ap = add_argument(ap, "--out_folder", help = "Path to output results file, where chunk folder are stored")
ap = add_argument(ap, "--prefix", help = "Prefix to save the results file name")
ap = add_argument(ap, "--expected_chunks", help = "Number of expected chunks to be merged. Folders should be named output_chunk_NUMBER")

argv = if (interactive()) 
    parse_args(ap, argv = c("--vcf_pos","results/example/vcf_pos.txt.gz",
                            "--gtf_annot","results/example/gencode_genes_annotation.txt",
                            "--out_folder","results/example/",
                            "--prefix","example",
                            "--expected_chunks","40")) else parse_args(ap)
if (interactive())  setwd("~/ad-omics/amp_pd/suezPD/")

vcf_pos <- argv$vcf_pos
gtf_annot <- argv$gtf_annot
prefix <- argv$prefix
out_folder <- argv$out_folder
expected_chunks <- argv$expected_chunks

gtf <- suppressMessages(read_tsv(gtf_annot)) %>% setNames(paste0('gene.', gsub("gene_","",names(.))))
gtf$gene.tss = ifelse(gtf$gene.strand == "+", gtf$gene.start, gtf$gene.end)

# Read variants position to add into final report
snp_df <- read_tsv(gzfile(vcf_pos), col_names = c("chr", "pos", "variant_id"), col_types = c("cnc") )

res_1 = foreach(f = paste0(out_folder, "/output_",1:expected_chunks,"/all_qtl.txt.gz"), .combine = bind_rows, .errorhandling = "remove") %dopar% {
    dat = suppressMessages(read_tsv(f))
}
res_2 = foreach(f = paste0(out_folder, "/output_",1:expected_chunks,"/all_qtl.nom_pval.txt.gz"), .combine = bind_rows, .errorhandling = "remove") %dopar% {
    dat = suppressMessages(read_tsv(f)) %>% pivot_wider(names_from = "test", values_from = c("pval","bonf_pval")) 
}
res_nominal1 <- res_1 %>% left_join(res_2)
res_nominal <- snp_df %>% 
    right_join(res_nominal1, by = c("variant_id"="cis_snp")) %>% 
    left_join(gtf, by = c("gene"="gene.id")) %>% 
    arrange(mixedrank(chr), pos, gene.tss, gene.start, gene.end)

res_geno1 = foreach(f = paste0(out_folder, "/output_",1:expected_chunks,"/all_qtl.geno.pval.txt.gz"), .combine = bind_rows, .errorhandling = "remove") %dopar% {
    dat = suppressMessages(read_tsv(f)) %>% dplyr::select(-nom_pval)
}
res_geno <- snp_df %>% 
    right_join(res_geno1, by = c("variant_id"="cis_snp")) %>% 
    left_join(gtf, by = c("gene"="gene.id")) %>% 
    arrange(mixedrank(chr), pos, gene.tss, gene.start, gene.end)

res_interaction1 = foreach(f = paste0(out_folder, "/output_",1:expected_chunks,"/all_qtl.interact.pval.txt.gz"), .combine = bind_rows, .errorhandling = "remove") %dopar% {
    dat = suppressMessages(read_tsv(f)) %>% dplyr::select(-nom_pval)
}
res_interaction <- snp_df %>% 
    right_join(res_interaction1, by = c("variant_id"="cis_snp")) %>% 
    left_join(gtf, by = c("gene"="gene.id")) %>% 
    arrange(mixedrank(chr), pos, gene.tss, gene.start, gene.end)

write_tsv(res_nominal, file = paste0(out_folder,"/", prefix, ".cis_qtl_nominal.tsv"))
write_tsv(res_geno, file = paste0(out_folder,"/", prefix, ".cis_qtl.tsv"))
write_tsv(res_interaction, file = paste0(out_folder,"/", prefix, ".cis_interaction_qtl.tsv"))

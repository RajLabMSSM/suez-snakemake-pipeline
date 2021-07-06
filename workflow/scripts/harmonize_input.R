# Generate RData with expression and genetic relationship matrix matching
# Remove samples with missing information in the metadata
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(require(tidyverse))
shhh(require(argparser))
shhh(require(S4Vectors))
shhh(require(Rsamtools))

## ----------------------------------------------------------------------------------------------------------------------------------------------------

ap = arg_parser(description = "Harmonize data inputs. Match samples/genes across input files.", name = "harmonize_input.R")
ap = add_argument(ap, "--sample_key", help = " ", default = "sample_key.txt")
ap = add_argument(ap, "--metadata", help = " ", default = "metadata.txt")
ap = add_argument(ap, "--vcf", help = " ", default = "genotypes.vcf.gz")
ap = add_argument(ap, "--grm", help = "Genetic relatedness matrix, calculated by plink. (e.g. plink2 --bfile plink_file --make-grm-bin -out prefix).")
ap = add_argument(ap, "--expression_matrix", help = "Gene expression matrix rdf file (in TPM).", default = "genes_tpm.rds", type = "character")
ap = add_argument(ap, "--tpm_threshold", help = "average TPM threshold.", default = 1)
ap = add_argument(ap, "--gtf", help = "Processed GTF file. Run process_gtf.R first.")
ap = add_argument(ap, "--interaction", default = "case_control", help = "What to look for interactions with.")
ap = add_argument(ap, "--out_folder", help = "Path to output results file, where chunk folder are stored")
ap = add_argument(ap, "--prefix", help = "Prefix to save the results file name")

argv = if (interactive()) 
    parse_args(ap, argv = c("--sample_key","example/sampleKey.txt",
                            "--metadata","example/metadata.txt",
                            "--grm","results/example/example_genotypes",
                            "--expression_matrix","example/genes_tpm.rds",
                            "--gtf","results/example/gencode_genes_annotation.txt",
                            "--interaction","neutrophils_percent_z",
                            "--out_folder","results/example/",
                            "--prefix","example")) else parse_args(ap)

argv = if (interactive()) 
    parse_args(ap, argv = c("--sample_key","/hpc/users/viallr01/ad-omics/amp_pd/suezPD/z_data_mic/migasti_lps_sample_key.tsv",
                            "--metadata","/hpc/users/viallr01/ad-omics/amp_pd/suezPD/z_data_mic/migasti_lps_metadata.tsv",
                            "--grm","results/LPS_run1/LPS_run1_genotypes",
                            "--expression_matrix","/hpc/users/viallr01/ad-omics/amp_pd/suezPD/z_data_mic/genes_tpm.rds",
                            "--gtf","results/LPS_run1/gencode_genes_annotation.txt",
                            "--interaction","LPS",
                            "--out_folder","results/LPS_run1/",
                            "--prefix","LPS_run1.LPS")) else parse_args(ap)

argv = if (interactive()) 
    parse_args(ap, argv = c("--sample_key","/sc/arion/projects/ad-omics/amp_pd/data/2020_v2release_1218/eQTL_mapping/GATK/cell_count_interaction/sampleKey.txt",
                            "--metadata","/sc/arion/projects/ad-omics/amp_pd/data/2020_v2release_1218/processed/all_metadata_021021.txt",
                            "--vcf","/sc/arion/projects/ad-omics/amp_pd/data/2020_v2release_1218/gatk_wgs_qc/output/chrAll_QCFinished_MAF0.01.vcf.gz",
                            "--grm","results/gatk_cell_count_interaction/gatk_cell_count_interaction_genotypes",
                            "--expression_matrix","/sc/arion/projects/ad-omics/amp_pd/data/2020_v2release_1218/processed/salmon_tpm.rds",
                            "--tpm_threshold","1",
                            "--gtf","results/gatk_cell_count_interaction/gencode_genes_annotation.txt",
                            "--interaction","neutrophils_percent_z",
                            "--out_folder","results/gatk_cell_count_interaction/",
                            "--prefix","results/gatk_cell_count_interaction/neutrophils_percent_z")) else parse_args(ap)

#print(argv)

## ----------------------------------------------------------------------------------------------------------------------------------------------------
#rootname=argv$grm
readGRM <- function(rootname){
    grm.file = read.table(paste0(rootname, ".grm"))
    grm.id.file = read.table(paste0(rootname, ".grm.id"))
    
    grm.mat = grm.file %>% pivot_wider(2, names_from = 1, values_from = 4) %>% dplyr::select(-1) %>% as.data.frame()
    grm.mat[lower.tri(grm.mat)] <- t(grm.mat)[lower.tri(grm.mat)]
    rownames(grm.mat) = grm.id.file$V2
    colnames(grm.mat) = grm.id.file$V2    
    
    ret <- list()
    ret$grm <- grm.file
    ret$grm_mat <- grm.mat
    ret$id <- grm.id.file
    return(ret)
}

## 1. Load metadata -----------------------------------------------------------------------------------------------------------------------------------
cat("\n\n### Loading metadata...")
meta = suppressMessages(read_tsv(argv$metadata))

cat(paste0("\n",nrow(meta), " entries."))

## Check interaction term is specified ----------------------------------------------------------------------------------------------------------------
condition = argv$interaction

if (condition %in% colnames(meta)){
    cat(paste0("\nInteraction term '", condition, "' is specified in the metadata provided. Proceeding..."))
    meta2 = meta[,c("sample_id","participant_id",condition)]
    # Remove samples with missing data
    toRemove = is.na(meta2[,condition])
    meta2 = meta2[!toRemove,]
    meta2 = meta2 %>% dplyr::rename("condition" = all_of(condition))
    cat(paste0("\n",sum(toRemove), " samples were removed due to missing information. ", nrow(meta2), " samples remaining." ))
}else{
    stop(paste0("Interaction term '", condition, "' is NOT specified in the metadata provided. Please check you inputs and try again."))
}

## 2. Load sampleKey ---------------------------------------------------------------------------------------------------------------------------------
cat("\n\n### Loading sampleKey...")
sample_key = suppressMessages(read_tsv(argv$sample_key))

cat(paste0("\n",length(unique(sample_key$sample_id)), " unique samples and ", length(unique(sample_key$participant_id)), " unique participants."))

## 3. Load relationship matrix -----------------------------------------------------------------------------------------------------------------------
cat("\n\n### Loading GRM...")
grm_obj <- readGRM(argv$grm)
grm_ids = grm_obj$id
grm = grm_obj$grm_mat
grm[upper.tri(grm)] <- t(grm)[upper.tri(grm)] 
diag(grm) = 1 # is this correct? 

cat(paste0("\n",nrow(grm), " participants included."))

## 4. Load expression --------------------------------------------------------------------------------------------------------------------------------
cat("\n\n### Loading expression data...")
genes_tpm_exp = readRDS(argv$expression_matrix)

cat(paste0("\n",nrow(genes_tpm_exp), " phenotypes (genes) and ", ncol(genes_tpm_exp), " samples."))

## 5. Load GTF ---------------------------------------------------------------------------------------------------------------------------------------
cat("\n\n### Loading GTF...")
gtf  = suppressMessages(read_tsv(argv$gtf)) %>% filter(! chr %in% c("chrX", "chrY", "chrM") ) # only handle autosomes at least for now
gtf = gtf %>% mutate(left = ifelse(strand == "+", start, end), right = left, geneid = gene_id) # just get TSS

cat(paste0("\nYour GTF contains ", nrow(gtf), " phenotypes (genes)."))

## 6. Match samples and donors -----------------------------------------------------------------------------------------------------------------------
cat("\n\n### Matching samples and donors across different datasets...")
meta_samples = unique(meta2$sample_id)
ge_samples = colnames(genes_tpm_exp)

meta_individuals = unique(meta2$participant_id)
genetic_individuals = colnames(grm)

matching_samples = Reduce(base::intersect, list(meta_samples, ge_samples, sample_key$sample_id))
matching_individuals = Reduce(base::intersect, list(meta_individuals, genetic_individuals, sample_key$participant_id))

meta_final = meta2 %>% filter((sample_id %in% matching_samples) & (participant_id %in% matching_individuals)) %>% 
    dplyr::rename("individual" = "participant_id") %>%
    distinct() 
samples = unique(meta_final$sample_id)
individuals = unique(meta_final$individual)

ge = genes_tpm_exp[,samples]

cat(paste0("\n",length(samples), " samples and ", length(individuals), " individuals are matching all input data."))

## 7. Filter GTF by VCF chrom ----------------------------------------------------------------------------------------------------------------------
chr_in_vcf = system(paste0("tabix -l ", argv$vcf), intern = T)
gtf_sub = gtf[ gtf$chr %in% chr_in_vcf, ]
cat(paste0("\n",length(unique(gtf_sub$chr)), " chromosomes intersected the GTF and the VCF provided. "))

## 8. Filter genes by GTF ----------------------------------------------------------------------------------------------------------------------------
genes_to_keep = intersect.Vector( rownames(ge), gtf_sub$gene_id )
ge_sub = ge[genes_to_keep, ]
cat(paste0("\n",length(genes_to_keep), " genes in the expression table intersected the GTF provided."))

## 9. Filter genes by TPM ----------------------------------------------------------------------------------------------------------------------------
cat("\nFiltering genes by median TPM...")
#rs = rowMeans(as.matrix(ge_sub))
rs = apply( as.matrix(ge_sub), MARGIN = 1, FUN = median )

ge_sub1 = ge_sub[rs > argv$tpm_threshold, ] 

ge_norm = ge_sub1 %>% t() 
ge_norm = log10(ge_norm + 0.01) %>% scale() %>% t() %>% as.matrix()

cat("\n", nrow(ge_norm), " phenotypes (genes) with median avg. TPM > ",argv$tpm_threshold," in ", ncol(ge_norm), " samples. Analysis will continue with that number.")

## 10. Write files ------------------------------------------------------------------------------------------------------------------------------------
cat("\n\n### Writing harmonized input files.\n\n")

ge = as.matrix(ge_sub)
grm = grm[meta_final$individual, meta_final$individual]
meta_final = as.data.frame(meta_final)
gtf = as.data.frame(gtf_sub[match(genes_to_keep,gtf_sub$gene_id), ])

save(ge,ge_norm,grm,meta_final,gtf, file = paste0(argv$prefix, "/harmonized.RData"))

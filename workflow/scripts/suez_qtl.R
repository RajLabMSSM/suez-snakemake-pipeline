shhh <- suppressPackageStartupMessages # It's a library, so shhh!

## -------------------------------------------------------------------------------------------------------------------------------------------------
## Data prep for PPMI (Hard coded -- remove in the future)
## -------------------------------------------------------------------------------------------------------------------------------------------------
DATADIR = "~/ad-omics/amp_pd/"
vcfFile = paste0(DATADIR, "WGS-QC-Pipeline/AMPPD_joint/output/chrAll_QCFinished_MAF0.01.vcf.gz")
tabix_file = paste0(DATADIR, "WGS-QC-Pipeline/AMPPD_joint/output/chrAll_QCFinished_MAF0.01.vcf.gz.tbi")

###########################################################################
gtf_file = paste0(DATADIR,"data/2020_v2release_1218/resources/gencode_v29_genes.txt")
# # Reading gencode annotations
# gencode_gtf = read.gtf(gzfile(paste0(DATADIR,"data/2020_v2release_1218/resources/gencode.v29.primary_assembly.annotation.gtf.gz")), attr = "split")
# gencode_genes = gencode_gtf[gencode_gtf$feature=="gene",c("seqname","start","end","strand","gene_id","gene_name")] %>% distinct() %>% dplyr::rename(chr = seqname)
# gencode_genes$gene_length = abs(gencode_genes$end-gencode_genes$start)
# write.table(gencode_genes, file = gtf_file, sep = "\t", quote = F, row.names = F, col.names = T)

###########################################################################
#"plink2 --bfile chrAll_QCFinished_MAF0.01 --pca --make-grm-bin -out chrAll_QCFinished_MAF0.01"
grm_file = paste0(DATADIR,"WGS-QC-Pipeline/AMPPD_joint/output/chrAll_QCFinished_MAF0.01")

###########################################################################
#load(paste0(DATADIR,"/data/2020_v2release_1218/eQTL_mapping/GATK/gene_matrix.RData"))
#saveRDS(genes_tpm, file = paste0(DATADIR,"/suezPD/test/matrix.genes.tpm.Salmon.rds"))
gene_expression_matrix_file = paste0(DATADIR,"suezPD/test/matrix.genes.tpm.Salmon.rds")

###########################################################################
sample_key_file = paste0(DATADIR,"suezPD/test/sample_key.txt")
metadata_file = paste0(DATADIR,"suezPD/test/metadata.txt")

# load(paste0(DATADIR,"/data/2020_v2release_1218/processed/all_metadata_021021.RData"))
# samples_in_rna_and_wgs = read.table(paste0(DATADIR,"/data/2020_v2release_1218/samples_in_rna_and_wgs.list"), header = F, col.names = "participant_id")
# genes_tpm <- readRDS(file = gene_expression_matrix_file)
# 
# clinical_and_tech_metadata$hasWGS = ifelse(clinical_and_tech_metadata$participant_id %in% samples_in_rna_and_wgs$participant_id, T, F)
# clinical_and_tech_metadata = clinical_and_tech_metadata[clinical_and_tech_metadata$sample_id %in% colnames(genes_tpm),]
# 
# # selecting only PPMI samples, and info from disease status, visit, and cell type proportions
# meta <- clinical_and_tech_metadata %>% filter(study == "PPMI", hasWGS==T) %>%
#     select(sample_id, participant_id, case_control_other_latest, visit_month,
#            has_known_gba_mutation_in_wgs,has_known_lrrk2_mutation_in_wgs,has_known_snca_mutation_in_wgs,has_known_pd_mutation_in_wgs,
#            neutrophils_percent,lymphocytes_percent,monocytes_percent,eosinophils_percent,basophils_percent) %>% distinct()
# # create column status with CASE and CONTROL only (OTHER change to NA to be removed)
# meta$status = meta$case_control_other_latest
# meta$status[!meta$case_control_other_latest %in% c("Case","Control")] <- NA
# # create progression column (status + visit)
# meta$progression = paste0(meta$visit_month,"_",meta$case_control_other_latest)
# # create standardized values for cell proportions
# meta$monocytes_percent_z = scale(as.numeric(meta$monocytes_percent))
# meta$lymphocytes_percent_z = scale(as.numeric(meta$lymphocytes_percent))
# meta$basophils_percent_z = scale(as.numeric(meta$basophils_percent))
# meta$eosinophils_percent_z = scale(as.numeric(meta$eosinophils_percent))
# meta$neutrophils_percent_z = scale(as.numeric(meta$neutrophils_percent))
# # create numeric mutation columns (1 = yes, 0 = no)
# meta$has_known_gba_mutation_in_wgs = as.numeric(as.factor(meta$has_known_gba_mutation_in_wgs))-1
# meta$has_known_lrrk2_mutation_in_wgs = as.numeric(as.factor(meta$has_known_lrrk2_mutation_in_wgs))-1
# meta$has_known_pd_mutation_in_wgs = as.numeric(as.factor(meta$has_known_pd_mutation_in_wgs))-1
# meta$has_known_snca_mutation_in_wgs = as.numeric(as.factor(meta$has_known_snca_mutation_in_wgs))-1
# 
# sample_key = meta %>% select(sample_id, participant_id) %>% distinct()
# 
# write_tsv(sample_key, sample_key_file)
# write_tsv(meta, metadata_file)

# # Toy example (30 donors, 150 samples)
# sampleKey = read_tsv(sample_key_file)
# metadata = read_tsv(metadata_file)
# genes_tpm = readRDS(gene_expression_matrix_file)
# 
# donors = sampleKey %>% group_by(participant_id) %>% summarise(n = n()) %>% arrange(-n) %>% head(n = 30)
# samplesKey_toy = sampleKey[sampleKey$participant_id %in% donors$participant_id,]
# write_tsv(samplesKey_toy, file = paste0(DATADIR,"suezPD/example/sampleKey.txt"))
# 
# rownames(metadata) = metadata$sample_id
# metadata_toy = metadata[samplesKey_toy$sample_id,]
# write_tsv(metadata_toy, file = paste0(DATADIR,"suezPD/example/metadata.txt"))
# 
# rs = rowMeans(as.matrix(genes_tpm))
# genes_tpm_toy = genes_tpm[rs > 10, samplesKey_toy$sample_id] # TPM < 0
# saveRDS(genes_tpm_toy, file = paste0(DATADIR,"suezPD/example/genes_tpm.rds"))

#print(R.version)
#cat(".libPaths:", .libPaths())

shhh(library(argparser))

## -------------------------------------------------------------------------------------------------------------------------------------------------

ap = arg_parser(description = "map (interaction) eQTLs with suez", name = "suez_qtl.R")
ap = add_argument(ap, "which_genes", default = 1, help = "seq(which_genes, num_genes, by=num_jobs) are the gene indices that will be processed (for splitting across compute cluster).")
ap = add_argument(ap, "interaction", default = "case_control", help = "What to look for interactions with.")
ap = add_argument(ap, "--expression_matrix", help = "Gene expression matrix rdf file (in TPM).", default = "genes_tpm.rds", type = "character")
ap = add_argument(ap, "--sample_key", help = " ", default = "sample_key.txt")
ap = add_argument(ap, "--metadata", help = " ", default = "metadata.txt")
ap = add_argument(ap, "--vcf", help = "vcf file.", default = "genotypes.vcf.gz")
ap = add_argument(ap, "--tabix", help = "tabix for vcf file. Default vcf file name followed by .tbi")
ap = add_argument(ap, "--grm", help = "Genetic relatedness matrix, calculated by plink. (e.g. plink2 --bfile plink_file --make-grm-bin -out prefix).")
ap = add_argument(ap, "--gtf", help = "Processed GTF file. Run process_gtf.R first.")
ap = add_argument(ap, "--num_jobs", help = "number of jobs being run on the cluster.", default = 40)
ap = add_argument(ap, "--threads", help = "number of cores to use.", default = 20)
ap = add_argument(ap, "--cisdist", help = "max distance from TSS to look for SNPs.", default = "100000")
ap = add_argument(ap, "--checkpoint_dir", help = "where to save results. Defaults to checkpoint_{interaction}_{cisdist}_{maf_threshold}_{normalization_approach}_{permutation_approach}",  flag = F)
ap = add_argument(ap, "--out_dir", help = " ", default = "output")
ap = add_argument(ap, "--maf_threshold", help = "MAF threshold.", default = 0.05)
ap = add_argument(ap, "--normalization_approach", help = "How to normalize GE, one of: qq, linear, log.", default = "qq")
ap = add_argument(ap, "--permutation_approach", help = "One of (none,permute,boot). We recommend boot (parametric bootstrap) since the permutation test is not really valid for interaction effect testing.", default = "boot")
ap = add_argument(ap, "--ge_PCs_to_remove", help = "Number of gene expression PCs to remove.", default = 20)
ap = add_argument(ap, "--build", help = "Genome build used for the genotype data.", default = "GRCh38")

argv = if (interactive()) 
    parse_args(ap, argv = c("1", 
                            "visit_month",
                            #"visit_month",
                            #"monocytes_percent_z",
                            "--expression_matrix",gene_expression_matrix_file,
                            "--metadata",metadata_file,
                            "--sample_key",sample_key_file,
                            "--cisdist","10000",
                            "--vcf",vcfFile,
                            "--tabix",tabix_file,
                            "--grm",grm_file,
                            "--gtf",gtf_file,
                            "--out_dir",DATADIR,
                            "--checkpoint_dir",paste0(DATADIR,"suezPD/"))) else parse_args(ap)

argv = if (interactive()) 
    parse_args(ap, argv = c("3", 
                            "neutrophils_percent_z",
                            "--expression_matrix","example/genes_tpm.rds",
                            "--metadata","example/metadata.txt",
                            "--vcf","example/genotypes.vcf.gz",
                            "--sample_key","example/sampleKey.txt",
                            "--cisdist","100000",
                            "--gtf","results/example/gencode_genes_annotation.txt",
                            "--grm","results/example/example_genotypes",
                            "--out_dir","results/example/neutrophils_percent_z/output_3")) else parse_args(ap)

if (is.na(argv$tabix))
    argv$tabix = paste0(argv$vcf,".tbi")

if (is.na(argv$checkpoint_dir))
    argv$checkpoint_dir = paste0(argv$out_dir,paste0(paste("/suez_checkpoints", 
                                                              argv$interaction, 
                                                              argv$cisdist, 
                                                              argv$maf_threshold, 
                                                              argv$normalization_approach, 
                                                              argv$permutation_approach, 
                                                              sep = "_"), "/")) else 
    argv$checkpoint_dir = paste0(argv$out_dir,"/",argv$checkpoint_dir)

#print(argv)

#if (interactive()) .libPaths(c("/usr/local/lib64/R/library","/sc/arion/projects/ad-omics/ricardo/Rlib4"))
#if (interactive()) setwd("/sc/arion/projects/ad-omics/amp_pd/suezPD/")
## -------------------------------------------------------------------------------------------------------------------------------------------------
shhh(library(plinkFile))
suppressWarnings(shhh(require(suez)))
shhh(require(tidyverse))
shhh(require(doMC))
shhh(require(magrittr))
shhh(require(VariantAnnotation))
shhh(require(irlba))
shhh(library(Rgb))

select = dplyr::select
rename = dplyr::rename

## -------------------------------------------------------------------------------------------------------------------------------------------------
#source("/sc/arion/projects/ad-omics/ricardo/MyRepo/suez/R/map_interaction_qtl.R")
#source("/sc/arion/projects/ad-omics/ricardo/MyRepo/suez/R/stanmodels.R")
## -------------------------------------------------------------------------------------------------------------------------------------------------
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
## -------------------------------------------------------------------------------------------------------------------------------------------------
registerDoMC(argv$threads)

## VCF file ----------------------------------------------------------------------------------------------------------------------------------------
tab = TabixFile(argv$vcf, argv$tabix)
## Load relationship matrix ------------------------------------------------------------------------------------------------------------------------
print("### Loading GRM...")
grm_obj <- readGRM(argv$grm)
grm_ids = grm_obj$id
grm = grm_obj$grm_mat
grm[upper.tri(grm)] <- t(grm)[upper.tri(grm)] 
diag(grm) = 1 # is this correct? 

print(paste0(nrow(grm), " participants included."))

## Load GTF ----------------------------------------------------------------------------------------------------------------------------------------
print("### Loading GTF...")
gtf  = suppressMessages(read_tsv(argv$gtf)) %>% filter(! chr %in% c("chrX", "chrY") ) # only handle autosomes at least for now
gtf = gtf %>% mutate(left = ifelse(strand == "+", start, end), right = left, geneid = gene_id) # just get TSS

print(paste0("Your GTF contains ", nrow(gtf), " phenotypes (genes)."))

## Load expression ---------------------------------------------------------------------------------------------------------------------------------
print("### Loading expression data...")
genes_tpm_exp = readRDS(argv$expression_matrix)

print(paste0(nrow(genes_tpm_exp), " phenotypes (genes) and ", ncol(genes_tpm_exp), " samples."))

## Load metadata -----------------------------------------------------------------------------------------------------------------------------------
print("### Loading metadata...")
meta = suppressMessages(read_tsv(argv$metadata))

print(paste0(nrow(meta), " entries."))

## Load sampleKey -----------------------------------------------------------------------------------------------------------------------------------
print("### Loading sampleKey...")
sample_key = suppressMessages(read_tsv(argv$sample_key))

print(paste0(length(unique(sample_key$sample_id)), " unique samples and ", length(unique(sample_key$participant_id)), " unique participants."))

## Build interaction term --------------------------------------------------------------------------------------------------------------------------
condition = argv$interaction
meta2 = meta[,c("sample_id","participant_id",condition)]
# Remove samples with missing data
meta2 = meta2[!is.na(meta2[,condition]),]
meta2 = meta2 %>% dplyr::rename("condition" = all_of(condition))

## Match samples and donors ------------------------------------------------------------------------------------------------------------------------
print("### Matching samples and donors across different datasets...")
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

print(paste0(length(samples), " samples and ", length(individuals), " individuals matching."))

## Filter genes by TPM and log transform -----------------------------------------------------------------------------------------------------------
rs = rowMeans(as.matrix(ge))
ge = ge[rs > 0, ] # TPM < 0

ge_norm = ge %>% t() 
ge_norm = log10(ge_norm + 0.01) %>% scale() %>% t()

## Regress covariates ------------------------------------------------------------------------------------------------------------------------------
X = model.matrix( ~ condition, data = meta_final)

b = solve( t(X) %*% X, t(X) %*% t(ge_norm) )
resi = as.matrix(ge_norm - t(X %*% b))

svd_ge = irlba(resi, argv$ge_PCs_to_remove)

recons = svd_ge$u %*% diag(svd_ge$d) %*% t(svd_ge$v)

ge_pc = (ge_norm - recons) # %>% t() %>% scale() %>% t() # mapping rescales anyway

## Filter genes by GTF -----------------------------------------------------------------------------------------------------------------------------
genes_to_keep = intersect.Vector( rownames(ge_pc), gtf$gene_id )
ge_pc_sub = ge_pc[genes_to_keep, ]

## run suez ----------------------------------------------------------------------------------------------------------------------------------------
nr = if (interactive()) 1 else nrow(ge_pc_sub)

genes_ind_here = seq(argv$which_genes, nrow(ge_pc_sub), argv$num_jobs)
#genes_ind_here = seq(1, nr, argv$num_jobs)
#genes_ind_here = (grep("ENSG00000142731.10",rownames(ge_pc_sub)))

all_eqtl = suppressMessages(map_interaction_qtl(input = as.matrix(ge_pc_sub[genes_ind_here,,drop=F]),
                                                genotype = tab,
                                                genome_build = argv$build,
                                                geneloc = as.data.frame(gtf),
                                                snploc = NULL,
                                                anno = meta_final,
                                                debug = F, # interactive(),
                                                normalization_approach = argv$normalization_approach,
                                                permutation_approach = argv$permutation_approach,
                                                cisdist = as.numeric(argv$cisdist), # 100000,
                                                sample_kernel = grm[meta_final$individual, meta_final$individual],
                                                checkpoint_dir = argv$checkpoint_dir,
                                                maf_threshold = argv$maf_threshold))

#print(head(all_eqtl))
all_eqtl %>% write_tsv(paste0(argv$out_dir,"/all_qtl.txt.gz"))

# # Outputs
# l0: loglikelihood for model with just the condition (e.g. case vs control)
# l_geno: loglikelihood with genotype term added
# l_interact: loglikelihood with interaction term added
# l_boot_geno and l_boot_interact are for parametric bootstrap

## P-values ----------------------------------------------------------------------------------------------------------------------------------------
## enable running this code from the cached results
P = length(unique(meta_final$condition)) # 15

# Get P-values for all pairs
res_nom = foreach(f = list.files(argv$checkpoint_dir), .combine = bind_rows, .errorhandling = "remove") %dopar% {
    dat = suppressMessages(read_tsv(paste0(argv$checkpoint_dir, f)))
    dat = dat %>% mutate( geno = pchisq( 2.0*(l_geno - l0), lower.tail = F , 1 ),
                          interact = pchisq( 2.0*(l_interact - l_geno), lower.tail = F , P-1 ))
    dat %>% select(gene, cis_snp, geno, interact) %>% 
        pivot_longer(c(geno, interact), names_to = "test", values_to = "pval") %>% 
        group_by(gene,test) %>% 
        mutate(bonf_pval = p.adjust(pval, n = n(), method = "bonferroni")) 
}
res_nom %>% write_tsv(paste0(argv$out_dir, "/all_qtl.nom_pval.txt.gz"))

# Get P-values for top associations
res = foreach(f = list.files(argv$checkpoint_dir), .combine = bind_rows, .errorhandling = "remove") %dopar% {
    dat = suppressMessages(read_tsv(paste0(argv$checkpoint_dir, f)))
    dat = dat %>% mutate( geno = pchisq( 2.0*(l_geno - l0), lower.tail = F , 1 ),
                          interact = pchisq( 2.0*(l_interact - l_geno), lower.tail = F , P-1 ))
    dat %>% select(gene, cis_snp, geno, interact) %>% 
        pivot_longer(c(geno, interact), names_to = "test", values_to = "pval") %>% 
        group_by(gene,test) %>% 
        mutate(bonf_pval = p.adjust(pval, n = n(), method = "bonferroni")) %>% 
        top_n(1, -pval)
}
res %>% write_tsv(paste0(argv$out_dir, "/all_qtl.pval.txt.gz"))

## regular eQTLs
res_geno = res %>% filter(test == "geno") %>% ungroup() %>% select(-test) %>% group_by(gene) %>% dplyr::slice(1) %>% ungroup() %>% 
    mutate(nom_pval = pval, bonf_pval = pmin(bonf_pval, 1), q = p.adjust(bonf_pval, method="BH"))
cat("eQTLs at 5% FDR:", sum(res_geno$q < 0.05), "\n")

res_geno %>% write_tsv(paste0(argv$out_dir, "/all_qtl.geno.pval.txt.gz"))

## interaction eqtl
res_interact = res %>% filter(test == "interact") %>% ungroup() %>% select(-test) %>% group_by(gene) %>% dplyr::slice(1) %>% ungroup() %>% 
    mutate(nom_pval = pval, bonf_pval = pmin(bonf_pval, 1), q = p.adjust(bonf_pval, method="BH"))
cat("eQTLs at 5% FDR:", sum(res_interact$q < 0.05), "\n")

res_interact %>% write_tsv(paste0(argv$out_dir, "/all_qtl.interact.pval.txt.gz"))

shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(library(argparser))

## -------------------------------------------------------------------------------------------------------------------------------------------------
ap = arg_parser(description = "map (interaction) eQTLs with suez", name = "suez_qtl.R")
ap = add_argument(ap, "which_genes", default = 1, help = "seq(which_genes, num_genes, by=num_jobs) are the gene indices that will be processed (for splitting across compute cluster).")
ap = add_argument(ap, "interaction", default = "case_control", help = "What to look for interactions with.")
ap = add_argument(ap, "--harmonizedRData", help = "Harmonized RData.", default = "harmonized.RData", type = "character")
ap = add_argument(ap, "--vcf", help = "vcf file.", default = "genotypes.vcf.gz")
ap = add_argument(ap, "--tabix", help = "tabix for vcf file. Default vcf file name followed by .tbi")
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
    parse_args(ap, argv = c("87", 
                            "case_control_other_latest",
                            "--harmonizedRData","results/example/example.case_control_other_latest.harmonized.RData",
                            "--vcf","example/genotypes.vcf.gz",
                            "--cisdist","10000",
                            "--num_jobs","200",
                            "--threads","1",
                            "--out_dir","results/example/case_control_other_latest/output_87")) else parse_args(ap)

argv = if (interactive()) 
    parse_args(ap, argv = c("1", 
                            "LPS",
                            "--harmonizedRData","results/LPS_run1/LPS/harmonized.RData",
                            "--vcf","z_data_mic/gsa_microglia_05-2021_allAncestry_MAF.vcf.gz",
                            "--cisdist","10000",
                            "--num_jobs","200",
                            "--threads","1",
                            "--out_dir","results/LPS_run1/LPS/output_1")) else parse_args(ap)

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

#if (interactive()) .libPaths(c("/usr/local/lib64/R/library","/sc/arion/projects/ad-omics/ricardo/Rlib4"))
if (interactive()) setwd("/sc/arion/projects/ad-omics/amp_pd/suezPD/")
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
registerDoMC(argv$threads)

## VCF file ----------------------------------------------------------------------------------------------------------------------------------------
tab = TabixFile(argv$vcf, argv$tabix)

## Using Harmonized data  --------------------------------------------------------------------------------------------------------------------------
load(argv$harmonizedRData) # expression_sub,grm_sub,meta_sub,gtf_sub2

## Filter genes by TPM and log transform -----------------------------------------------------------------------------------------------------------
samples = unique(meta_final$sample_id)
individuals = unique(meta_final$individual)
print(paste0(length(samples), " samples and ", length(individuals), " individuals matching."))

print("### Filtering genes by TPM...")
rs = rowMeans(as.matrix(ge))
ge = ge[rs > 5, ] # TPM < 0

ge_norm = ge %>% t() 
ge_norm = log10(ge_norm + 0.01) %>% scale() %>% t()

print(paste0(nrow(ge_norm), " phenotypes (genes) with avg. TPM > 5 in ", ncol(ge_norm), " samples."))

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
genes_ind_here = seq(argv$which_genes, nrow(ge_pc_sub), argv$num_jobs)

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

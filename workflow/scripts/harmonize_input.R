# Generate RData with expression and genetic relationship matrix matching
# Remove samples with missing information in the metadata
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(require(tidyverse))
shhh(require(argparser))
shhh(require(S4Vectors))

## ----------------------------------------------------------------------------------------------------------------------------------------------------

ap = arg_parser(description = "Harmonize data inputs. Match samples/genes across input files.", name = "harmonize_input.R")
ap = add_argument(ap, "--sample_key", help = " ", default = "sample_key.txt")
ap = add_argument(ap, "--metadata", help = " ", default = "metadata.txt")
ap = add_argument(ap, "--grm", help = "Genetic relatedness matrix, calculated by plink. (e.g. plink2 --bfile plink_file --make-grm-bin -out prefix).")
ap = add_argument(ap, "--expression_matrix", help = "Gene expression matrix rdf file (in TPM).", default = "genes_tpm.rds", type = "character")
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

print(argv)

## ----------------------------------------------------------------------------------------------------------------------------------------------------
readGRM <- function(rootname){
    bin.file.name <- paste(rootname, ".grm.bin", sep="")
    n.file.name <- paste(rootname, ".grm.N.bin", sep="")
    id.file.name <- paste(rootname, ".grm.id", sep="")
    
    cat("Reading IDs\n")
    id <- read.table(id.file.name)
    n <- dim(id)[1]
    cat("Reading GRM\n")
    bin.file <- file(bin.file.name, "rb")
    grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
    close(bin.file)
    cat("Reading N\n")
    n.file <- file(n.file.name, "rb")
    N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
    close(n.file)
    
    cat("Creating data frame")
    l <- list()
    for(i in 1:n)
    {
        l[[i]] <- 1:i
    }
    col1 <- rep(1:n, 1:n)
    col2 <- unlist(l)
    grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)
    
    grm_mat <- matrix(NA, nrow = n, ncol = n)
    for(i in 1:nrow(grm)){
        grm_mat[grm$id1[i],grm$id2[i]] <- grm$grm[i]
    }
    dimnames(grm_mat) = list(id$V2, id$V2)
    
    ret <- list()
    ret$grm <- grm
    ret$grm_mat <- grm_mat
    ret$id <- id
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
gtf  = suppressMessages(read_tsv(argv$gtf)) %>% filter(! chr %in% c("chrX", "chrY") ) # only handle autosomes at least for now
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

## 7. Filter genes by GTF ----------------------------------------------------------------------------------------------------------------------------
genes_to_keep = intersect.Vector( rownames(ge), gtf$gene_id )
ge_sub = ge[genes_to_keep, ]
cat(paste0("\n",length(genes_to_keep), " genes in the expression table intersected the GTF provided. Analysis will continue with that number."))

## 8. Write files ------------------------------------------------------------------------------------------------------------------------------------
cat("\n\n### Writing harmonized input files.\n\n")

expression_sub = as.matrix(ge_sub)
grm_sub = grm[meta_final$individual, meta_final$individual]
meta_sub = as.data.frame(meta_final)
gtf_sub = as.data.frame(gtf[match(genes_to_keep,gtf$gene_id), ])

save(expression_sub,grm_sub,meta_sub,gtf_sub, file = paste0(argv$out_folder,"/", argv$prefix, ".harmonized.RData"))

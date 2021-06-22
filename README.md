# *suez-snakemake-pipeline*: Snakemake pipeline for calling QTLs using suez

suez extends the PANAMA eQTL mapping framework (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002330) to handle

multiple conditions (e.g. stimulated/unstimulated), including response-eQTL mapping
latent factor correction
known kinship between individuals

### Install requirements:

First install suez R package.
```
# Using R (ml R/4.0.3)
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("ricardovialle/suez")
```

Then install snakemake. Using a separated environment is recommended.
```
conda create -n snakemake_env -c bioconda snakemake
conda activate snakemake_env # Command to activate the environment. To deactivate use "conda deactivate"
```

### Inputs:

Input files must be described in a config file (`config.yaml`). Description of each file is described bellow.
```
# ID to identify the run folder and result files
dataCode: example 

# Tab separated file with two columns named: sample_id,	participant_id
# Each row is a sample (participant_id can be repeated, in case of multiple samples for the same donor)
sampleKey: "example/sampleKey.txt"

# R data.frame matrix with gene expression in TPM for each sample (indicated in the sampleKey/sample_id column)
countMatrixRData: "example/genes_tpm.rds"

# Table with metadata information for each sample. Must have two columns named: sample_id,	participant_id. IDs will be used for matching with the other tables.
metadataRData: "example/metadata.txt"

# Flag to run interaction eQTL (for now it needs to be True)
interaction: True
# Column names in metadata to be used as interaction term (samples with NA will be droped)
interaction_name:
    - 'case_control_other_latest'
    #- 'neutrophils_percent_z'
    #- 'visit_month'

# VCF file with IDs corresponding to sampleKey/participant_id column. Incompatible IDs will not be included.
VCF: "example/genotypes.vcf.gz"

# Gene annotation file for collect gene positions (use the same version as the one used in the RNAseq mapping)
GTF: "ref/gencode.v29.annotation.gtf.gz"

# For now only eQTL. Maybe in the future we can add other models. 
mode: "eQTL"
```

### Run:
```
snakemake --configfile config/config.yaml -pr --cores 20 --use-conda
```

### Outputs:

Results will be saved in a folder named `results/{dataCode}/{interaction_name}`.

eQTL resutls files will be generated in the end.

1. nominal results (`results/{dataCode}/{interaction_name}/{dataCode}.cis_qtl_nominal.tsv`):

**chr**|**pos**|**variant\_id**|**gene**|**l0**|**l\_geno**|**l\_interact**|**l\_boot\_geno**|**l\_boot\_interact**|**pval\_geno**|**pval\_interact**|**bonf\_pval\_geno**|**bonf\_pval\_interact**|**gene.chr**|**gene.start**|**gene.end**|**gene.strand**|**gene.name**|**gene.length**|**gene.tss**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
chr1|1275912|rs6685064|ENSG00000175756.13|-34.96882145519742|-34.757420777367074|-34.748307195058686|-21.516580738109866|-21.149373257868078|0.5155425218133911|0.990927820511392|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1277269|rs59225223|ENSG00000175756.13|-34.96882145519742|-34.757420777367074|-34.748307195058686|-20.193890580695758|-20.158508131786952|0.5155425218133911|0.990927820511392|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1279694|rs74046664|ENSG00000175756.13|-34.96882145519742|-34.757420777367074|-34.748307195058686|-25.72567495847092|-25.72538249822327|0.5155425218133911|0.990927820511392|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1279698|rs74046665|ENSG00000175756.13|-34.96882145519742|-34.757420777367074|-34.748307195058686|-9.33206894248778|-8.748470338926012|0.5155425218133911|0.990927820511392|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1279728|rs79179971|ENSG00000175756.13|-34.96882145519742|-34.757420777367074|-34.748307195058686|-32.95859554278356|-31.895592011565164|0.5155425218133911|0.990927820511392|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1280044|rs1268339|ENSG00000175756.13|-34.96882145519742|-33.08131181802449|-32.198863839941126|-41.673529359597104|-41.184382735676806|0.05202330152608478|0.4137687740042784|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1282558|rs77709654|ENSG00000175756.13|-34.96882145519742|-34.919917384231255|-33.807176369130794|-42.34123549885604|-41.79207226292382|0.7544759749583494|0.3286568718712499|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1282706|rs6603788|ENSG00000175756.13|-34.96882145519742|-34.757420777367074|-34.748307195058686|-41.54259001123501|-41.22000429579967|0.5155425218133911|0.990927820511392|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|1283275|rs114027205|ENSG00000175756.13|-34.96882145519742|-33.15137927271206|-32.813045209343535|-34.61485239937534|-33.56196150369829|0.05658053488253116|0.7129570752692607|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495


2. top eQTL per phenotype (`results/{dataCode}/{interaction_name}/{dataCode}.cis_qtl.tsv`):

**chr**|**pos**|**variant\_id**|**gene**|**pval**|**bonf\_pval**|**q**|**gene.chr**|**gene.start**|**gene.end**|**gene.strand**|**gene.name**|**gene.length**|**gene.tss**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
chr1|1366511|rs144090606|ENSG00000175756.13|0.0029108147048559916|0.8848876702762214|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|2286690|rs59476564|ENSG00000157933.9|0.008409021938056092|1|1|chr1|2228695|2310119|+|SKI|81424|2228695
chr1|2384241|rs2843128|ENSG00000157916.19|0.0013997496829830747|0.6816780956127574|1|chr1|2391775|2405444|+|RER1|13669|2391775
chr1|6757900|rs10864637|ENSG00000007923.15|0.11968310873907757|1|1|chr1|6634168|6701924|-|DNAJC11|67756|6701924
chr1|7922756|rs17229081|ENSG00000116288.12|0.0067034587057322875|1|1|chr1|7954291|7985505|+|PARK7|31214|7954291
chr1|8877167|rs12752143|ENSG00000074800.15|0.006357460028843177|1|1|chr1|8861000|8879250|-|ENO1|18250|8879250
chr1|10382809|rs11121556|ENSG00000142657.20|8.164516809670452e-5|0.020901163032756357|0.08987500104085233|chr1|10398592|10420144|+|PGD|21552|10398592
chr1|11693901|rs2038131|ENSG00000177674.15|8.878139188241366e-4|0.47498044657091304|1|chr1|11736084|11754802|+|AGTRAP|18718|11736084
chr1|12068441|rs12029016|ENSG00000028137.18|0.0029251314238791922|0.7722346959041068|1|chr1|12167003|12209228|+|TNFRSF1B|42225|12167003


3. interaction results (`results/{dataCode}/{interaction_name}/{dataCode}.cis_interaction_qtl.tsv`):

**chr**|**pos**|**variant\_id**|**gene**|**pval**|**bonf\_pval**|**q**|**gene.chr**|**gene.start**|**gene.end**|**gene.strand**|**gene.name**|**gene.length**|**gene.tss**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
chr1|1360502|rs76580958|ENSG00000175756.13|0.059790518605557605|1|1|chr1|1373730|1375495|-|AURKAIP1|1765|1375495
chr1|2285267|rs2971436|ENSG00000157933.9|0.09787023522453538|1|1|chr1|2228695|2310119|+|SKI|81424|2228695
chr1|2300156|rs1809822|ENSG00000157916.19|0.014506799200068858|1|1|chr1|2391775|2405444|+|RER1|13669|2391775
chr1|6736962|rs75272952|ENSG00000007923.15|0.1937720962540602|1|1|chr1|6634168|6701924|-|DNAJC11|67756|6701924
chr1|7911018|rs679563|ENSG00000116288.12|0.004938830511635669|1|1|chr1|7954291|7985505|+|PARK7|31214|7954291
chr1|8922405|rs6690365|ENSG00000074800.15|0.0250960208424478|1|1|chr1|8861000|8879250|-|ENO1|18250|8879250
chr1|10465784|rs34132699|ENSG00000142657.20|0.008133530739828056|1|1|chr1|10398592|10420144|+|PGD|21552|10398592
chr1|11764246|rs4845877|ENSG00000177674.15|0.01707645517419329|1|1|chr1|11736084|11754802|+|AGTRAP|18718|11736084
chr1|12151199|rs12409061|ENSG00000048707.15|0.23938244621519544|1|1|chr1|12230030|12512047|+|VPS13D|282017|12230030

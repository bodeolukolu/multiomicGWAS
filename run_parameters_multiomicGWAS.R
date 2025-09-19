#!/usr/bin/env Rscript

source("https://github.com/bodeolukolu/multiomicGWAS/raw/refs/heads/main/multiomicGWAS.R")

multiomicGWAS <- (
    wdir = "./",
    projname = "GWAS",
    ploidy_levels = c("2","4","6","8"),
    trait_names = c("trait1","trait2"),
    model_effect = c("Add","Dom"),
    fdr = TRUE,
    bonferroni = TRUE,
    suggestive = "5",
    perm = "1",
    cores = "1",
    genofile_2x = "SNP.txt",
    genofile_4x = "SNP.txt",
    genofile_6x = "SNP.txt",
    genofile_8x = "SNP.txt",
    phenofile = "traits.txt",
    method = "MLM",
    covariate_pheno = "trait",
    covariate_metag = FALSE,
    maf = "0.02",
    LOCO = FALSE,
    pheno_taxa_strain = "taxa1",
    pheno_taxa_species = "taxa2",
    metag_data_strains = "metag.txt",
    metag_data_species = "metag.txt"
)

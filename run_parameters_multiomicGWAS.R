#!/usr/bin/env Rscript

# # List of required packages
# pkgs <- c("AGHmatrix", "qqplotr", "ggplot2", "dplyr", "data.table", "stringr",
#           "heatmaply", "ppcor", "zoo", "GGally", "reshape2", "compositions", "sommer", "mice", "qvalue")
# # Install missing packages
# for (p in pkgs) {
#   if (!requireNamespace(p, quietly = TRUE)) {
#     message("Installing package: ", p)
#     install.packages(p, repos = "https://cloud.r-project.org")
#   }
# }
# install.packages("remotes")
# library("remotes")
# remotes::install_github("jendelman/GWASpoly")


source("https://github.com/bodeolukolu/multiomicGWAS/raw/refs/heads/main/multiomicGWAS.R")

multiomicGWAS (
    wdir = "./",
    projname = "GWAS",
    ploidy_levels = c(2,4,6,8),
    trait_names = c("trait1","trait2"),
    model_effect = c("Add","Dom"),
    fdr = TRUE,
    bonferroni = TRUE,
    suggestive = "5",
    perm = "1",
    cores = "1",
    genofile_2x = NULL,
    genofile_4x = NULL,
    genofile_6x = NULL,
    genofile_8x = NULL,
    phenofile = "traits.txt",
    method = c("MLM", "GLM"),
    covariate_pheno = "trait",
    covariate_metag = FALSE,
    maf = "0.02",
    LOCO = FALSE,
    pheno_taxa_strain = "taxa1",
    pheno_taxa_species = "taxa2",
    metag_data_strains = "metag.txt",
    metag_data_species = "metag.txt"
)

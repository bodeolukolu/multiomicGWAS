
# Introduction
multiomicGWAS is a flexible GWAS pipeline leveraging GWASpoly to perform genome-wide association analyses across ploidy levels from 2 to 8. The pipeline can integrate secondary traits and host-associated metagenome/microbiome data as covariate, and it can also allow each microbial taxon to be treated as an independent phenotypic trait for GWAS. Relationship matrices are computed automatically: genomic relationship matrices (GRM) using AGHmatrix and microbiome/metagenome kernels using the Aitchison method for compositional data, enabling streamlined multi-omic association analyses.

For questions, bugs, and suggestions, please contact bolukolu@utk.edu.

## How to run multiomicGWAS
Download the <run_parameters_multiomicGWAS.R>, edit parameters and run the pipeline.

## Description of parameters


## Dependencies
- GWASpoly
- AGHmatrix
- qqplotr
- ggplot2
- dplyr
- data.table
- stringr
- heatmaply
- ppcor
- zoo
- GGally
- reshape2
- compositions
- sommer
- mice
- qvalue


## Select Article Referencing GBSapp


## Acknowledgment


## Troubleshooting

## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License
<a href="https://github.com/bodeolukolu/GBSapp/blob/master/LICENSE">Apache License Version 2.0</a>

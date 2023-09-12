# snpboost
An algorithm to apply statistical boosting on genotype data via a batch-based approach as described in https://www.biorxiv.org/content/10.1101/2022.04.29.489836v1.

snpboost includes R functions in which the boosting algorithm is implemented.

Simulations provides R Code and bash scripts that were used to run the simulation on phenotypes based on UKBB genotype data and fit statistical boosting models on those.

UK_biobank_application provides R scripts that were used to fit snpboost and snpnet models on UKBB data.

Installation:
The following requirements of snpboost are available from CRAN:
  - tidyverse
  - data.table
  - Rfast
  - parallel
  - dataPreparation
  - glmnet
    
Like snpnet, it also depends on the pgenlibr package. One can install it by running the following commands in R. Notice that the installation of pgenlibr requires zstd(>=1.4.4). It can be built from source or simply available from conda, pip or brew.

```r
library(devtools)
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```

Additionally PLINK 2.0 is required. It can be installed from https://www.cog-genomics.org/plink/2.0/.

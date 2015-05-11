<!-- README.md is generated from README.Rmd. Please edit that file -->




argyle
======

An `R` package for import, QC and (some) analysis of genotyping and hybridization-intensity data from Illumina Infinium arrays.

Usage
-----

``` {.r}
library(argyle)

data(snps)
geno <- read.beadstudio("sample", snps, in.path = "./")
```

Interface to `PLINK`
--------------------

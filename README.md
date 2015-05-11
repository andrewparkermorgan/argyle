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

## see marker map and sample metadata
markers(geno)
samples(geno)

## subset operations: equivalent pairs
subset(geno, chr == "chrX")
geno[ markers(geno)$chr == "chrX", ]

subset(geno, sex == 2, by = "samples")
geno[ ,samples(geno)$sex == 2 ]

## run QC checks and flag samples above thresholds
geno <- run.qc.checks(geno, max.H = 3000, max.N = 5000)
# how many samples fail QC?
summarize.filters(geno)
# remove samples failing QC
filt <- apply.filters(geno)

## compute BAF/LRR for copy-number analysis
data(clusters)
geno.norm <- tQN(geno, clusters = clusters)

## grab intensity data for specific marker(s)
ii <- get.intensity(geno, c("JAX00240610","JAX00240636"))

## grab BAF/LRR
blr <- get.baf(geno.norm, c("JAX00240610","JAX00240636"))
```

Interface to `PLINK`
--------------------

Computation on large SNP array genotyping datasets is not a new problem. Many common operations -- frequency statistics (sample-wise and marker-wise), differentiation statistics (\(F_st\) et al), homozygosity checks, association testing, multivariate clustering by PCA and MDS -- are implemented efficiently in the `PLINK` package. The input formats popularized by `PLINK` are now used by other software in population genetics.

This package provides functions to read and write **binary** `PLINK` filesets. The binary fileset consists of three files:

-   `*.fam`: the \`\`family file'' describing samples (6 columns): family ID, sample ID, mom ID, dad ID, sex (0=unknown, 1=male, 2=female), phenotype (-9=missing)
-   `*.bim`: a file describing the marker map (6 columns, at least): chromosome, marker ID, genetic position (cM), physical position (bp), allele 1, allele 2
-   `*.bed`: compact binary representation of genotypes using 2 bits per genotype

Note that order matters: genotypes from the `*.bed` file are mapped to samples and markers using order of appearance in the `*.fam` and `*.bim` files.

``` {.r}
## this command produces files 'sample.bed', 'sample.bim' and 'sample.fam'
write.plink(geno, "./sample")
## ... and this one reads them back in
geno <- read.plink("/sample")
```

Also included are thin wrappers around some `PLINK` utilities. These of course require a working executable named `plink` in the user's path. [**TODO**]

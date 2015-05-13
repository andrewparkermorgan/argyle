<!-- README.md is generated from README.Rmd. Please edit that file -->




argyle
======

An `R` package for import, QC and (some) analysis of genotyping and hybridization-intensity data from Illumina Infinium arrays.

Dependencies
------------

Effort has been made to keep to a minimum the number of package dependencies, subject to the constraint that I don't want to re-implement from scratch what others have done better.

-   `data.table`: really fast and efficient handling of big (multi-GB scale) table-style data with low overhead
-   `preprocessCore` (Biodoncuctor): robust quantile normalization routine written in `C`
-   `plyr`: generalizations of base-`R`'s `apply()` family
-   `reshape2`: easy "flattening" of matrices to dataframes
-   `digest`: for computing MD5 checksums to check data integrity

The following are required for some functions but one could get by without them:

-   `ggplot2`: required for the plotting functions
-   `corpcor`: required for "fast"-mode PCA

Installation
------------

Installation of the package directly from Github requires `devtools`.

``` {.r}
library(devtools)

## allow R to look for pacakges in both CRAN and Bioconductor
setRepositories(1:2)

## install from Github source
install_github("andrewparkermorgan/argyle")
```

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

Performance
-----------

For a realistic test of `argyle`'s performance, I use an Illumina BeadStudio dataset from the MegaMUGA array for mouse, containing 77808 markers x 96 samples (available from [do.jax.org](http://churchill.jax.org/research/cc/do_data/megamuga/raw/MegaMUGA_22Oct2012/)). The relevant files, compressed with ZIP, have total size 149 MB.

Testing on my Mac Pro desktop (OS X 10.6.8.10K549, 2x2.26 GHz quad-core Xeon, 24 GB DDR3 RAM) under `R` 3.1.3:

    #> Loading argyle
    #> Loading required package: data.table
    #> Loading required package: digest
    #> Loading required package: preprocessCore
    #> Loading required package: reshape2
    #> Loading required package: plyr

``` {.r}
system.time( geno <- read.beadstudio("", snps, "./data/MM_sample") )
#> Reading sample manifest from < ./data/MM_sample/Sample_Map.zip > ...
#> Reading genotypes and intensities from < ./data/MM_sample/FinalReport.zip > ...
#> 
Read 0.0% of 7469568 rows
Read 12.3% of 7469568 rows
Read 24.6% of 7469568 rows
Read 36.9% of 7469568 rows
Read 48.9% of 7469568 rows
Read 61.2% of 7469568 rows
Read 73.5% of 7469568 rows
Read 85.8% of 7469568 rows
Read 98.0% of 7469568 rows
Read 7469568 rows and 11 (of 11) columns from 0.478 GB file in 00:00:12
#> Constructing genotype matrix...
#> Constructing intensity matrices...
#>   77725 sites x 96 samples
#> Done.
#>    user  system elapsed 
#>  33.959   4.771  41.859

summary(geno)
#> --- geno ---
#> A genotypes object with 77725 sites x 96 samples
#> Intensity data: no  
#> Sample metadata: no
#> Filters set: 0 sites / 0 samples 
#> File source: /Users/apm/Dropbox/pmdvlab/argyle/data/MM_sample (on 2015-05-13 14:48:54 )
#> Checksum: bd98520fb7f0dae385eda5d1fa7143c4

print(object.size(geno), units = "Mb")
#> 77.5 Mb

sessionInfo()
#> R version 3.1.3 (2015-03-09)
#> Platform: x86_64-apple-darwin10.8.0 (64-bit)
#> Running under: OS X 10.6.8 (Snow Leopard)
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] argyle_0.1            plyr_1.8.1            reshape2_1.4.1       
#> [4] preprocessCore_1.28.0 digest_0.6.8          data.table_1.9.4     
#> [7] devtools_1.7.0       
#> 
#> loaded via a namespace (and not attached):
#>  [1] chron_2.3-45    evaluate_0.5.5  formatR_1.0     htmltools_0.2.6
#>  [5] knitr_1.9       Rcpp_0.11.5     rmarkdown_0.5.1 roxygen2_4.1.1 
#>  [9] stringr_0.6.2   tools_3.1.3     yaml_2.1.13
```

The final object is, at 78 MB, about half the size of the input. Peak memory usage for parsing BeadStudio output is about 150 MB (as determined by line-by-line profiling not shown here), or on the order of the size of the input. It is worth noting that if the input files were decompressed they would occupy about 500 MB on disk.

Interface to [`R/DOQTL`](http://cgd.jax.org/apps/doqtl/DOQTL.shtml)
-------------------------------------------------------------------

Genotypes processed with `argyle` can be packaged into a set of `R` objects (bundled in an `*.Rdata` file) suitable for use as input to Dan Gatti's `DOQTL` software. `DOQTL` performs haplotype reconstruction and genetic mapping (under both linkage/composite-interval and single-marker association models) in multifounder advanced intercross populations. Its namesake is the Diversity Outbred (DO) mouse population (see [do.jax.org](http://do.jax.org/)).

``` {.r}
## export for DOQTL
export.doqtl(geno, "./doqtl.objects.Rdata")
```

Interface to [`PLINK`](https://www.cog-genomics.org/plink2/)
------------------------------------------------------------

Computation on large SNP array genotyping datasets is not a new problem. Many common operations -- frequency statistics (sample-wise and marker-wise), differentiation statistics (\(F_{st}\) et al), homozygosity checks, association testing, multivariate clustering by PCA and MDS -- are implemented efficiently in the `PLINK` package. The input formats popularized by `PLINK` are now used by other software in population genetics.

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

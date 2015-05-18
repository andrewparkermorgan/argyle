<!-- README.md is generated from README.Rmd. Please edit that file -->




<img src="figs/argyle.png" alt="argyle logo" height=50 /> argyle
================================================================

**An R package for GenotYpes from ILlumina Et al.** Utilities for import, QC and (some) analysis of genotyping and hybridization-intensity data from Illumina Infinium arrays using `R`.

Dependencies
------------

Effort has been made to keep to a minimum the number of package dependencies, subject to the constraint that I don't want to re-implement from scratch what others have done better.

-   `data.table`: really fast and efficient handling of big (multi-GB scale) table-style data with low overhead
-   `preprocessCore` (Biodoncuctor): robust quantile normalization routine written in `C`
-   `plyr`: generalizations of base-`R`'s `apply()` family
-   `reshape2`: easy "flattening" of matrices to dataframes
-   `digest`: for computing MD5 checksums to check data integrity

The following are required for some functions but one could get by without them:

-   `ggplot2` (and friends): required for the plotting functions
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
Read 18.6% of 7469568 rows
Read 37.8% of 7469568 rows
Read 55.8% of 7469568 rows
Read 73.5% of 7469568 rows
Read 91.6% of 7469568 rows
Read 7469568 rows and 11 (of 11) columns from 0.478 GB file in 00:00:09
#> Constructing genotype matrix...
#> Constructing intensity matrices...
#>   77725 sites x 96 samples
#> Done.
#>    user  system elapsed 
#>  25.437   2.732  28.793

summary(geno)
#> --- geno ---
#> A genotypes object with 77725 sites x 96 samples
#> Allele encoding: native 
#> Intensity data: yes (raw) 
#> Sample metadata: yes ( 0 male / 0 female / 96 unknown )
#> Filters set: 0 sites / 0 samples 
#> File source: /Users/apm/Dropbox/pmdvlab/argyle/data/MM_sample (on 2015-05-17 21:40:19 )
#> Checksum: 5dec3557df628860ca53c7192faa838f

print(object.size(geno), units = "Mb")
#> 202.1 Mb

sessionInfo()
#> R version 3.1.2 (2014-10-31)
#> Platform: x86_64-apple-darwin10.8.0 (64-bit)
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] argyle_0.1            plyr_1.8.1            reshape2_1.4         
#> [4] preprocessCore_1.28.0 digest_0.6.4          data.table_1.9.4     
#> [7] devtools_1.6.1       
#> 
#> loaded via a namespace (and not attached):
#>  [1] chron_2.3-45     evaluate_0.5.5   formatR_1.0      htmltools_0.2.6 
#>  [5] knitr_1.8        Rcpp_0.11.3      rmarkdown_0.3.10 stringr_0.6.2   
#>  [9] tools_3.1.2      yaml_2.1.13
```

The final object is 202.1 Mb, and is completely "self-contained" in that sample and marker metadata are stored alongside the genotypes and hybridization intensities. It is worth noting that if the input files were decompressed they would occupy about 500 MB on disk.

Usage
-----

``` {.r}
data(ex)

summary(ex)

## see marker map and sample metadata
head( markers(ex) )
head( samples(ex) )

## subset operations: hard brackets or subset()
ex[ 1:10,1:2 ]

## or stuff like
x <- subset(ex, chr == "chrM")
x <- subset(ex, sex == 2, by = "samples")

## run QC checks and flag samples above thresholds
ex <- run.qc.checks(ex, max.H = 5e3, max.N = 500)
# how many samples fail QC?
summarize.filters(ex)

## peek at a couple of markers
plot.clusters(ex, c("JAX00282035","JAX00282036"))
```

Interface to [`R/DOQTL`](http://cgd.jax.org/apps/doqtl/DOQTL.shtml)
-------------------------------------------------------------------

Genotypes processed with `argyle` can be packaged into a set of `R` objects (bundled in an `*.Rdata` file) suitable for use as input to Dan Gatti's `DOQTL` software. `DOQTL` performs haplotype reconstruction and genetic mapping (under both linkage/composite-interval and single-marker association models) in multifounder advanced intercross populations. Its namesake is the Diversity Outbred (DO) mouse population (see [do.jax.org](http://do.jax.org/)).

``` {.r}
## export for DOQTL
export.doqtl(ex, "./doqtl.objects.Rdata")
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
## this command produces files 'sample.bed', 'sample.bim' and 'sample.fam' in the R sessions temporary directory
ff <- file.path(tempdir(), "sample")
ex <- recode(ex, "native")
write.plink(ex, ff)

## ... and this one reads it back in
summary( read.plink(ff) )
```

Also included are thin wrappers around some `PLINK` utilities. These of course require a working executable named `plink` in the user's path. [**TODO**]

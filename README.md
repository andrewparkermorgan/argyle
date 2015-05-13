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
Read 12.9% of 7469568 rows
Read 25.7% of 7469568 rows
Read 38.6% of 7469568 rows
Read 51.3% of 7469568 rows
Read 64.1% of 7469568 rows
Read 76.8% of 7469568 rows
Read 89.7% of 7469568 rows
Read 7469568 rows and 11 (of 11) columns from 0.478 GB file in 00:00:11
#> Constructing genotype matrix...
#> Constructing intensity matrices...
#>   77725 sites x 96 samples
#> Done.
#>    user  system elapsed 
#>  35.265   5.125  43.969

summary(geno)
#> --- geno ---
#> A genotypes object with 77725 sites x 96 samples
#> Intensity data: yes (raw) 
#> Sample metadata: no
#> Filters set: 0 sites / 0 samples 
#> File source: /Users/apm/Dropbox/pmdvlab/argyle/data/MM_sample (on 2015-05-13 17:32:03 )
#> Checksum: 5dec3557df628860ca53c7192faa838f

print(object.size(geno), units = "Mb")
#> 202.1 Mb

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
#> [7] devtools_1.7.0        knitr_1.9            
#> 
#> loaded via a namespace (and not attached):
#>  [1] chron_2.3-45    evaluate_0.5.5  formatR_1.0     htmltools_0.2.6
#>  [5] Rcpp_0.11.5     rmarkdown_0.5.1 roxygen2_4.1.1  stringr_0.6.2  
#>  [9] tools_3.1.3     yaml_2.1.13
```

The final object is 202.1 Mb, and is completely "self-contained" in that sample and marker metadata are stored alongside the genotypes and hybridization intensities. It is worth noting that if the input files were decompressed they would occupy about 500 MB on disk.

Usage
-----

``` {.r}
## overview of object's contents
summary(geno)
#> --- geno ---
#> A genotypes object with 77725 sites x 96 samples
#> Intensity data: yes (raw) 
#> Sample metadata: no
#> Filters set: 0 sites / 0 samples 
#> File source: /Users/apm/Dropbox/pmdvlab/argyle/data/MM_sample (on 2015-05-13 17:32:03 )
#> Checksum: 5dec3557df628860ca53c7192faa838f

## see marker map and sample metadata
head( markers(geno) )
#>              chr      marker     cM     pos A1 A2
#> UNC6        chr1        UNC6 1.4995 3000355  T  C
#> JAX00000010 chr1 JAX00000010 1.5620 3125499  A  G
#> JAX00240603 chr1 JAX00240603 1.6075 3242877  C  T
#> JAX00240610 chr1 JAX00240610 1.6111 3256689  C  T
#> JAX00240613 chr1 JAX00240613 1.6260 3313481  T  C
#> JAX00240636 chr1 JAX00240636 1.6433 3379644  C  A
head( samples(geno) )
#>                   fid        iid mom dad sex pheno
#> 1-670268     1-670268   1-670268   0   0   0    -9
#> 11-670268   11-670268  11-670268   0   0   0    -9
#> 111-670268 111-670268 111-670268   0   0   0    -9
#> 112-670268 112-670268 112-670268   0   0   0    -9
#> 113-670268 113-670268 113-670268   0   0   0    -9
#> 114-670268 114-670268 114-670268   0   0   0    -9

## subset operations: hard brackets or subset()
geno[ markers(geno)$chr == "chrM",1:3 ]
#>             1-670268 11-670268 111-670268
#> Mit002      "G"      "G"       "T"       
#> Mit001      "G"      "G"       "A"       
#> Mit003      "A"      "A"       "A"       
#> Mit004      "C"      "C"       "C"       
#> Mit005      "C"      "C"       "C"       
#> repMit005   "C"      "C"       "C"       
#> Mit006      "T"      "T"       "T"       
#> Mit007      "C"      "C"       "C"       
#> Mit008      "A"      "A"       "A"       
#> Mit027      "G"      "G"       "A"       
#> Mit009      "C"      "C"       "C"       
#> Mit010      "A"      "A"       "A"       
#> Mit011      "G"      "A"       "A"       
#> JAX00725096 "C"      "C"       "T"       
#> Mit012      "N"      "N"       "N"       
#> JAX00725100 "C"      "C"       "T"       
#> Mit028      "C"      "C"       "T"       
#> JAX00725105 "A"      "A"       "A"       
#> Mit029      "C"      "C"       "T"       
#> Mit013      "C"      "A"       "C"       
#> Mit014C     "C"      "A"       "C"       
#> Mit014G     "A"      "N"       "A"       
#> Mit016      "C"      "C"       "C"       
#> Mit023      "T"      "T"       "T"       
#> Mit024      "T"      "T"       "T"       
#> Mit025      "T"      "T"       "T"       
#> Mit026      "N"      "T"       "T"       
#> Mit017      "G"      "A"       "A"       
#> Mit018      "G"      "A"       "A"       
#> Mit019      "G"      "A"       "G"       
#> Mit020      "T"      "H"       "T"       
#> Mit031      "T"      "H"       "T"       
#> Mit021      "T"      "C"       "N"       
#> Mit022      "G"      "A"       "A"       
#> Mit0040     "T"      "T"       "T"       
#> Mit0041     "C"      "C"       "C"       
#> Mit0042     "T"      "T"       "A"       
#> Mit0043     "C"      "C"       "N"       
#> Mit0044     "C"      "C"       "C"       
#> Mit0045     "H"      "H"       "T"       
#> Mit0046     "A"      "A"       "A"       
#> Mit0047     "T"      "T"       "T"       
#> Mit0048     "C"      "C"       "C"       
#> Mit0049     "T"      "T"       "T"       
#> Mit0050     "G"      "A"       "G"       
#> Mit0051     "T"      "T"       "C"       
#> attr(,"map")
#>              chr      marker cM   pos A1 A2
#> Mit002      chrM      Mit002 NA    54  G  T
#> Mit001      chrM      Mit001 NA    55  G  A
#> Mit003      chrM      Mit003 NA   165  A  C
#> Mit004      chrM      Mit004 NA   348  C  T
#> Mit005      chrM      Mit005 NA   350  C  G
#> repMit005   chrM   repMit005 NA   350  C  G
#> Mit006      chrM      Mit006 NA   572  T  C
#> Mit007      chrM      Mit007 NA   817  C  T
#> Mit008      chrM      Mit008 NA  1584  A  G
#> Mit027      chrM      Mit027 NA  1590  A  G
#> Mit009      chrM      Mit009 NA  1756  C  T
#> Mit010      chrM      Mit010 NA  2185  A  G
#> Mit011      chrM      Mit011 NA  2340  G  A
#> JAX00725096 chrM JAX00725096 NA  2525  T  C
#> Mit012      chrM      Mit012 NA  2615  T  C
#> JAX00725100 chrM JAX00725100 NA  2840  T  C
#> Mit028      chrM      Mit028 NA  2934  T  C
#> JAX00725105 chrM JAX00725105 NA  3220  A  G
#> Mit029      chrM      Mit029 NA  4123  T  C
#> Mit013      chrM      Mit013 NA  4947  C  A
#> Mit014C     chrM     Mit014C NA  7351  A  C
#> Mit014G     chrM     Mit014G NA  7351  A  C
#> Mit016      chrM      Mit016 NA  9461  T  C
#> Mit023      chrM      Mit023 NA  9829  T  G
#> Mit024      chrM      Mit024 NA  9830  T  G
#> Mit025      chrM      Mit025 NA  9831  T  G
#> Mit026      chrM      Mit026 NA  9832  T  G
#> Mit017      chrM      Mit017 NA  9985  G  A
#> Mit018      chrM      Mit018 NA  9985  G  A
#> Mit019      chrM      Mit019 NA 12522  G  A
#> Mit020      chrM      Mit020 NA 12889  T  C
#> Mit031      chrM      Mit031 NA 12889  C  T
#> Mit021      chrM      Mit021 NA 13840  T  C
#> Mit022      chrM      Mit022 NA 15123  G  A
#> Mit0040     chrM     Mit0040 NA 15357  T  C
#> Mit0041     chrM     Mit0041 NA 15396  C  T
#> Mit0042     chrM     Mit0042 NA 15448  T  A
#> Mit0043     chrM     Mit0043 NA 15552  C  T
#> Mit0044     chrM     Mit0044 NA 15843  C  A
#> Mit0045     chrM     Mit0045 NA 15866  T  C
#> Mit0046     chrM     Mit0046 NA 15964  A  G
#> Mit0047     chrM     Mit0047 NA 16010  C  T
#> Mit0048     chrM     Mit0048 NA 16089  T  C
#> Mit0049     chrM     Mit0049 NA 16107  T  C
#> Mit0050     chrM     Mit0050 NA 16201  A  G
#> Mit0051     chrM     Mit0051 NA 16221  C  T
#> attr(,"ped")
#>                   fid        iid mom dad sex pheno
#> 1-670268     1-670268   1-670268   0   0   0    -9
#> 11-670268   11-670268  11-670268   0   0   0    -9
#> 111-670268 111-670268 111-670268   0   0   0    -9
#> attr(,"intensity")
#> attr(,"intensity")$x
#>             1-670268 11-670268 111-670268
#> Mit002         0.130     0.133      2.480
#> Mit001         0.169     0.146      2.593
#> Mit003         2.463     2.668      2.265
#> Mit004         0.323     0.389      0.266
#> Mit005         1.931     2.435      2.172
#> repMit005      1.674     2.303      2.015
#> Mit006         3.307     3.690      3.110
#> Mit007         0.118     0.135      0.112
#> Mit008         1.830     1.827      1.798
#> Mit027         0.069     0.062      1.488
#> Mit009         0.052     0.036      0.053
#> Mit010         3.073     3.441      2.833
#> Mit011         0.138     2.926      2.379
#> JAX00725096    0.112     0.127      2.142
#> Mit012         2.902     3.053      2.642
#> JAX00725100    0.370     0.371      2.986
#> Mit028         0.145     0.188      2.101
#> JAX00725105    2.489     2.616      2.174
#> Mit029         0.231     0.248      2.316
#> Mit013         0.092     1.166      0.064
#> Mit014C        0.167     3.461      0.213
#> Mit014G        2.772     0.630      2.676
#> Mit016         0.156     0.130      0.149
#> Mit023         1.846     1.952      1.657
#> Mit024         1.649     1.859      1.602
#> Mit025         0.400     1.894      1.598
#> Mit026         0.117     1.034      0.971
#> Mit017         0.105     2.911      2.580
#> Mit018         0.087     2.902      2.369
#> Mit019         0.244     2.799      0.304
#> Mit020         2.557     0.431      2.011
#> Mit031         2.620     0.409      2.045
#> Mit021         2.738     0.470      0.630
#> Mit022         0.407     3.532      2.650
#> Mit0040        2.475     2.731      2.406
#> Mit0041        0.040     0.061      0.056
#> Mit0042        2.626     2.711      0.260
#> Mit0043        0.109     0.161      0.094
#> Mit0044        0.370     0.525      0.291
#> Mit0045        1.151     1.052      2.685
#> Mit0046        2.941     3.248      2.588
#> Mit0047        2.693     2.980      2.539
#> Mit0048        0.064     0.065      0.041
#> Mit0049        1.377     1.592      1.221
#> Mit0050        0.170     2.730      0.193
#> Mit0051        2.869     2.539      0.118
#> 
#> attr(,"intensity")$y
#>             1-670268 11-670268 111-670268
#> Mit002         1.937     1.968      0.104
#> Mit001         1.841     1.972      0.097
#> Mit003         0.000     0.000      0.054
#> Mit004         1.984     2.138      1.962
#> Mit005         0.631     0.586      0.713
#> repMit005      0.551     0.483      0.599
#> Mit006         0.129     0.222      0.115
#> Mit007         1.766     1.797      1.742
#> Mit008         0.067     0.071      0.069
#> Mit027         1.125     1.398      0.073
#> Mit009         1.681     1.685      1.674
#> Mit010         0.288     0.377      0.300
#> Mit011         1.922     0.164      0.116
#> JAX00725096    1.499     1.604      0.080
#> Mit012         0.120     0.229      0.122
#> JAX00725100    2.200     2.481      0.163
#> Mit028         1.604     1.731      0.202
#> JAX00725105    0.212     0.285      0.245
#> Mit029         1.741     1.958      0.108
#> Mit013         1.065     0.016      0.889
#> Mit014C        0.580     0.652      0.599
#> Mit014G        0.422     0.390      0.460
#> Mit016         1.608     1.562      1.598
#> Mit023         0.018     0.011      0.025
#> Mit024         0.006     0.000      0.006
#> Mit025         0.031     0.007      0.023
#> Mit026         0.028     0.017      0.022
#> Mit017         2.111     0.535      0.432
#> Mit018         1.913     0.550      0.329
#> Mit019         1.717     0.396      1.979
#> Mit020         0.077     1.701      0.143
#> Mit031         0.087     1.631      0.171
#> Mit021         0.192     1.846      0.427
#> Mit022         2.026     0.208      0.148
#> Mit0040        0.009     0.000      0.009
#> Mit0041        1.541     1.694      1.641
#> Mit0042        0.223     0.201      2.076
#> Mit0043        1.576     1.789      0.028
#> Mit0044        1.451     1.540      1.327
#> Mit0045        1.724     1.878      0.257
#> Mit0046        0.480     0.543      0.456
#> Mit0047        0.060     0.180      0.106
#> Mit0048        1.104     1.148      0.974
#> Mit0049        0.030     0.000      0.000
#> Mit0050        2.031     0.282      2.009
#> Mit0051        0.073     0.070      1.244
#> 
#> attr(,"normalized")
#> [1] FALSE
#> attr(,"alleles")
#> [1] "native"
#> attr(,"filter.sites")
#>      Mit002      Mit001      Mit003      Mit004      Mit005   repMit005 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>      Mit006      Mit007      Mit008      Mit027      Mit009      Mit010 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>      Mit011 JAX00725096      Mit012 JAX00725100      Mit028 JAX00725105 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>      Mit029      Mit013     Mit014C     Mit014G      Mit016      Mit023 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>      Mit024      Mit025      Mit026      Mit017      Mit018      Mit019 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>      Mit020      Mit031      Mit021      Mit022     Mit0040     Mit0041 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>     Mit0042     Mit0043     Mit0044     Mit0045     Mit0046     Mit0047 
#>       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
#>     Mit0048     Mit0049     Mit0050     Mit0051 
#>       FALSE       FALSE       FALSE       FALSE 
#> attr(,"filter.samples")
#>   1-670268  11-670268 111-670268 
#>      FALSE      FALSE      FALSE 
#> attr(,"class")
#> [1] "genotypes" "matrix"
## or stuff like
x <- subset(geno, chr == "chrM")
x <- subset(geno, sex == 2, by = "samples")

## run QC checks and flag samples above thresholds
geno <- run.qc.checks(geno, max.H = 60e3, max.N = 5e3)
#> Performing QC checks on genotype calls...
#> Recoding to 0/1/2 using reference alleles.
#> Performing QC checks on hybridization intensities...
#> 0 markers and 16 samples now flagged as low-quality.
# how many samples fail QC?
summarize.filters(geno)
#>   sites samples 
#>       0      16
## graphical summary
qcplot(geno)
```

![](README-unnamed-chunk-5-1.png)

``` {.r}

## peek at a couple of markers
plot.clusters(geno, 1:3)
#> Recoding to 0/1/2 using reference alleles.
```

![](README-unnamed-chunk-5-2.png)

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
geno <- read.plink("./sample")
```

Also included are thin wrappers around some `PLINK` utilities. These of course require a working executable named `plink` in the user's path. [**TODO**]

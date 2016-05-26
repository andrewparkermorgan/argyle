## genetics.R
## implement simple genetics logic on 'genotypes' object

## for input checking, but not used to avoid repeated copying of big objects
.check.and.warn <- function(x) {
	
	if (is.matrix(x)) {
		if (!inherits(x, "genotypes")) {
			warning(paste("Input is not a 'genotypes' object; assuming it is a matrix of",
						  "genotypes with shape (markers x samples)."))
			if (!is.numeric(x)) {
				stop("Input is not numeric.")
			}
			else {
				return(TRUE)
			}
		}
		else {
			if (!is.numeric(x)) {
				stop("Input is not numeric.")
			}
			else {
				return(TRUE)
			}
		}
	}
	else {
		stop("Input is not a matrix.")
	}
	
}

#' Calculate allele frequencies by marker or sample
#' 
#' @param gty a \code{genotypes} object with numeric allele encoding
#' @param by margin over which to calculate frequencies (usually \code{"markers"})
#' @param na.rm skip over missing genotypes in frequency calculations
#' @param counts logical; if \code{TRUE}, report allele counts instead of relative frequencies
#' @param ... ignored
#' 
#' @return a named vector of allele frequencies
#' 
#' @export
freq <- function(gty, by = c("markers","samples"), na.rm = TRUE, counts = FALSE, ...) {
	
	#if (!.check.and.warn(gty))
	#	stop("Can't understand input.")
	
	if (!is.numeric(gty))
		stop("Input is not numeric.")
	
	by <- match.arg(by)
	if (by == "markers")
		if (!counts)
			return( rowMeans(gty, na.rm = na.rm)/2 )
		else
			return( rowSums(gty, na.rm = na.rm) )
	else
		if (!counts)
			return( colMeans(gty, na.rm = na.rm)/2 )
		else
			return( colSums(gty, na.rm = na.rm) )
	
}



#' Calculate minor-allele frequency (MAFs) from a genotypes matrix
#' 
#' Just an alias or \code{freq(x, "markers")}.  Sticklers for nomenclature will be aware that this
#' function does not truly return the MAF unless the allele encoding is \code{"relative"}.
#' 
#' @param ... passed through to call to \code{freq()}
#' 
#' @seealso \code{\link{freq}}
#'
#' @export
maf <- function(...) freq(..., by = "markers")

#' Count number of missing genotypes by marker or sample
#'
#' @param gty a \code{genotypes} object
#' @param by margin over which to count missing calls
#' @param ... ignored
#'
#' @return a named vector of missing-genotype *counts*
#'
#' @export
nmiss <- function(gty, by = c("markers","samples"), ...) {
	
	#if (!.check.and.warn(gty))
	#	stop("Can't understand input.")
	
	by <- match.arg(by)
	if (by == "markers")
		dim.along <- 1
	else
		dim.along <- 2
	
	apply(gty, dim.along, function(x) sum(is.na(x)))
	
}

#' Calculate rate of missing genoypes by marker or sample
#' 
#' @param gty a \code{genotypes} object
#' @param by margin over which to count missing calls
#' @param ... ignored
#' 
#' @return a named vector of missing-genotype *proportions*
#' 
#' @export
missingness <- function(gty, by = c("markers","samples"), ...) {
	
	rez <- nmiss(gty, by = by, ...)
	return( rez/length(rez) )
	
}

#' Identify the consensus genotype call among a group of samples
#' 
#' @param gty a \code{genotypes} object
#' @param nas.allowed maximum proportion of missing genotypes to ignore 
#' @param ... ignored
#' 
#' @return a named vector of consensus genotypes, with length equal to number of rows in \code{gty}
#' 
#' @details Define the "consensus" genotype to be the genotype with greatest frequency among a group
#' 	of samples, with ties broken in favor of the major allele.  This function computes that consensus
#' 	ignoring missing values, and then masks sites with missingness greater than \code{nas.allowed}.
#' 
#' @export
consensus <- function(gty, nas.allowed = 0.0, ...) {
	
	#if (!.check.and.warn(gty))
	#	stop("Can't understand input.")
	
	rez <- apply(gty, 1, function(x) which.max(tabulate(x+1, nbins = 3))-1)
	nas <- (missingness(gty, "markers") > nas.allowed)
	rez[nas] <- NA
	
	return(rez)
	
}

## helper function to determine if a marker is segregating
is.segregating <- function(x, allow.het = TRUE, ...) {
	
	tbl <- tabulate(x+1, nbins = 3)
	if (allow.het) {
		flag <- sum(tbl > 0) > 1
		flag <- flag | any(x == 1, na.rm = TRUE)
	}	
	else {
		flag <- sum(tbl[c(1,3)] > 0) > 1
	}
	return(flag)
	
}

#' Identify segregating sites among a group of samples
#' 
#' @param gty a \code{genotypes} object
#' @param nas.allowed maximum proportion of missing genotypes to ignore
#' @param allow.het logical; if \code{FALSE}, ignore samples with heterozygous genotype
#' @param ... ignored
#' 
#' @return a named logical vector, of length equal to the number of rows in \code{gty},
#' 	indicating whether a site is segregating
#' 	
#' @details A site is defined to be "segregating" if two or more alleles are present in one or
#' 	more copies across all samples.  Sites at which missingness is greater than \code{nas.missing}
#' 	are set to \code{NA}.
#' 	
#' 	When \code{allow.het = FALSE}, this is equivalent to asking whether at least one pair of samples
#' 	has a fixed difference at each site.
#' 	
#' @seealso \code{\link{fixed.diffs}}
#' 
#' @export
segregating <- function(gty, nas.allowed = 0.0, allow.het = TRUE, ...) {
	
	#if (!.check.and.warn(gty))
	#	stop("Can't understand input.")
	
	rez <- apply(gty, 1, is.segregating, allow.het = allow.het)
	nas <- (missingness(gty, "markers") > nas.allowed)
	rez[nas] <- NA
	
	return(rez)
	
}

#' Identify fixed differences between a pair of samples
#' 
#' @param gty a \code{genotypes} object
#' @param cols indexing vector for selecting columns from \code{gty}
#' @param ... ignored
#' 
#' @return a named logical vector, of length equal to the number of rows in \code{gty},
#' 	indicating whether a site is a fixed difference between the two samples
#' 	
#' @details A site is defined to have a "fixed difference" between two samples if they have opposite
#' 	homozygous genotypes.  In contrast to \code{segregating()}, this function rejects sites at which either
#' 	sample is heterozygous or has a missing genotype, rather than marking them with \code{NA}s.
#' 	
#' 	If more than two samples are present in the input, an error occurs.
#' 	
#' @seealso \code{\link{segregating}}
#' 
#' @export
fixed.diffs <- function(gty, cols = TRUE, ...) {
	
	#if (!.check.and.warn(gty))
	#	stop("Can't understand input.")
	
	gty <- gty[ ,TRUE, drop = FALSE ]
	if (ncol(gty) != 2)
		stop("A 'fixed difference' can only exist between exactly two samples.")
	
	diffs <- (gty[ ,1 ] != gty[ ,2 ])
	ns <- is.na(rowSums(gty))
	hs <- (gty[ ,1 ] == 1) | (gty[ ,2 ] == 1)
	
	return(setNames( as.vector(diffs & !ns & !hs), rownames(gty) ))
	
}

#' Find markers informative between sample groups
#' 
#' @param gty a \code{genotypes} object
#' @param between expression for grouping samples (see Details)
#' @param ... other parameters to be passed to \code{consensus()} (such as \code{nas.allowed})
#' 
#' @return a named logical vector, of length equal to the number of rows in \code{gty},
#' 	indicating whether a marker is informative in the given samples
#' 	
#' @details A marker is defined to be "informative" if it is not monomorphic across a group of samples.
#' 	In the context of a genetic cross (eg. an F2 intercross), an informative marker is one with a fixed
#' 	difference between the two parental lines.  This function generalizes that idea to more than two groups:
#' 	a marker is informative if it has a fixed difference between any pair of groups, as defined in the
#' 	\code{between} expression.
#' 	
#'  The expression \code{between} can be protected for evaluation in the environment of the sample metadata
#'  in \code{attr(gty, "ped")} (view it with \code{samples(gty)}) by wrapping it in \code{.()}.  To group
#'  on more than one variable, provide multiple expressions in a list or separate them with commas like
#'  \code{.(expr1, expr2)}.  (This syntax will be familiar to users of the \code{plyr} package.)
#'  
#'  Consensus genotypes are first computed within the groups defined by unique values of \code{between},
#'  and then markers with fixed differences between any two groups are idenfitied.
#' 	
#' @seealso \code{\link{segregating}}, \code{\link{fixed.diffs}}
#' 
#' @export
informative <- function(gty, between, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	#message("Calculating consensus genotypes...")
	cons <- do.call("cbind", genoapply(gty, 2, between, consensus, ..., strip = TRUE))
	
	#message("Iterating over groups...")
	seg <- segregating(cons, allow.het = FALSE, ...)
	
	return(seg)
	
}

#' Predict genotype of an F1 individual given genotypes of its parents
#' 
#' @param object a \code{genotypes} object
#' @param na.rm logical; should missing genotypes be ignored?
#' @param ... ignored
#' 
#' @return a new \code{genotypes} object containing predicted genotypes for all possible F1s between
#' 	the parents in \code{object}
#' 	
#' @details Each column in the input is assumed to be either a representative individual from an inbred
#' 	line, or a consensus genotype across several such individuals.  No checks are performed on the homzygosity
#' 	of the parents, but parents with high heterozygosity will result in F1 predictions with mostly missing
#' 	genotypes.
#' 	
#' 	When \code{na.rm = TRUE}, a missing genotype in either parent causes an \code{NA} to be emitted for
#' 	the F1; when \code{na.rm = FALSE}, the non-missing parental genotype (if any) is emitted. Heterozygous
#' 	calls in either parent always generate an \code{NA} in the F1.
#' 
#' 	F1s are named like "parent1::parent2".
#' 
#' 	Sex chromosomes and mitochonria receive no special treatment, so genotypes there will probably be bogus
#' 	(except for the case of chrX in F1 females.)
#' 	
#' @seealso \code{\link{segregating}}, \code{\link{fixed.diffs}}
#' 
#' @export predict.f1
predict.f1 <- function(object, na.rm = FALSE, ...) {
	
	if (!inherits(object, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	pairs <- combn(ncol(object), 2)
	f1.names <- paste(colnames(object)[ pairs[1,] ], colnames(object)[ pairs[2,] ],
					  sep = "::")
	rez <- matrix(NA, ncol = ncol(pairs), nrow = nrow(object),
				  dimnames = list(rownames(object), f1.names))
	
	map <- attr(object, "map")
	object <- .copy.matrix.noattr(object)
	message(paste("Predicting F1 genotypes for", ncol(pairs), "pairs of parents..."))
	for (i in seq_len(ncol(pairs))) {
		f1 <- as.vector(rowMeans(object[ ,pairs[,i] ], na.rm = na.rm))
		het <- apply(object, 1, function(x) any(x == 1, na.rm = na.rm))
		f1[het] <- NA
		rez[,i] <- f1
	}
	
	rez <- genotypes(rez, map = map, alleles = "01")
	return(rez)
	
}

pairwise.diffs <- function(gty, ...) {
	
	pairs <- combn(ncol(gty), 2)
	pnames <- paste(colnames(gty)[ pairs[1,] ], colnames(gty)[ pairs[2,] ],
					 sep = "::")
	rez <- matrix(NA, ncol = ncol(pairs), nrow = nrow(gty),
				  dimnames = list(rownames(gty), pnames))
	
	map <- attr(gty, "map")
	gty <- .copy.matrix.noattr(gty)
	message(paste("Calculating pairwise differences for", ncol(pairs), "pairs of parents..."))
	for (i in seq_len(ncol(pairs))) {
		is.diff <- as.integer(gty[ pairs[1,i] ] == gty[ pairs[2,i] ])
		rez[,i] <- is.diff
	}
	
	rez <- genotypes(rez, map = map, alleles = "01")
	return(rez)
	
}

#' Calculate the proportion of sites at which each sample is heterozygous
#' 
#' @param gty a \code{genotypes} object
#' @param na.rm logical; ignore missing genotypes when counting heterozygous sites?
#' @param ... ignored
#' 
#' @return a named vector with the proportion of sites of heterozygous for each sample
#' 
#' @export
prop.het <- function(gty, na.rm = TRUE, ...) {
	
	gty[ gty != 1 ] <- 0
	colMeans(gty, na.rm = na.rm)
	
}

#' Calculate heterozygosity by marker
#' 
#' @param gty a \code{genotypes} object
#' @param hwe if \code{TRUE}, compute the MLE \code{hat{h}} assuming Hardy-Weinberg equilibrium
#' @param na.rm logical; ignore missing genotypes?
#' @param ... ignored
#' 
#' @return a named vector with the heterozygosity of each marker
#'
#' @export
heterozygosity <- function(gty, hwe = FALSE, na.rm = TRUE, ...) {
	
	
	if (!hwe) {
		gty[ gty != 1 ] <- 0
		rowMeans(gty, na.rm = na.rm)
	}
	else {
		apply(gty, 1, function(x) 1 - mean(x == 0, na.rm = na.rm)^2 - mean(x == 2, na.rm = na.rm))
	}
		
	
}

## cheap estimate of pairwise LD via genotypic correlation
ld <- function(gty, force = FALSE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (prod(dim(gty)) > 1e6 && !force)
		stop("Resulting matrix will have >1M entries. Use force = TRUE to proceed anyway.")
	
	ldmat <- cor(t(gty))^2
	class(ldmat) <- c("LD.matrix", class(ldmat))
	return(ldmat)
	
}

#' Compute rate of Mendelian inconsistency between an individual and possible mother-father pairs.
#' 
#' @param gty a \code{genotypes} object containing exactly one target individual
#' @param parents a \code{genotypes} object containing possible parents for the target individual
#' @param verbose logical; if \code{TRUE}, send progress messages to terminal
#' @param ... ignored
#' 
#' @return a dataframe of mother-father pairs and corresponding Mendel distances
#' 
#' @details Define the "Mendel distance" between an individual and two possible parents as the mean
#' 	rate of Mendelian inconsistencies -- that is, markers at which zero alleles are shared IBS between
#' 	one parent and one offspring -- across the pairs (father, offspring) and (mother, offspring). Given
#' 	a set of possible parents (in \code{parents}), this function creates all possible male-female pairs
#' 	and computes the Mendel distance between each pair and a single target individual (in \code{gty}).
#' 	
#' 	If the sex of the target individual is non-missing, the sex chromosomes will be handled specially:
#' 	for males, chrY is compared to the father's chrY and chrX to the mother's chrX; for females, chrX
#' 	is compared to both parents and chrY is ignored.  If the sex of the target individual is missing,
#' 	the sex chromosomes will be ignored, and mitochondrial markers are always ignored because they are
#' 	so few in number and tend to be homoplasic.
#' 
#' @seealso \code{\link{guess.parents}}
#' 
#' @export
mendel.distance <- function(gty, parents, verbose = TRUE, ...) {
	
	if (!inherits(gty, "genotypes") || !inherits(parents, "genotypes"))
		stop("Both input (gty) and possible parents (parents) must be of class 'genotypes'.")
	
	if (attr(gty, "alleles") != attr(parents, "alleles") || !is.numeric(gty))
		stop("Alleles should be in same numeric encoding in unknown (gty) and parents (parents).")
	
	if (ncol(gty) > 1)
		stop("Can only handle one sample of unknown parentage at a time.")
	
	## extract parents of each sex
	moms <- parents[ ,attr(parents, "ped")$sex == 2 ]
	dads <- parents[ ,attr(parents, "ped")$sex == 1 ]
	if (nrow(moms) < 1 || nrow(dads) < 1)
		stop("Not enough parents.")
	
	## keep only markers shared between parents and offspring...
	keep <- intersect(rownames(gty), rownames(moms))
	gty <- gty[ keep, ]
	moms <- moms[ keep, ]
	dads <- dads[ keep, ]
	
	## generate possible pairs
	#pairs <- as.matrix(expand.grid( colnames(moms), colnames(dads), stringsAsFactors = FALSE ))
	if (verbose)
		message("Investigating ", nrow(pairs), " possible parents...")
	
	## get chromosomes by inheritance pattern
	autos <- autosomes(gty[,1])
	chrx <- xchrom(gty[,1])
	chry <- ychrom(gty[,1])
	
	## weighting factors
	A <- nrow(autos)
	X <- nrow(chrx)
	Y <- nrow(chry)
	
	## autosomal distances
	dm <- apply(autosomes(moms), 2, ibs0, autos)
	dd <- apply(autosomes(dads), 2, ibs0, autos)
	
	if (sex(gty) == 1) {
		## males: chrX distances to mom only, plus chrY to dad
		dm <- (A/(A+X))*dm + (X/(A+X))*apply(xchrom(moms), 2, ibs0, chrx)
		dd <- (A/(A+Y))*dd + (Y/(A+Y))*apply(ychrom(dads), 2, ibs0, chry) 
	}
	else if (sex(gty) == 2) {
		## females: chrX distances to both parents
		dm <- (A/(A+X))*dm + (X/(A+X))*apply(xchrom(moms), 2, ibs0, chrx)
		dd <- (A/(A+X))*dd + (X/(A+X))*apply(xchrom(dads), 2, ibs0, chrx)
	}
	
	rez <- data.frame()
	if (ncol(moms))
		rez <- rbind(rez, data.frame(iid = colnames(gty), parent = colnames(moms), sex = 2, score = dm[ colnames(moms) ]))
	
	if (ncol(dads))
		rez <- rbind(rez, data.frame(iid = colnames(gty), parent = colnames(dads), sex = 2, score = dd[ colnames(dads) ]))
		
	return(rez)
	
}

ibs0 <- function(x,y) sum(abs(x-y) == 2, na.rm = TRUE)/sum(!is.na(x+y))

ibs0.fancy <- function(x, y, ...) {
	
	if (!all(inherits(x, "genotypes"), inherits(y, "genotypes")))
		stop("Both input objects must be of class 'genotypes'.")
	
	## make sure only one sample in each object
	x <- x[,1]
	y <- y[,1]
	
	## get sexes
	sx <- sex(x)
	sy <- sex(y)
	
	## autosomal genotypes
	xa <- autosomes(x)
	ya <- autosomes(y)
	
	## chrX genotypes
	xx <- x[ grepl("X", attr(x, "map")$chr), ]
	yx <- y[ grepl("X", attr(x, "map")$chr), ]
	
	## chrY genotypes
	xy <- x[ grepl("Y", attr(x, "map")$chr), ]
	yy <- y[ grepl("Y", attr(x, "map")$chr), ]
	
	## chrM genotypes
	xm <- x[ grepl("M", attr(x, "map")$chr), ]
	ym <- y[ grepl("M", attr(x, "map")$chr), ]
	
	## autosomal distance
	da <- sum(abs(xa-ya) == 2, na.rm = TRUE)
	na <- sum(!is.na(xa+ya))
	
	## sex chr distances
	dx <- sum(abs(xx-yx))
	nx <- sum(!is.na(xx+yx))
	if (sx == 1 && sy == 1) {
		## only check chrY if both samples male
		dy <- sum(abs(xy-yy) == 2, na.rm = TRUE)
		ny <- sum(!is.na(xy+yy))
	}
	else {
		dy <- 0
		ny <- 0
	}
	
	## mitochondria distances
	if (sx == 2 && sy == 2) {
		## only check mitochondria if both samples female
		dm <- sum(abs(xm-ym) == 2, na.rm = TRUE)
		nm <- sum(!is.na(xm+ym))
	}
	else {
		dm <- 0
		nm <- 0
	}
	
	return( (dx+dy+dm+da)/(nx+ny+nm+na) )

}

#' Attempt to guess the mother-father pair corresponding to offspring
#' 
#' @param gty a \code{genotypes} object containing target individuals
#' @param parents a \code{genotypes} object containing possible parents for the target individuals
#' @param ... ignored
#' 
#' @return the sample metedata from the input object, with \code{mom} and \code{dad} values replaced
#' 	by the function's best guesses
#' 
#' @details This wrapper calls \code{mendel.distance()} for all mother-father pairs available from
#' 	\code{parents} and all offspring in \code{gty}, and for each offspring reports the mother-father
#' 	pair with the smallest Mendel distance.  No guarantees are made with respect to near-misses.
#' 
#' @export
guess.parents <- function(gty, parents, ...) {
	## TODO: fix this
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	message("Attempting to guess parents of ", ncol(gty), " samples.")
	
	fam <- attr(gty, "ped")
	for (i in colnames(gty)) {
		message("\t... ", i)
		scores <- mendel.distance(gty[,i], parents, verbose = FALSE)
		guess <- scores[ which.min(scores$score),c("mom","dad") ]
		fam[ i,"mom" ] <- guess[1,"mom"]
		fam[ i,"dad" ] <- guess[1,"dad"]
	}
	
	message("Done.\n")
	return(fam)
	
}

#' Calculate a simple genetic distance: proportion of alleles shared IBS
#' 
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' 
#' @return a \code{dist} object whose entires are the proportion of alleles shared
#' 	identical-by-state between samples
#' 	
#' @details To get a plain distance matrix, do \code{as.matrix(dist(x))}.
#' 
#' @seealso \code{\link[stats]{dist}}
#' 
#' @export
dist.genotypes <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (!is.numeric(gty))
		gty <- recode(gty, "01")
	
	message("Computing distance matrix...")
	## compute distance matrix with Rcpp, super fast
	d <- dist_ibs(gty)
	dimnames(d) <- list(colnames(gty), colnames(gty))
	return( as.dist(d) )
	
}
#' @export
dist <- function(x, ...) UseMethod("dist")

weir.fst <- function(gty, subpop = NULL, per.locus = FALSE, ...) {
	
	## Copyright Eva Chan 2008
	## eva@evachan.org
	##
	## A script to estimate the variance components and fixation indices as described in
	## Weir & Cockerham 1984 Evolution 38(6): 1358-1370.
	##
	## Arguments
	## =========
	## geno:   matrix of genotypes with rows corresp. to markers and columns to individuals; 
	##         notation for genotyeps are {0,1,2} indicating the number of one of the two alleles
	## subpop: vector indicting the sub-popln to which the individuals belong to
	##
	## Output
	## ======
	## list of two objects: perloc and global
	## perloc: matrix of 6 columns and as many rows as markers in geno
	##         the 6 columns contain the estiamted variance components and fixation indices per locus
	##         a = component of variance between subpops
	##         b = component of variance between individuals within subpops
	##         c = component of variance between gametes within individuals
	## global: numeric vector of three values corresponding to the esimated F (Fit), theta (Fst), & 
	##         f (Fis) across all loci
	##
	## Note
	## ====
	## R/HIERFSTAT also estimate F-statistics using variance component estimation.
	## Results from that package is not too different to those from this function; I suspect
	## there are two sources of differences:
	## 1) all estimates of variance component from HIERFSTAT are doubled in magnitude to those
	##    from this function (i.e. scaled by factor of 2);
	## 2) rounding off variations may also be present.
	## The scaled difference in the estimates of variance components poses no problem when
	## calcualting fixation indicies as the values scaling factor is cancelled out in the
	## calculation of the ratios.
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	if (!is.numeric(gty))
		geno <- recode(gty, "01")
	else
		geno <- gty
	
	if (is.null(subpop))
		subpop <- as.character(attr(geno, "ped")$fid)
	else
		if (length(subpop) != ncol(geno))
			stop("Length of population labels doesn't match dimension of genotypes matrix.")
	
	spop <- unique(as.character(subpop))    ## unique spops
	r <- length(spop)
	message("Calculating Fst with following population labels:\n",
			paste(spop, collapse = ", "))
	
	n11 <- n12 <- n22 <- matrix(NA, ncol=r, nrow=nrow(geno))
	for(i in 1:r) {
		inds <- which(subpop == spop[i])
		n11[,i] <- rowSums(geno[,inds] == 0,na.rm = TRUE)
		n12[,i] <- rowSums(geno[,inds] == 1,na.rm = TRUE)
		n22[,i] <- rowSums(geno[,inds] == 2,na.rm = TRUE)
	}
	ni <- n11 + n12 + n22
	pi_tilda <- ((2 * n11) + n12) / (2 * ni)
	hi_tilda <- n12 / ni
	n_bar <- rowSums(ni, na.rm = TRUE)/r
	#      C_square <- ( apply(ni*ni,1,sum,na.rm=T) - (n_bar*n_bar*r) ) / ( (n_bar*n_bar) * (r-1) )    ## mod 2/3/2008
	#      nc <- n_bar * (1 - (C_square/r))
	nc <- ((r*n_bar) - rowSums(((ni*ni)/(r*n_bar)), na.rm = TRUE)) / (r - 1)
	p_bar <- rowSums( (ni*pi_tilda)/(r*n_bar), na.rm = TRUE )
	s_square <- rowSums( (ni*((pi_tilda-p_bar)^2)) / ((r-1)*n_bar), na.rm = TRUE )
	h_bar <- rowSums((ni*hi_tilda)/(r*n_bar), na.rm = TRUE )
	
	#      F_hat = 1 - ( (h_bar*(1-(C_square/r))) / ( (2*p_bar*(1-p_bar)*(1-((n_bar*C_square)/(r*(n_bar-1))))) + (2*(s_square/r)*(1+(((r-1)*(n_bar*C_square))/(r*(n_bar-1))))) + ((h_bar/2)*(C_square/(r*(n_bar-1)))) ))
	#      theta_hat <- (s_square - ((1/(n_bar-1))*((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - (h_bar/4)))) / (((1-((n_bar*C_square)/(r*(n_bar-1))))*p_bar*(1-p_bar)) + ((1+(((r-1)*n_bar*C_square)/(r*(n_bar-1))))*(s_square/r)) + ((C_square/(r*(n_bar-1)))*(h_bar/4)))
	#      f_hat <- 1 - (h_hat / ((((2*n_bar)/(n_bar-1))*p_bar*(1-p_bar)) - (((2*n_bar*(r-1))/(r*(n_bar-1)))*s_square) - ((1/(n_bar-1))*(h_bar/2))))
	
	a_hat <- (n_bar/nc) * ( s_square - ((1/(n_bar-1))*((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - ((1/4)*h_bar))) )
	b_hat <- (n_bar/(n_bar-1)) * ((p_bar*(1-p_bar)) - (((r-1)/r)*s_square) - ((((2*n_bar)-1)/(4*n_bar))*h_bar))
	c_hat <- h_bar/2
	
	F_hat <- 1 - (c_hat / (a_hat + b_hat + c_hat))
	theta_hat <- a_hat / (a_hat + b_hat + c_hat)
	f_hat <- 1 - (c_hat / (b_hat + c_hat))
	
	F_hat_w <- 1 - (sum(c_hat,na.rm=T) / sum((a_hat + b_hat + c_hat),na.rm=T))
	theta_hat_w <- sum(a_hat,na.rm=T) / sum((a_hat + b_hat + c_hat),na.rm=T)
	f_hat_w <- 1 - (sum(c_hat,na.rm=T) / sum((b_hat + c_hat),na.rm=T))
	
	if (per.locus)
		list( perloc=cbind(a_hat=a_hat, b_hat=b_hat, c_hat=c_hat, F_hat=F_hat, theta_hat=theta_hat, f_hat=f_hat),
		 	 global=c(F_hat=F_hat_w, theta_hat=theta_hat_w, f_hat=f_hat_w) )
	else
		return(theta_hat_w)
	
}


#' Thin markers to a specified density on the genetic map
#'
#' @param gty a \code{genotypes} object
#' @param spacing numeric; target distance between retained markers, in cM
#' @param ... ignored
#' 
#' @return a new \code{genotypes} object with fewer markers
#' 
#' @aliases thin
#' @export
thin.genotypes <- function(gty, spacing = 5.0, ...) {
	
	## start borrowing from ggplot2
	find.zero <- function (x_range, width, boundary) {
		shift <- floor((x_range[1] - boundary)/width)
		boundary + shift * width
	}
	fixed.bins <- function (x, width, center = NULL, boundary = NULL, closed = c("right","left")) {
		x <- as.numeric(x)
		width <- as.numeric(width)
		closed <- match.arg(closed)
		x_range <- range(x, na.rm = TRUE, finite = TRUE)
		if (length(x_range) == 0) {
			return(x)
		}
		if (!is.null(boundary) && !is.null(center)) {
			stop("Only one of 'boundary' and 'center' may be specified.")
		}
		if (is.null(boundary)) {
			if (is.null(center)) {
				boundary <- width/2
			}
			else {
				boundary <- center - width/2
			}
		}
		boundary <- as.numeric(boundary)
		min_x <- find.zero(x_range, width, boundary)
		max_x <- max(x, na.rm = TRUE) + (1 - 1e-08) * width
		breaks <- seq(min_x, max_x, width)
		cut(x, breaks, include.lowest = TRUE, right = (closed == "right"))
	}
	## end borrowing
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (!.has.valid.map(gty))
		stop("Input object must have a valid genetic map.")
	
	# subsample by chromosome
	message("Starting with ", nrow(gty), " markers...")
	suppressMessages({
		thinned <- genoapply(subset(gty, !is.na(cM)), 1, .(chr), function(d) {
			binned <- fixed.bins(markers(d)$cM, spacing)
			.keep1 <- function(z) z[1]
			keep <- sapply(split(seq_len(nrow(d)), binned), .keep1)
			keep <- union(keep, c(1,nrow(d)))
			d[ sort(keep), ]
		})
		thinned <- Reduce(rbind, thinned)
	})
	
	dropped <- setdiff(attr(gty, "map")$chr, attr(thinned, "map")$chr)
	
	message("Thinned to ", nrow(thinned), " markers.")
	if (length(dropped))
		message("Dropped these chromosomes without genetic positions: ", paste(dropped, collapse = ", "))
	return(thinned)
	
}
thin <- function(gty, ...) UseMethod("thin")
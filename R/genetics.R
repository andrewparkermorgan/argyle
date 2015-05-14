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
#' 
#' @return a named vector of allele frequencies
#' 
#' @export
freq <- function(gty, by = c("markers","samples"), na.rm = TRUE, ...) {
	
	#if (!.check.and.warn(gty))
	#	stop("Can't understand input.")
	
	if (!is.numeric(gty))
		stop("Input is not numeric.")
	
	by <- match.arg(by)
	if (by == "markers")
		return( rowMeans(gty, na.rm = na.rm)/2 )
	else
		return( colMeans(gty, na.rm = na.rm)/2 )
	
}

#' Calculate minor-allele frequency (MAFs) from a genotypes matrix
#' 
#' Just an alias or \code{freq(x, "markers")}.  Sticklers for nomenclature will be aware that this
#' function does not truly return the MAF unless the allele encoding is \code{"relative"}.
#' 
#' @seealso \code{\link{freq}}
#'
#' @export
maf <- function(...) freq(..., by = "markers")

#' Count number of missing genotypes by marker or sample
#'
#' @param gty a \code{genotypes} object
#' @param by margin over which to count missing calls
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
#' 
#' @return a named vector of missing-genotype *proportions*
#' 
#' @export
missingness <- function(gty, by = by, ...) {
	
	rez <- nmiss(gty, by = by, ...)
	return( rez/length(rez) )
	
}

#' Identify the consensus genotype call among a group of samples
#' 
#' @param gty a \code{genotypes} object
#' @param nas.allowed maximum proportion of missing genotypes to ignore 
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
	cons <- do.call("cbind", genoapply(gty, between, consensus, ..., strip = TRUE))
	
	#message("Iterating over groups...")
	seg <- segregating(cons, allow.het = FALSE, ...)
	
	return(seg)
	
}

#' Predict genotype of an F1 individual given genotypes of its parents
#' 
#' @param gty a \code{genotypes} object
#' @param na.rm logical; should missing genotypes be ignored?
#' 
#' @return a new \code{genotypes} object containing predicted genotypes for all possible F1s between
#' 	the parents in \code{gty}
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
#' @export
predict.f1 <- function(gty, na.rm = FALSE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	pairs <- combn(ncol(gty), 2)
	f1.names <- paste(colnames(gty)[ pairs[1,] ], colnames(gty)[ pairs[2,] ],
					  sep = "::")
	rez <- matrix(NA, ncol = ncol(pairs), nrow = nrow(gty),
				  dimnames = list(rownames(gty), f1.names))
	
	map <- attr(gty, "map")
	gty <- .copy.matrix.noattr(gty)
	message(paste("Predicting F1 genotypes for", ncol(pairs), "pairs of parents..."))
	for (i in seq_len(ncol(pairs))) {
		f1 <- as.vector(rowMeans(gty[ ,pairs[,i] ], na.rm = na.rm))
		het <- apply(gty, 1, function(x) any(x == 1, na.rm = na.rm))
		f1[het] <- NA
		rez[,i] <- f1
	}
	
	rez <- genotypes(rez, map = map, alleles = "01")
	return(rez)
	
}

#' Calculate the proportion of sites at which each sample is heterozygous
#' 
#' @param gty a \code{genotypes} object
#' @param na.rm logical; ignore missing genotypes when counting heterozygous sites?
#' 
#' @return a named vector with the proportion of sites of heterozygous for each sample
#' 
#' @export
heterozygosity <- function(gty, na.rm = TRUE, ...) {
	
	gty[ gty != 1 ] <- 0
	colMeans(gty, na.rm = na.rm)
	
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


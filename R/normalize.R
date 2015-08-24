## NB: functions in this file are adapted with modification from John Didion's CLASP package:
##	Didion JP et al. (2014) BMC Genomics. doi:10.1186/1471-2164-15-847.
##	https://github.com/jdidion/clasp/
## Source code from CLASP is provided under the CC-BY-SA-4.0 license.  As required by the terms
## of that license, please note that adaptation of CLASP code in this package does not necessarily
## constitute an endorsement by CLASP's author.

#' Perform quantile normalization of intensity data.
#'
#' @param gty a \code{genotypes} object
#' @param weights a vector of column weights
#' @param force re-run the normalization procedure even if \code{attr(,"normalized")} is set to \code{TRUE}
#'
#' @return A copy of the input object, with raw intensities replaced by the normalized ones.
#'
#' @details A simple wrapper around the quantile normalization originally described in Bolstad et al. (2003),
#' 	as implemented in \code{preprocessCore::normalize.quantiles.robust()}.
#'
#' @references
#' Bolstad BM et al. (2003) A comparison of normalization methods for high density oligonucleotide
#' 	array data based on bias and variance. Bioinformatics 19(2): 185-193.
#'
#' @export
quantile.normalize <- function(gty, weights = NULL, force = FALSE, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes' with intensity information attached.")
	
	if (is.null(attr(gty, "normalized")))
		attr(gty, "normalized") <- FALSE
	
	if (attr(gty, "normalized")) {
		if (!force) {
			message("Intensities seem to already be normalized")
			return(gty)
		}
		else {
			warning("Intensities seem to already be normalized; doing it again, but it won't help.")
		}
	}	
	
	if (!is.null(weights) && length(weights) != ncol(gty))
		stop("Dimensiosn of weight vector and intensity matrices don't match.")
	
	message(paste("Performing robust quantile normalization with", ifelse(is.null(weights), "no weights","weights"), "..."))
	x.norm <- preprocessCore::normalize.quantiles.robust( attr(gty, "intensity")$x, weights = weights)
	y.norm <- preprocessCore::normalize.quantiles.robust( attr(gty, "intensity")$y, weights = weights)
	
	colnames(x.norm) <- colnames(attr(gty, "intensity")$x)
	rownames(x.norm) <- rownames(attr(gty, "intensity")$x)
	colnames(y.norm) <- colnames(attr(gty, "intensity")$y)
	rownames(y.norm) <- rownames(attr(gty, "intensity")$y)
	
	attr(gty, "intensity") <- list(x = x.norm, y = y.norm)
	attr(gty, "normalized") <- TRUE
	return(gty)
	
}

## now supersded by reworked tQN(..., adjust.lrr = TRUE)
.calc.lrr <- function(gty, means, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes' with intensity information attached.")
	
	## keep only markers which are in reference set AND in the data
	m <- intersect(rownames(gty), names(means))
	
	message("Found ", length(m), " markers in intersection of reference values and this dataset.")
	message("Calculating raw LRR (log2 ratio of total intensitiy)...")
	
	## initialize LRR matrix
	lrr <- matrix(NA, nrow = nrow(gty), ncol = ncol(gty),
				  dimnames = dimnames(gty))
	
	## spoof the BAF matrix, unless it already exists
	baf <- matrix(-1, nrow = nrow(gty), ncol = ncol(gty),
				  dimnames = dimnames(gty))
	if (is.null(attr(gty, "baf")))
		attr(gty, "baf") <- baf
	
	## calculate LRRs in parallel over markers, but serially over samples
	## TODO: progress bar
	for (s in seq_len(ncol(gty))) {
		lrr[ m,s ] <- log2( with(attr(gty, "intensity"), x[m,s]+y[m,s])/means[m] )
	}
	attr(gty, "lrr") <- lrr
	
	message("Done.")
	
	return(gty)
	
}

#' Perform tQN normalization of intensity data.
#'
#' @param gty a \code{genotypes} object
#' @param thresholds thresholds for scaling of x- and y-intensities; defaults recommended in Staaf et al. (2008)
#' @param clusters a pre-computed matrix of cluster means
#' @param prenorm logical; if \code{TRUE}, perform quantile normalization on whole dataset before tQN procedure
#' @param xynorm logical; if \code{TRUE}, perform within-sample normalization of x- vs y-intensities
#' @param adjust.lrr logical; if \code{TRUE}, normalize post-tQN LRRs against population mean (\code{clusters$Rmean})
#'
#' @return A copy of the input object, with raw intensities replaced by the normalized ones.  Two additional attributes
#' \code{baf} and \code{lrr} store the BAF (B-allele frequency) and LRR (log2 intensity ratio).
#'
#' @details Implements thresholded quantile normalization (tQN) as described in Staaf et al. (2008).  Quantile
#' 	normalization as originally described in Bolstad et al. (2003) matches quantiles across multiple samples
#' 	so that all samples' intensities have the same empirical distribution.  The tQN instead matches the quantiles
#' 	of the x- and y-intensities sample-wise in order to reduce noise in the B-allele frequency (BAF) calculation
#' 	proposed by Peiffer et al. (2006).  The quantile-normalized intensities are then subject to a threshold to limit
#' 	on the ratio between the transformed and raw values.
#' 	NB: the quality of the result of tQN depends strongly on the reference clusters provided in \code{clusters},
#' 	so beware.
#' 	
#' 	The object \code{clusters} should be a dataframe with one row per marker and at least the following six columns:
#' 	\code{A.R, A.T}, the values of R and theta, respectively, for the centroid of the AA homozygous cluster;
#' 	\code{B.R, B.T}, likewise for the BB homozygous cluster; and \code{H.R, H.T}, likewise for the AB heterozygous
#' 	cluster.
#' 	
#' 	The transformations proposed by Peiffer et al. (2006) assume that most samples will fall into three well-defined
#' 	clusters at each marker, save for a relatively small proportion of abberrantly-hybridizing samples.  Indeed the BAF
#' 	is only well-defined in this case.  However, these assumptions are a bit too restrictive for many arrays, and
#' 	in particular for arrays which include copy-number probes.  It may be possible to obtain a tighter distribution
#' 	of LRR values by choosing \code{adjust.lrr = TRUE} and re-computing the LRR values against a reference distribution
#' 	independent of BAF or underlying clustering pattern.  For this, an additional column \code{Rmean} is required in
#' 	\code{clusters} which gives the "mean" (or other appropriately-chosen central value) of R across *all* possible
#' 	clusters at this marker.
#'
#' @references
#' Adapted from code provided by Johan Staaf to John Didion.
#' 
#' Staaf J et al. (2008) BMC Bioinformatics. doi:10.1186/1471-2105-9-409.
#' 
#' Bolstad BM et al. (2003) A comparison of normalization methods for high density oligonucleotide
#' 	array data based on bias and variance. Bioinformatics 19(2): 185-193.
#' 
#' Peiffer DA et al. (2006) Genome Res 16(9): 1136-1148. doi:10.1101/gr.5402306.
#'
#' Didion JP et al. (2014) BMC Genomics. doi:10.1186/1471-2164-15-847.
#'
#' @export
tQN <- function(gty, thresholds = c(1.5, 1.5), clusters = NULL,
				prenorm = TRUE, xynorm = TRUE, adjust.lrr = TRUE, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes' with intensity information attached.")
	
	if (is.null(clusters))
		stop("Must provide reference clusters.")
	
	markers <- rownames(gty)
	baf <- matrix(NA, ncol = ncol(gty), nrow = nrow(gty),
				  dimnames = list(rownames(gty), colnames(gty)))
	xnorm <- ynorm <- lrr <- baf
	intens.raw <- attr(gty, "intensity")
	intens.raw$x[ is.na(intens.raw$x) ] <- 0
	intens.raw$y[ is.na(intens.raw$y) ] <- 0
	
	if (prenorm) {
		message("Performing initial quantile normalization...")
		intens.norm <- list(x = preprocessCore::normalize.quantiles.robust(intens.raw$x),
							y = preprocessCore::normalize.quantiles.robust(intens.raw$y))
	}
	else {
		intens.norm <- intens.raw
	}
		
	message(paste("Performing tQN normalization", ifelse(xynorm, "with","without"),
				  "additional within-sample normalization."))
	if (interactive())
		pb <- txtProgressBar(min = 0, max = ncol(gty), style = 3)
	for (i in seq_len(ncol(gty))) {
		
		## add intermediate garbage-collect, in case of big objects (?)
		if ((i %% 10) && i > 0)
			gc()
		
		rez <- tQN.sample(markers, QN.thresholds = thresholds, clusters = clusters,
						  intens.raw$x[,i], intens.raw$y[,i],
						  intens.norm$x[,i], intens.norm$y[,i],
						  adjust.lrr = adjust.lrr,
						  xynorm = xynorm )
		
		baf[ rownames(rez),i ] <- rez$BAF
		lrr[ rownames(rez),i ] <- rez$LRR
		
		xnorm[ rownames(rez),i ] <- rez$X
		ynorm[ rownames(rez),i ] <- rez$Y
		if (interactive())
			setTxtProgressBar(pb, i)
	}
		
	message("Done.")
	attr(gty, "intensity") <- list(x = xnorm, y = ynorm)
	attr(gty, "normalized") <- TRUE
	attr(gty, "baf") <- baf
	attr(gty, "lrr") <- lrr
	return(gty)
	
}

## from JPD: <https://github.com/jdidion/sci/blob/master/projects/megamugaQC/R/normalize.R>
tQN.sample <- function(markers, X, Y, Xnorm, Ynorm, adjust.lrr = FALSE, mask = character(0),
					   QN.thresholds = c(1.5, 1.5), clusters, xynorm = FALSE, ...) {
	
	result <- data.frame(BAF = rep(NA, length(X)), LRR = rep(NA, length(X)), 
						 X = rep(NA, length(X)), Y = rep(NA, length(Y)),
						 row.names = markers, check.names = FALSE)
	
	if (all(is.na(X) | is.na(Y))) {
		# special handling of masked samples
		return(result)
	}
	
	data <- cbind(X = X, Y = Y)
	rownames(data) <- markers
	inorm <- cbind(X = Xnorm, Y = Ynorm)
	rownames(inorm) <- markers
	
	# Remove any rows that aren't in the cluster file
	i <- intersect(rownames(clusters), markers)
	ref <- clusters[ i, ]
	data <- data[ i, ]
	inorm <- inorm[ i, ]
	
	# normalize X vs Y for *this sample*
	if (xynorm)
		inorm <- preprocessCore::normalize.quantiles.robust(inorm)
	
	# Calculate R
	R <- inorm[,1] + inorm[,2]
	na <- is.na(data[,1]) | is.na(data[,2]) | data[,1] <= 0 | data[,2] <= 0
	
	# Thesholding of QN effect
	aff.x <- !na & (inorm[,1] / data[,1]) > QN.thresholds[1]
	inorm[ aff.x,1 ] <- QN.thresholds[1] * data[ aff.x,1 ]
	aff.y <- !na & (inorm[,2] / data[,2]) > QN.thresholds[2]
	inorm[ aff.y,2 ] <- QN.thresholds[2] * data[ aff.y,2 ]
	
	Th <- theta(inorm[,1], inorm[,2])
	Th[na] <- NA
	Th[!na] <- bound(Th[!na], 0 ,1)
	
	# Calculate tQN X and Y to fit theta and R
	Y <- R * tan(Th * pi/2) / (1 + tan(Th * pi/2))
	X <- R - Y
	
	# Calculate BAF and LRR from corrected R and T
	#baflrr <- baf.lrr(R, Th, ref)
	## call fast Rcpp version
	baflrr <- tQN_Cpp(R, Th, ref$A.T, ref$A.R, ref$B.T, ref$B.R, ref$H.T, ref$H.R)
	
	#result[rownames(data),] <- data.frame(baflrr, X = X, Y = Y)
	if (adjust.lrr) {
		baflrr$LRR <- log2( (X+Y)/ref$Rmean )
	}
	result[rownames(data),] <- data.frame(BAF = baflrr$BAF, LRR = baflrr$LRR, X = X, Y = Y)
	return(result)
	
}

## from JPD: <https://github.com/jdidion/sci/blob/master/projects/megamugaQC/R/normalize.R>
## NB: left here for reference, but has been superseded by much faster Rcpp version in tQN_Cpp()
baf.lrr <- function(R, Th, ref) {
	
	print("I am the walrus")
	BAF <- Th
	LRR <- R
	
	med.tAA <- median(ref$A.T, na.rm = TRUE)
	med.tBB <- median(ref$B.T, na.rm = TRUE)
	
	## loop on markers
	for (i in seq_along(Th)) {
		
		if (!is.num(Th[i])) {
			BAF[i] <- NaN
			LRR[i] <- NaN
			next
		}
		
		th <- Th[i]
		rr <- R[i]
		rAA <- ref$A.R[i]
		rAB <- ref$H.R[i]
		rBB <- ref$B.R[i]
		tAA <- ref$A.T[i]
		tAB <- ref$H.T[i]
		tBB <- ref$B.T[i]
		
		## check existence of thetas
		e.tAA <- is.num(tAA)
		e.tAB <- is.num(tAB)
		e.tBB <- is.num(tBB)
		
		## check for A/B allele swaps
		if (e.tAA & e.tBB & tAA > tBB) {
			r.tmp <- rBB
			t.tmp <- tBB
			rBB <- rAA
			tBB <- tAA
			rAA <- r.tmp
			tAA <- t.tmp
		}
		
		# 0: Test for inconsistencies between tAA/tAB/tBB and rAA/rAB/rBB
		if (((e.tAA & e.tAB) && tAA > tAB) ||
				((e.tAA & e.tBB) && tAA > tBB) ||
				((e.tAB & e.tBB) && tAB > tBB)) {
			BAF[i] <- LRR[i] <- NaN
		}
		# 1: Triple blank SNP
		else if (!(e.tAA|e.tAB|e.tBB)) {
			BAF[i] <- LRR[i] <- NaN
		}
		# 2: Blank for AB, AA, while positive for BB
		else if (!(e.tAA|e.tAB) & e.tBB) {
			if (th >= tBB) {
				BAF[i] <- 1
				LRR[i] <- ifelse(rBB[i] <= 0, NaN, log2(rr / rBB[i]))
			}
			else {
				BAF[i] <- LRR[i] <- NaN
			}
		}
		# 3: Blank for AB, BB, while positive for AA
		else if (e.tAA & !(e.tAB|e.tBB)) {
			if (th <= tAA) {
				BAF[i] <- 0
				LRR[i] <- ifelse(tAA[i] <= 0, NaN, log2(rr / tAA[i]))
			}
			else {
				BAF[i] <- LRR[i] <- NaN
			}
		}
		# 4: Blank for AB while positive for AA & BB
		else if (e.tAA & !e.tAB & e.tBB) {
			# no AB cluster exist for this SNP, while AA & BB exists. Set it to the closest of AA or BB
			min.index <- which.min(c(abs(tAA - th), abs(tBB - th)))
			if (min.index == 1 && th < tAA) {
				BAF[i] <- 0
				LRR[i] <- ifelse(tAA[i] <= 0, NaN, log2(rr / tAA[i]))
			}
			else if (min.index != 1 && th >= tBB) {
				BAF[i] <- 1
				LRR[i] <- ifelse(rBB[i] <= 0, NaN, log2(rr / rBB[i]))
			}
			else {
				BAF[i] <- LRR[i] <- NaN
			}
		}
		# 5: Blank for AA while positive for AB & BB
		else if (!e.tAA & e.tAB & e.tBB) {
			if (th >= tBB) {
				BAF[i] <- 1
				LRR[i] <- NaN
			}
			# 5.1: SNP is "correctly between" ref$AB_T_Mean and ref$BB_T_Mean
			else if (th >= tAB) {
				# interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
				BAF[i] <- 0.5 + 0.5 * (th - tAB) / (tBB - tAB)
				eR <- rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB))
				LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
			}
			# 5.2: Heterozygous SNP is subjected to deletion or UPD of allele B making it unexectedly to be 
			# between ref$AA_T_Mean and ref$AB_T_Mean where it normally should not NOT BE.
			else {
				BAF[i] <- ifelse(th < med.tAA, 0, 0.5 * (th - med.tAA) / (tAB - med.tAA))
				LRR[i] <- NaN                   
			}
		}
		# 6: Blank for BB while positive for AA & AB
		else if (e.tAA & e.tAB & !e.tBB) {
			if (th < tAA) {
				BAF[i] <- 0
				LRR[i] <- NaN
			}
			# 6.1: SNP is "correctly between" ref$AA_T_Mean and ref$AB_T_Mean
			else if (th <= tAB) {
				#interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
				BAF[i] <- 0.5* (th - tAA) / (tAB - tAA)
				eR <- rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA))
				LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
			}
			# 2: Heterozygous SNP is subjected to deletion or UPD of allele A making it unexectedly to be 
			# between ref$AB_T_Mean and ref$BB_T_Mean where it normally should not NOT BE.
			else {
				BAF[i] <- ifelse(th > med.tBB, 1, 0.5 + 0.5 * (th - tAB) / (med.tBB[i] - tAB))
				LRR[i] <- NaN
			}
		}
		# 7: positive for AA & BB & AB, Illumina style calculation
		else {
			if (th < tAB) {
				BAF[i] <- ifelse(th < tAA, 0, 0.5 * (th - tAA) / (tAB - tAA))
				eR <- rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA))
			}
			else {
				BAF[i] <- ifelse(th >= tBB, 1, 0.5 + 0.5 * (th - tAB) / (tBB - tAB))
				eR <- rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB))
			}
			LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
		}
	}
	
	return( cbind(BAF = bound(BAF, 0, 1), LRR = LRR) )
}

## check that value "is a number" in the sense of being numeric, real and finite
is.num <- function(x) !(is.na(x) | is.nan(x) | is.infinite(x) | is.complex(x))

## calculate allelic intensity ratio
theta <- function(x, y) 2/pi*atan(x/y)

## clip values to fall within given bounds
bound <- function(x, lo, hi) {
	n <- is.num(x)
	x[ n & x < lo ] <- lo
	x[ n & x > hi ] <- hi
	return(x)
}
## qc.R
## functions for sample-level QC

## helper function for computing quantiles of columns of a matrix
column.quantiles <- function(x, q = seq(0,1,0.1), ..., .progress = "none") {
	
	.col.quant <- function(y) {
		data.frame(q = q, value = quantile(y, q, na.rm = TRUE))
	}
	
	rez <- plyr::adply(x, 2, .col.quant, .progress = .progress)
	colnames(rez)[1] <- "iid"
	rez$iid <- reorder(factor(rez$iid), rez$value, median)
	return(rez)
	
}

#' Summarize hybridization intensity by sample
#' 
#' @param gty a \code{genotypes} object with intensity data attached
#' @param q a vector of quantiles (in [0,1])
#' @param .progress show a progress bar; passed through to \code{plyr}, see \code{\link{plyr::ddply}}
#' 
#' @return a dataframe containing intensity quantiles for each sample, merged with
#' 	sample metadata (if present)
#' 	
#' @details The \code{q}th quantiles of "sum-intensity" are computed for each sample, across all markers
#' 	in the input dataset.  Missing values are silently ignored.  We define "sum-intensity" as
#' 	sum(sqrt(x_i^2 + y_i^2)), rather than sum(x_i + y_i), based on the intuition that distance from
#' 	the origin represents total hybridization intensity after application of Illumina's proprietary
#' 	affine-transformation scheme to the raw fluorescences.
#'
#' @seealso \code{\link{summarize.calls}}, \code{\link{intensity.vs.ref}}, \code{\link{run.qc.checks}}
#' 	 	
#' @export
summarize.intensity <- function(gty, q = seq(0,1,0.1), ..., .progress = "none") {
	
	if (!inherits(gty, "genotypes") && .has.valid.intensity(gty))
		stop("Please supply an object of class 'genotypes' with valid intensity information.")
	
	si <- with(attr(gty, "intensity"), sqrt(x^2 + y^2))
	rez <- column.quantiles(si, q = q, .progress = .progress)
	
	if (.has.valid.ped(gty))
		rez <- merge(rez, attr(gty, "ped"))
	
	return(rez)
	
}

#' Summarize genotype calls by sample or marker
#' 
#' @param gty a \code{genotypes} object with intensity data attached
#' @param by get call rates by sample or by marker
#' 
#' @return a dataframe with counts of A (reference or major allele), B (alternate or minor allele),
#' 	H (heterozygous) and N (missing/no-call) by either sample or marker.
#' 	
#' @details Any metadata associated with markers or samples is merged into the final result.
#'
#' @seealso \code{\link{summarize.intensity}}, \code{\link{intensity.vs.ref}}, \code{\link{run.qc.checks}}
#' 	 	
#' @export
summarize.calls <- function(gty, by = c("samples","markers"), ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	by <- match.arg(by)
	if (by == "samples")
		which.dim <- 2
	else
		which.dim <- 1
	
	## convert genotypes to numeric for faster summaries
	gty <- .copy.matrix.noattr(recode.genotypes(gty, "01"))
	## convert NAs to special code
	gty[ is.na(gty) ] <- 3
	## fast summary with tabulate()
	counts <- t(apply(gty, which.dim, function(x) tabulate(x+1, nbins = 4)))
	
	if (by == "samples")
		return(data.frame(iid = colnames(gty),
						  A = counts[,1], B = counts[,3], H = counts[,2], N = counts[,4]))
	else
		return(data.frame(marker = rownames(gty),
						  A = counts[,1], B = counts[,3], H = counts[,2], N = counts[,4]))
	
}

#' KS-test for difference in intensity distributions
#' 
#' @param gty a \code{genotypes} object with intensity data attached
#' @param ref a matrix, vector or object coercible to such, containing sum-intensities from reference samples
#' 
#' @return a named vector of D_j, the Kolmogorov-Smirnov test statistic for each sample j
#' 	
#' @details This function detects potential failed arrays by performing the Kolmogorov-Smirnov test
#' 	for difference between the "sum-intensities" (see \code{\link{summarize.intensity}}) of each sample
#' 	and some reference distribution of "sum-intensities" of known good samples.  This test should (obviously)
#' 	be performed *before* any normalizations are applied.  As such it may be useful for detecting batch
#' 	effects, although that possibility has not been systematically explored.
#' 	
#' 	The distribution of "sum-intensity" across an array is expected to be approximately normal. Outliers for
#' 	the D statistic come in two flavours (cf. Didion et al. (2014)): samples which fail completely, having a
#' 	heavily right-skewed intensity distribution; and samples which are genetically diverged from the reference
#' 	sample/species used in array design.  Divergent samples have a spike of intensities near zero, representing
#' 	failed hybridization due to off-target variation, but an otherwise normal intensity distribution. 
#'
#' @references
#' Didion JP et al. (2014) SNP array profiling of mouse cell lines identifies their strains of origin
#' 	and reveals cross-contamination and widespread aneuploidy. BMC Genomics 15(1): 847. doi:10.1186/1471-2164-15-847.
#'
#' @seealso \code{\link{summarize.intensity}}, \code{\link{summarize.calls}}, \code{\link{run.qc.checks}}
#' 	 	
#' @export
intensity.vs.ref <- function(gty, ref, ...) {
	
	if (!inherits(gty, "genotypes") && .has.valid.intensity(gty))
		stop("Please supply an object of class 'genotypes' with valid intensity information.")
	
	.colwise.ks <- function(x, y) {
		ks.test(x, y, alternative = "less")$statistic
	}
	
	si <- with(attr(gty, "intensity"), sqrt(x^2+y^2))
	apply(si, 2, .colwise.ks, y = as.vector(ref))
	
}

#' Perform basic sample-wise QC on genotype calls and intensities
#'
#' @param gty a \code{genotypes} object
#' @param ref.intensity a vector of "sum-intensities" from known good reference samples
#' @param max.H threshold for count of heterozygous calls, above which a sample is flagged
#' @param max.N threshold for count of no-calls, above which a sample is flagged
#' @param max.D upper threshold for D-statistic (see \code{\link{intensity.vs.ref}}) above which a
#' 	sample is flagged
#' @param min.D lower threshold for D-statistic (see \code{\link{intensity.vs.ref}}) above which a
#' 	sample is flagged
#' 
#' @return a copy of the input with sample filters set, and an object of class \code{QC.result} in
#' 	attr(,"qc")
#' 	
#' @details A wrapper for the sample-level QC functions, applied to genotype calls (always) and
#' 	hybridization intensities (if present.)  Samplies which fail are flagged but not actually
#' 	dropped from the result.
#' 	
#' @seealso \code{\link{summarize.calls}}, \code{\link{intensity.vs.ref}}, \code{\link{apply.filters}}
#' 
#' @export
run.qc.checks <- function(gty, ref.intensity = NULL,
						  max.H = Inf, max.N = Inf, max.D = Inf, min.D = Inf, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	sites <- rep(FALSE, nrow(gty))
	samples <- rep(FALSE, ncol(gty))
	if (!is.null(attr(gty, "filter.sites")))
		sites <- attr(gty, "filter.sites")
	if (!is.null(attr(gty, "filter.samples")))
		samples <- attr(gty, "filter.samples")
	
	if (length(sites) != nrow(gty))
		warning("Site filters don't match dimensions of genotype matrix.")
	if (length(samples) != ncol(gty))
		warning("Sample filters don't match dimensions of genotype matrix.")
	
	qc.rez <- list()
	message("Performing QC checks on genotype calls...")
	qc.rez$calls <- summarize.calls(gty, "samples")
	samples <- samples | with(qc.rez$calls, H > max.H | N > max.N)
	
	if (.has.valid.intensity(gty)) {
		message("Performing QC checks on hybridization intensities...")
		qc.rez$intensity <- summarize.intensity(gty, q = seq(0.05, 0.95, 0.05))
		if (!is.null(ref.intensity)) {
			qc.rez$D <- intensity.vs.ref(gty, ref.intensity)
			samples <- samples | (qc.rez$D > max.D | qc.rez$D < min.D)
		}
	}
	
	message(paste(sum(sites),"markers and", sum(samples), "samples now flagged as low-quality."))
	class(qc.rez) <- c("QC.result", class(qc.rez))
	
	attr(gty, "qc") <- qc.rez
	attr(gty, "filter.sites") <- sites
	attr(gty, "filter.samples") <- samples
	
	return(gty)
	
}

#' Drop samples and/or markers flagged as low-quality
#'
#' @param gty a \code{genotypes} object
#' @param apply.to dimensions along which to apply filters (samples, sites or both)
#' 
#' @return a copy of the input with flagged markers and/or samples dropped
#' 	
#' @details In the output object, all sample and marker filters are set to \code{FALSE}.
#' 	
#' @seealso \code{\link{run.qc.checks}}, \code{\link{apply.filters}}
#' 
#' @export
apply.filters <- function(gty, apply.to = c("both","samples","markers"), ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	do.samples <- FALSE
	do.sites <- FALSE
	apply.to <- match.arg(apply.to)
	if (apply.to == "both") {
		do.samples <- TRUE
		do.sites <- TRUE
	}
	else if (apply.to == "samples")
		do.samples <- TRUE
	else
		do.sites <- TRUE
	
	if (!any(do.sites, do.samples)) {
		message("Nothing to do; this code should never be reached.")
		return(gty)
	}
	
	sites <- FALSE
	samples <- FALSE
	if (do.sites) {
		if (!is.null(attr(gty, "filter.sites"))) {
			sites <- attr(gty, "filter.sites")
			if (length(sites) != nrow(gty))
				warning("Site filters don't match dimensions of genotype matrix.")
		}
	}
	if (do.samples) {
		if (!is.null(attr(gty, "filter.samples"))) {
			samples <- attr(gty, "filter.samples")
			if (length(samples) != ncol(gty))
				warning("Sample filters don't match dimensions of genotype matrix.")
		}
	}
	
	message(paste("Dropping", sum(sites), "markers and", sum(samples), "samples..."))
	return( gty[ !sites,!samples ] )
	
}

## just return the filter flags
get.filters <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	sites <- rep(FALSE, nrow(gty))
	samples <- rep(FALSE, ncol(gty))
	if (!is.null(attr(gty, "filter.samples")))
		samples <- attr(gty, "filter.samples")
	if (!is.null(attr(gty, "filter.sites")))
		sites <- attr(gty, "filter.sites")
	
	names(samples) <- colnames(gty)
	names(sites) <- rownames(gty)
	
	return(list(sites = sites, samples = samples))
	
}

#' @export
summarize.filters <- function(gty, ...) {
	
	filters <- get.filters(gty)
	sapply(filters, sum, na.rm = TRUE)
	
}
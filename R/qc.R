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

## helper function to compute sum-intensity
## NB: no sanity checks here!!
.si <- function(gty, ...) {
	with(attr(gty, "intensity"), sqrt(x^2 + y^2))
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
#' 	sqrt(x_i^2 + y_i^2), rather than sum(x_i + y_i), based on the intuition that distance from
#' 	the origin represents total hybridization intensity after application of Illumina's proprietary
#' 	affine-transformation scheme to the raw fluorescences.
#'
#' @seealso \code{\link{summarize.calls}}, \code{\link{intensity.vs.ref}}, \code{\link{run.sample.qc}}
#' 	 	
#' @export
summarize.intensity <- function(gty, q = seq(0,1,0.1), ..., .progress = "none") {
	
	if (!inherits(gty, "genotypes") && .has.valid.intensity(gty))
		stop("Please supply an object of class 'genotypes' with valid intensity information.")
	
	si <- .si(gty)
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
#' @details Any metadata associated with markers or samples is *not* merged into the final result, in order
#' 	to preserve parallelism between rows of the result and of the parent object.
#'
#' @seealso \code{\link{summarize.intensity}}, \code{\link{intensity.vs.ref}}, \code{\link{run.sample.qc}}
#' 	 	
#' @export
summarize.calls <- function(gty, by = c("samples","markers"), counts = TRUE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	by <- match.arg(by)
	if (by == "samples")
		which.dim <- 2
	else
		which.dim <- 1
	
	if (counts)
		denom <- 1
	else
		denom <- dim(gty)[ setdiff(c(2,1), which.dim) ]
	
	## convert genotypes to numeric for faster summaries
	if (!(attr(gty, "alleles") %in% c("01","relative")))
		gty <- .copy.matrix.noattr(recode.genotypes(gty, "01"))
	## convert NAs to special code
	gty[ is.na(gty) ] <- 3
	## fast summary with tabulate()
	counts <- t(apply(gty, which.dim, function(x) tabulate(x+1, nbins = 4)))
	
	if (by == "samples")
		return(data.frame(iid = colnames(gty),
						  A = counts[,1]/denom, B = counts[,3]/denom,
						  H = counts[,2]/denom, N = counts[,4]/denom))
	else
		return(data.frame(marker = rownames(gty),
						  A = counts[,1]/denom, B = counts[,3]/denom,
						  H = counts[,2]/denom, N = counts[,4]/denom))
	
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
#' @seealso \code{\link{summarize.intensity}}, \code{\link{summarize.calls}}, \code{\link{run.sample.qc}}
#' 	 	
#' @export
intensity.vs.ref <- function(gty, ref, ...) {
	
	if (!inherits(gty, "genotypes") && .has.valid.intensity(gty))
		stop("Please supply an object of class 'genotypes' with valid intensity information.")
	
	.colwise.ks <- function(x, y) {
		ks.test(x, y, alternative = "less")$statistic
	}
	
	si <- .si(gty)
	apply(si, 2, .colwise.ks, y = as.vector(ref))
	
}

#' Perform basic sample-wise QC on genotype calls and intensities
#'
#' @param gty a \code{genotypes} object
#' @param ref.intensity a vector of "sum-intensities" from known good reference samples
#' @param max.H threshold for count of heterozygous calls, above which a sample is flagged;
#' 	OR a named list of thresholds, with names to match "family ID" (column \code{fid}) in pedigree
#' @param max.N threshold for count of no-calls, above which a sample is flagged; OR a named list
#' 	as above
#' @param max.D upper threshold for D-statistic (see \code{\link{intensity.vs.ref}}) above which a
#' 	sample is flagged; OR a named list as above
#' @param min.D lower threshold for D-statistic (see \code{\link{intensity.vs.ref}}) above which a
#' 	sample is flagged; OR a named list as above
#' @param apply logical; if \code{TRUE}, remove samples failing the filters, rather than flagging them
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
run.sample.qc <- function(gty, ref.intensity = NULL,
						  max.H = Inf, max.N = Inf, max.D = Inf, min.D = Inf, 
						  apply = FALSE, hits = 0, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	fl <- get.filters(gty)
	
	if (length(fl$sites) != nrow(gty))
		warning("Site filters don't match dimensions of genotype matrix.")
	if (length(fl$samples) != ncol(gty))
		warning("Sample filters don't match dimensions of genotype matrix.")
	
	qc.rez <- list()
	message("Performing QC checks on genotype calls...")
	qc.rez$calls <- summarize.calls(gty, "samples")
	
	sm <- as.character(qc.rez$calls$iid)
	fid <- as.character(attr(gty, "ped")[ sm, "fid" ])
	if (is.list(max.H))
		max.H <- as.vector(max.H[fid])
	if (is.list(max.N))
		max.N <- as.vector(max.N[fid])
	if (is.list(max.D))
		max.D <- as.vector(max.D[fid])
	if (is.list(min.D))
		min.D <- as.vector(max.D[fid])
	
	fail.hets <- as.vector(qc.rez$calls$H > max.H)
	fail.ns <- as.vector(qc.rez$calls$N > max.N)
	fail.ds <- rep(FALSE, ncol(gty))
	
	if (.has.valid.intensity(gty)) {
		message("Performing QC checks on hybridization intensities...")
		qc.rez$intensity <- summarize.intensity(gty, q = seq(0.05, 0.95, 0.05))
		if (!is.null(ref.intensity)) {
			qc.rez$D <- intensity.vs.ref(gty, ref.intensity)
			fail.ds <- as.vector(qc.rez$D > max.D | qc.rez$D < min.D)
		}
	}
	
	fl$samples[fail.hets] <- paste0(fl$samples[fail.hets], "H")
	fl$samples[fail.ns] <- paste0(fl$samples[fail.ns], "N")
	fl$samples[fail.ds] <- paste0(fl$samples[fail.ds], "I")
	isfl <- lapply(fl, function(x) {
		x[ is.na(x) ] <- ""
		nchar(x) > hits
	})
	qc.rez$calls$filter <- isfl$samples
	
	message(paste(sum(isfl$sites),"markers and", sum(isfl$samples), "samples now flagged as low-quality."))
	class(qc.rez) <- c("QC.result", class(qc.rez))
	
	attr(gty, "qc") <- qc.rez
	attr(gty, "filter.sites") <- fl$sites
	attr(gty, "filter.samples") <- fl$samples
	
	if (apply)
		gty <- apply.filters(gty)
	
	return(gty)
	
}

#' Perform basic marker-wise QC on genotype calls
#'
#' @param gty a \code{genotypes} object
#' @param max.H threshold for count of heterozygous calls, above which a site is flagged
#' @param max.N threshold for count of no-calls, above which a site is flagged
#' @param min.hom threshold for count of homozygous calls, at or below which a site is flagged
#' @param apply logical; if \code{TRUE}, remove samples failing the filters, rather than flagging them
#' 
#' @return a copy of the input with sample filters set, and an object of class \code{marker.QC.result} in
#' 	attr(,"marker.qc")
#' 	
#' @details A wrapper for the marker-level QC functions, applied to genotype calls.  Samplies which
#' 	fail are flagged but not actually dropped from the result (unless \code{apply = TRUE}.)
#' 	
#' @seealso \code{\link{summarize.calls}}, \code{\link{run.sample.qc}}, \code{\link{apply.filters}}
#' 
#' @export
run.marker.qc <- function(gty, max.H = Inf, max.N = Inf, min.hom = 0, hits = 0, apply = FALSE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	fl <- get.filters(gty)
	
	if (length(fl$sites) != nrow(gty))
		warning("Site filters don't match dimensions of genotype matrix.")
	if (length(fl$samples) != ncol(gty))
		warning("Sample filters don't match dimensions of genotype matrix.")
	
	qc.rez <- list()
	message("Performing QC checks on genotype calls per marker...")
	qc.rez$calls <- summarize.calls(gty, "markers")
	
	fail.hets <- as.vector(qc.rez$calls$H > max.H)
	fail.ns <- as.vector(qc.rez$calls$N > max.N)
	fail.freq <- as.vector(with(qc.rez$calls, A+B <= min.hom))
	
	fl$sites[fail.hets] <- paste0(fl$sites[fail.hets], "H")
	fl$sites[fail.ns] <- paste0(fl$sites[fail.ns], "N")
	fl$sites[fail.freq] <- paste0(fl$sites[fail.ns], "F")
	isfl <- lapply(fl, function(x) {
		x[ is.na(x) ] <- ""
		nchar(x) > hits
	})
	
	qc.rez$calls$filter <- isfl$sites
	
	message(paste(sum(isfl$sites),"markers and", sum(isfl$samples), "samples now flagged as low-quality."))
	class(qc.rez) <- c("marker.QC.result", class(qc.rez))
	
	attr(gty, "marker.qc") <- qc.rez
	attr(gty, "filter.sites") <- fl$sites
	attr(gty, "filter.samples") <- fl$samples
	
	if (apply)
		gty <- apply.filters(gty)
	
	return(gty)
	
}

## given a list of filters, update it with 'data' a list(filter_name = logical(to_filter))
.update.filters <- function(fl, what = c("samples","markers"), data, ...) {
	
	for (d in names(data)) {
		ii <- which(data[[d]])
		for (i in ii) {
			these <- fl[[what]][[i]]
			these <- unique(c(these, d))
			fl[[what]][[i]] <- these
		}
	}
	
	return(fl)
	
}

## initialize a filter list from vector of names
.init.filters <- function(nn,...) {
	
	setNames( rep("", length(nn)), nn )
	
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
#' @seealso \code{\link{run.sample.qc}}, \code{\link{apply.filters}}
#' 
#' @export
apply.filters <- function(gty, apply.to = c("both","samples","markers"), hits = 0, ...) {
	
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
	fl <- is.filtered(gty, hits = hits)
	
	if (do.sites)
		sites <- fl$sites
		if (length(sites) != nrow(gty))
			warning("Site filters don't match dimensions of genotype matrix.")
	
	if (do.samples)
		samples <- fl$samples
		if (length(samples) != ncol(gty))
			warning("Sample filters don't match dimensions of genotype matrix.")
	
	message(paste("Dropping", sum(sites), "markers and", sum(samples), "samples..."))
	return( gty[ !sites,!samples ] )
	
}


#' Check if markers or samples are marked with filters
#' 
#' @param gty a \code{genotypes} object
#' @param hits integer; maximum number of filters which can be set before a marker or sample is flagged
#' 
#' @return a list with two elements: \code{$sites}, logical vector of filter status for markers;
#' 	and \code{$samples}, logical vector of filter status for samples
#' 
#' @export
is.filtered <- function(gty, hits = 0, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	fl <- get.filters(gty)
	samples <- nchar(fl$samples) > hits
	sites <- nchar(fl$sites) > hits
	
	names(samples) <- colnames(gty)
	names(sites) <- rownames(gty)
	
	return(list(sites = sites, samples = samples))
	
}

## return the filter lists
get.filters <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	fsites <- attr(gty, "filter.sites")
	if (is.null(fsites))
		fsites <- setNames( rep("", nrow(gty)), rownames(gty) )
		
	fsamples <- attr(gty, "filter.samples")
	if (is.null(fsamples))
		fsamples <- setNames( rep("", ncol(gty)), colnames(gty) )
	
	fsamples[ is.na(fsamples) ] <- ""
	fsites[ is.na(fsites) ] <- ""
	
	return(list(sites = fsites, samples = fsamples))
	
}

#' @export
summarize.filters <- function(gty, filter.codes = "NHIF", ...) {
	
	fl <- get.filters(gty)
	codes <- unlist(strsplit(filter.codes, ""))
	sapply(fl, function(x) sapply(codes, function(y) sum(grepl(y, x))))
	
}
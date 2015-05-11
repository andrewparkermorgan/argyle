## qc.R
## functions for sample-level QC

## helper function for computing quantiles of columns of a matrix
column.quantiles <- function(x, q = seq(0,1,0.1), ..., .progress = "none") {
	
	.col.quant <- function(y) {
		data.frame(q = q, value = quantile(y, q, na.rm = TRUE))
	}
	
	rez <- adply(x, 2, .col.quant, .progress = .progress)
	colnames(rez)[1] <- "id"
	rez$id <- reorder(factor(rez$id), rez$value, median)
	return(rez)
	
}

## compute quantiles of sum-intensity by sample and report as dataframe
summarize.intensity <- function(gty, q = seq(0,1,0.1), ..., .progress = "none") {
	
	if (!inherits(gty, "genotypes") && .has.valid.intensity(gty))
		stop("Please supply an object of class 'genotypes' with valid intensity information.")
	
	si <- with(attr(gty, "intensity"), sqrt(x^2 + y^2))
	rez <- column.quantiles(si, q = q, .progress = .progress)
	
	return(rez)
	
}

## count Hs, Ns by sample or by marker
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
		return(data.frame(id = colnames(gty),
						  A = counts[,1], B = counts[,3], H = counts[,2], N = counts[,4]))
	else
		return(data.frame(marker = rownames(gty),
						  A = counts[,1], B = counts[,3], H = counts[,2], N = counts[,4]))
	
}

## apply Kolmogorov-Smirnov test for equality of distribution between sum-intensity of each sample
## and some appropriate reference (NB: choice of reference is important)
intensity.vs.ref <- function(gty, ref, ...) {
	
	if (!inherits(gty, "genotypes") && .has.valid.intensity(gty))
		stop("Please supply an object of class 'genotypes' with valid intensity information.")
	
	.colwise.ks <- function(x, y) {
		ks.test(x, y, alternative = "less")$statistic
	}
	
	si <- with(attr(gty, "intensity"), sqrt(x^2+y^2))
	apply(si, 2, .colwise.ks, y = ref)
	
}

## wrapper function for call- and intensity-level QC, done sample-wise
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

## apply filters to input object, dropping samples and/or sites flagged as low-quality
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

## count how many sites and samples are flagged for removal
summarize.filters <- function(gty, ...) {
	
	filters <- get.filters(gty)
	sapply(filters, sum, na.rm = TRUE)
	
}
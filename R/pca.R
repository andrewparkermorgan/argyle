## pca.R
## functions for performing PCA on genotypes and hybridization intensities

## impute allele freq in place of missing values
.fudge.missing.genotypes <- function(g, check.variance = TRUE, eps = 1e-6, ...) {
	
	X <- t(.copy.matrix.noattr(g))
	## replace missing values with column mean
	if (any(is.na(colMeans(X)))) {
		message("\treplacing missing values with minor-allele frequency...")
		X <- apply(X, 2, function(x) {
			x[ is.na(x) ] <- mean(x, na.rm = TRUE)
			return(x)
		})
	}
	
	## remove zero-variance columns
	if (check.variance) {
		const <- abs(colVars(X) - 0) <= eps
		const[ is.na(const) ] <- TRUE
		X <- X[ ,!const ]
	}
	
	return(X)
	
}

## internal function for actually doing PCA
.do.pca.genotypes <- function(gty, extras = NULL, what = c("genotypes","intensity","norm"),
							  fast = FALSE, scale = TRUE, center = TRUE, eps = 1e-6, K = 1:3, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (!is.null(extras) && !inherits(extras, "genotypes"))
		stop(paste("If supplying extra samples to project onto pre-computed PCs",
				   "they should also be suppled as class 'genotypes'."))
	
	what <- match.arg(what)
	.make.pca.input <- function(g) {
		if (what == "genotypes") {
			if (!is.numeric(g)) {
				g <- recode.genotypes(g, "01")
			}
			X <- t(.copy.matrix.noattr(g))
			## replace missing values with column mean
			if (any(is.na(colMeans(X)))) {
				message("\treplacing missing values with minor-allele frequency...")
				X <- apply(X, 2, function(x) {
					x[ is.na(x) ] <- mean(x, na.rm = TRUE)
					return(x)
				})
			}
		}
		else if (what == "intensity") {
			if (.has.valid.intensity(g)) {
				x <- with(attr(g, "intensity"), .copy.matrix.noattr(x))
				y <- with(attr(g, "intensity"), .copy.matrix.noattr(y))
				rownames(x) <- paste0(rownames(x), "_x")
				rownames(y) <- paste0(rownames(y), "_y")
				x[ is.na(x) ] <- 0
				y[ is.na(y) ] <- 0
				X <- t(rbind(x,y))	
			}
			else {
				stop("Need valid intensity matrices in order to do PCA on intensity.")
			}
		}
		else if (what == "norm") {
			if (.has.valid.baflrr(g)) {
				x <- attr(g, "baf")
				#y <- attr(g, "lrr")
				rownames(x) <- paste0(rownames(x), "_x")
				#rownames(y) <- paste0(rownames(y), "_y")
				x[ is.na(x) ] <- 0
				#y[ is.na(y) ] <- 0
				#X <- t(rbind(x,y))	
				X <- t(x)
			}
			else {
				stop("Need valid BAF and LRR matrices in order to do PCA on normalized intensity.")
			}
		}
		return(X)
	}
	
	## remove zero-variance columns, if scaling requested
	.check.variance <- function(x) {
		const <- abs(colVars(x) - 0) <= eps
		const[ is.na(const) ] <- TRUE
		return(!const)
	}

	## finally do PCA
	message("Preparing input matrices...")
	
	if (!is.null(extras)) {
		suppressMessages(both <- merge(gty, extras))
		both <- .make.pca.input(both)
		keep <- .check.variance(both[ colnames(gty), ]) & .check.variance(both[ colnames(extras), ])
		features <- both[ colnames(gty),keep ]
		topredict <- both[ colnames(extras),keep ]
	}
	else {
		features <- .make.pca.input(gty)
		keep <- .check.variance(features)
		features <- features[ ,keep ]
	}
		
	message(paste("Computing principal components of", what, "matrix (", paste(dim(features), collapse = " x ") ,")..."))
	if (fast) {
		message("\t(using corpcor::fast.svd() ...)")
		pc <- .fast.pca(features, K = K, .center = center, .scale = scale)	
	}
	else {
		message("\t(using base::prcomp() ...)")
		pc <- prcomp(features, center = center, scale. = scale)
	}
	proj <- predict(pc)
	if (!is.null(extras)) {
		message("Projecting extra samples onto existing PCs...")
		predicted <- scale(topredict, pc$center, pc$scale) %*% pc$rotation
		#predicted <- predict(pc, newdata = topredict)
		proj <- rbind(proj, predicted)
	}
	
	## add %variance explained
	attr(proj, "effective.dim") <- dim(features)
	attr(proj, "pca") <- pc
	attr(proj, "explained") <- pc$sdev^2/sum(pc$sdev^2)
	
	message("Done.")
	return(proj)
	
}

## do (faster?) PCA via SVD using corpcor::fast.svd()
.fast.pca <- function(x, K = 1:3, .scale = TRUE, .center = TRUE, ...) {
	
	xx <- scale(x, scale = .scale, center = .center)
	y <- xx/(ncol(xx)-1)
	udv <- corpcor::fast.svd(y)
	
	proj <- udv$u
	rownames(proj) <- rownames(x)
	colnames(proj) <- paste0("PC", 1:ncol(proj)) 
	rot <- udv$v
	rownames(rot) <- rownames(rot)
	colnames(rot) <- paste0("PC", 1:ncol(rot)) 
	
	rez <- list(sdev = udv$d,
				rotation = rot,
				x = proj,
				center = .center, scale = .scale)
	class(rez) <- c("prcomp", class(rez))
	return(rez)
	
}

#' Perform PCA on a \code{genotypes} object
#'
#' Performs principal components analysis (PCA) on either (numerically-coded) genotypes or
#' on the underliny 2D hybridization-intensity matrices.
#'
#' @param x an object of class \code{genotypes}
#' @param extras a second dataset of class \code{genotypes}, to be projected onto PCs computed
#' 	from the first. NB: this function *does not* verify that the underlying features (ie. markers) match.
#' @param what \code{"genotypes"} to do PCA on genotypes (coded 0/1/2); \code{"intensity"} to do PCA on
#' 	underlying 2D intensities.  The latter triggers an error if intensity matrices are absent.
#' @param K how many PCs to return.
#' @param fast if \code{TRUE}, use \code{corpcor::fast.svd()} to speed up calculations
#' @param ... other parameters for call to \code{prcomp}, such as \code{center} and \code{scale}
#' 	(both \code{TRUE} by default)
#'
#' @details Uses base-\code{R}'s \code{prcomp} under the hood (unless \code{fast = TRUE}). By default,
#' 	columns are centered and scaled; not doing so will likely produce a strange result.  When doing PCA on genotypes,
#' 	missing values are replaced with the column mean, which in many circumstances can be interpreted as the
#' 	minor-allele frequency. (This is very similar to the behavior of PLINK.)  When doing PCA on intensities, missing
#' 	values are set to zero -- but even no-call genotypes have nonmissing intensities on Illumina arrays, so this is unlikely
#' 	to have any effect in practice.
#'
#' @return a dataframe with as many rows as samples, in which the first columns are sample IDs and any associated
#' 	metadata (as returned by \code{samples(x)}), followed by the first K PCs.  Scaled eigenvalues (ie. percent
#' 	of variance explained) are provided as \code{attr(,"explained")}.
#'
#' @seealso \code{pca.plink()} for using PLINK's (much faster and more powerful) implementation
#'
#' @export
pca.genotypes <- function(x, extras = NULL, what = c("genotypes","intensity","norm"), K = 3, fast = FALSE, ...) {
	
	if (!inherits(x, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (!is.numeric(K) || K < 1)
		stop("You probably don't want to do PCA with <1 dimension.")
	
	## run the PCA
	rez <- .do.pca.genotypes(x, extras = extras, what = what, fast = fast, K = K, ...)
	
	## add sample metatdata, if any
	meta <- data.frame(iid = rownames(rez))
	rownames(meta) <- rownames(rez)
	if (.has.valid.ped(x)) {
		fam <- attr(x, "ped")
		if (!is.null(extras))
			fam <- rbind(fam, attr(extras, "ped"))
		meta <- merge(meta, fam, by = "iid")
		rownames(meta) <- as.character(meta$iid)
		if (!(nrow(meta) == nrow(rez) && all(rownames(rez) %in% meta$iid)))
			stop("Sample metadata and PCA result are out of sync.")
	}
	rez.df <- data.frame(meta, rez[ rownames(meta),1:K ])
	
	## copy over the %variance explained
	class(rez.df) <- c("pca.result", class(rez.df))
	attr(rez.df, "explained") <- attr(rez, "explained")[1:K]
	if ("effective.dim" %in% names(attributes(rez)))
		attr(rez.df, "effective.dim") <- attr(rez, "effective.dim")
	
	return(rez.df)
	
}
#' @export
pca <- function(x, ...) UseMethod("pca")
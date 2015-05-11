## pca.R
## functions for performing PCA on genotypes and hybridization intensities

## internal function for actually doing PCA
.do.pca.genotypes <- function(gty, extras = NULL, what = c("genotypes","intensity"),
							  scale = TRUE, center = TRUE, eps = 1e-6, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (!is.null(extras) && !inherits(extras, "genotypes"))
		stop(paste("If supplying extra samples to project onto pre-computed PCs",
				   "they should also be suppled as class 'genotypes'."))
	
	what <- match.arg(what)
	.make.pca.input <- function(g) {
		if (what == "genotypes") {
			if (!is.numeric(gty)) {
				gty <- recode.genotypes(gty, "01")
			}
			X <- t(.copy.matrix.noattr(gty))
			## replace missing values with column mean
			if (any(is.na(colMeans(X)))) {
				message("\treplacing missing values with minor-allele frequency...")
				X <- apply(X, 2, function(x) {
					x[ is.na(x) ] <- mean(x, na.rm = TRUE)
					return(x)
				})
			}
		}
		else if (.has.valid.intensity(gty)) {
			x <- with(attr(gty, "intensity"), .copy.matrix.noattr(x))
			y <- with(attr(gty, "intensity"), .copy.matrix.noattr(y))
			rownames(x) <- paste0(rownames(x), "_x")
			rownames(y) <- paste0(rownames(x), "_y")
			x[ is.na(x) ] <- 0
			y[ is.na(y) ] <- 0
			X <- t(rbind(x,y))	
		}
		else {
			stop("Need valid intensity matrices in order to do PCA on intensity.")
		}
		
		## remove zero-variance columns, if scaling requested
		if (scale) {
			const <- abs(colVars(X) - 0) <= eps
			const[ is.na(const) ] <- TRUE
			X <- X[ ,!const ]
		}
		return(X)
	}

	## finally do PCA
	message(paste("Computing principal components of", what, "matrix..."))
	pc <- prcomp(.make.pca.input(gty), center = center, scale. = scale)
	proj <- predict(pc)
	if (!is.null(extras)) {
		message("Projecting extra samples onto existing PCs...")
		topredict <- .make.pca.input(extras)
		predicted <- predict(pc, newdata = topredict)
		proj <- rbind(proj, predicted)
	}
	
	## add %variance explained
	attr(proj, "pca") <- pc
	attr(proj, "explained") <- pc$sdev^2/sum(pc$sdev^2)
	
	message("Done.")
	return(proj)
	
}

#' Perform PCA on a \code{genotypes} object
#'
#' Performs principal components analysis (PCA) on either (numerically-coded) genotypes or
#' on the underliny 2D hybridization-intensity matrices.
#'
#' @param gty an object of class \code{genotypes}
#' @param extras a second dataset of class \code{genotypes}, to be projected onto PCs computed
#' 	from the first. NB: this function *does not* verify that the underlying features (ie. markers) match.
#' @param what \code{"genotypes"} to do PCA on genotypes (coded 0/1/2); \code{"intensity"} to do PCA on
#' 	underlying 2D intensities.  The latter triggers an error if intensity matrices are absent.
#' @param K how many PCs to return.
#' @param \code{...} other parameters for call to \code{prcomp}, such as \code{center} and \code{scale}
#' 	(both \code{TRUE} by default)
#'
#' @details Uses base-\code{R}'s \code{prcomp} under the hood. By default, columns are centered and scaled;
#' 	not doing so will likely produce a strange result.  When doing PCA on genotypes, missing values are
#' 	replaced with the column mean, which in many circumstances can be interpreted as the minor-allele frequency.
#' 	(This is very similar to the behavior of PLINK.)  When doing PCA on intensities, missing values are set
#' 	to zero -- but even no-call genotypes have nonmissing intensities on Illumina arrays, so this is unlikely
#' 	to have any effect in practice.
#'
#' @return a dataframe with as many rows as samples, in which the first columns are sample IDs and any associated
#' 	metadata (as returned by \code{samples(gty)}), followed by the first K PCs.  Scaled eigenvalues (ie. percent
#' 	of variance explained) are provided as \code{attr(,"explained")}.
#'
#' @seealso \code{pca.plink()} for using PLINK's (much faster and more powerful) implementation
#'
#' @examples
#' pca.genotypes(geno)
#' # equivalently, as shortcut:
#' pca(geno)
#'
#' @export
pca.genotypes <- function(gty, extras = NULL, what = c("genotypes","intensity"), K = 3, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (!is.numeric(K) || K < 1)
		stop("You probably don't want to do PCA with <1 dimension.")
	
	## run the PCA
	rez <- .do.pca.genotypes(gty, extras = extras, what = what, ...)
	
	## add sample metatdata, if any
	meta <- data.frame(iid = rownames(rez))
	rownames(meta) <- rownames(rez)
	if (.has.valid.ped(gty)) {
		meta <- merge(meta, attr(gty, "ped"), by = "iid")
		rownames(meta) <- as.character(meta$iid)
		if (!(nrow(meta) == nrow(rez) && all(rownames(rez) %in% meta$iid)))
			stop("Sample metadata and PCA result are out of sync.")
	}
	rez.df <- data.frame(meta, rez[ rownames(meta),1:K ])
	
	## copy over the %variance explained
	attr(rez.df, "explained") <- attr(rez, "explained")[1:K]
	
	return(rez.df)
	
}
pca <- function(x, ...) UseMethod("pca")
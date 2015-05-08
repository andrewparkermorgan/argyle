## normalize.R

## use preprocessCore::normalize.quantiles() to do quantile normalization
quantile.normalize <- function(geno, weights = NULL, force = FALSE, ...) {
	
	require(preprocessCore)
	
	if (!(inherits(geno, "genotypes") && is.list(attr(geno, "intensity"))))
		stop("Please supply an object of class 'genotypes' with intensity information attached.")
	
	if (is.null(attr(geno, "normalized")))
		attr(geno, "normalized") <- FALSE
	
	if (attr(geno, "normalized") && !force) {
		message("Intensities seem to already be normalized")
		return(geno)
	}
	else
		warning("Intensities seem to already be normalized; doing it again won't help.")
	
	if (!is.null(weights) && length(weights) != ncol(geno))
		stop("Dimensiosn of weight vector and intensity matrices don't match.")
	
	message(paste("Performing robust quantile normalization with", ifelse(is.null(weights), "no weights","weights"), "..."))
	x.norm <- normalize.quantiles.robust(attr(geno, "intensity")$x, weights = weights)
	y.norm <- normalize.quantiles.robust(attr(geno, "intensity")$x, weights = weights)
	
	attr(geno, "intensity") <- list(x = x.norm, y = y.norm)
	attr(geno, "normalized") <- TRUE
	return(geno)
	
}
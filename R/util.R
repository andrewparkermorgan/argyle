## util.R
## mostly internal utility functions


## copy matrix and keep row/column names, but drop other attributes (including class!)
.copy.matrix.noattr <- function(x, ...) {
	
	if (!is.matrix(x))
		stop("Input not a matrix.")
	
	rez <- matrix(as.vector(x), ncol = ncol(x), nrow = nrow(x))
	colnames(rez) <- colnames(x)
	rownames(rez) <- rownames(x)
	return(rez)
	
}

## fast computation of variance column-wise, without looping and using fast internals
## thanks to TS Lumley
colVars <- function(x, na.rm = TRUE, warn.missing = FALSE, ...) {
	
	if (!is.matrix(x))
		stop("Input not a matrix.")
	
	n <- nrow(x)
	sigmasq <- ( n/(n-1) * (colMeans(x*x, na.rm = na.rm)-colMeans(x, na.rm = na.rm)^2) )
	if (warn.missing && is.na(sigmasq))
		return(0)
	else
		return(sigmasq)
	
}

## declare an object a 'genotypes'
bless <- function(x) {
	if (!inherits(x, "genotypes"))
		class(x) <- c("genotypes", class(x))
	return(x)
}
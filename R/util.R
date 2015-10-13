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
colVars <- function(x, na.rm = TRUE, ...) {
	
	if (!is.matrix(x))
		stop("Input not a matrix.")
	
	n <- nrow(x)
	return( n/(n-1) * (colMeans(x*x, na.rm = na.rm)-colMeans(x, na.rm = na.rm)^2) )
	
}

## declare an object a 'genotypes'
bless <- function(x) {
	if (!inherits(x, "genotypes"))
		class(x) <- c("genotypes", class(x))
	return(x)
}
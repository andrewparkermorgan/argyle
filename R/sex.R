## sex.R
## predict/check sex of samples
## relies on some stuff in package mclust

## very simple threshold prediction of sex, based on count of good calls on chrY
.predict.sex.ycalls <- function(gty, max.N = 30, min.AB = 40, ...) {
	
	## get calls from chrY markers
	gty <- subset(gty, chr == "chrY" | toupper(chr) == "Y")
	cc <- summarize.calls(gty, verbose = FALSE)
	
	## apply thresholds
	sexes <- setNames( rep(0, ncol(gty)), colnames(gty) )
	sexes[ (cc$A + cc$B >= min.AB) & (cc$N <= max.N) ] <- 1
	sexes[ (cc$A + cc$B < min.AB) & (cc$N > max.N) ] <- 2
	
	return(sexes)
	
}

#' Predict sample sexes based on genotype and intensity data
#' 
#' @param object a \code{genotypes} object
#' @param method how to go about making sex predictions (see Details)
#' @param ... other parameters passed to underlying prediction functions
#' 
#' @return a dataframe with 4 columns: individual ID, nominal sex (0 if unknown), predicted sex
#' 	(0 if ambiguous), probability (NA if not a model-based prediction)
#' 	
#' @details Implements several (soon) methods for predicting the sex of a sample given genotype calls
#' 	and hybridization intensity on the sex chromosomes.  Assumes the sex-chromosome system of eutherian
#' 	mammals (female karyotype XX, male karyotype XY).  Sex chromosomes should be named "X" ("chrX") and
#' 	"Y" ("chrY") respectively.
#' 	
#' 	Method \code{"ycalls"} simply counts the number of missing and of non-missing, non-heterozygous calls
#' 	at Y-linked markers, and applies a threshold to both.  The defaults are calibrated to the GigaMUGA
#' 	array for mouse.  Females have mostly missing calls, males have mostly non-missing calls, and no
#' 	sample should have many heterozygous calls.
#' 
#' @export predict.sex
predict.sex <- function(object, method = c("ycalls"), ...) {
	
	if (!(inherits(object, "genotypes") && .has.valid.map(object)))
		stop("Please supply an object of class 'genotypes' with valid marker map.")
	
	if (method == "ycalls") {
		message("Predicting sex using count of good calls on chrY...")
		sexes <- .predict.sex.ycalls(object, ...)
		prob <- setNames( rep(NA, ncol(object)), colnames(object) )
	}
	else {
		stop("Other sexing methods not implemented yet.")
	}
	
	if (.has.valid.ped(object) && "sex" %in% colnames(attr(object, "ped"))) {
		nominal <- setNames( attr(object, "ped")$sex, rownames(attr(object, "ped")) )
	}
	else {
		nominal <- setNames( rep(0, ncol(object)), colnames(object) )
	}
	
	rez <- data.frame(iid = colnames(object), nominal = nominal[ colnames(object) ],
					  predicted = sexes[ colnames(object) ], prob = prob[ colnames(object) ])
	return(rez)
	
}

sexysum <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	gty <- subset(gty, chr == "chrX" | chr == "chrY")
	bychr <- genoapply(gty, 1, .(chr), summarize.calls)
	calls <- plyr::ldply(bychr)
	sexing <- plyr::ddply(calls, plyr::.(iid), plyr::summarise,
						  Ygood = A[ chr == "chrY" ]+B[ chr == "chrY" ],
						  Xhet = H[ chr == "chrX" ])
	sexing <- merge(sexing, attr(gty, "ped"))
	rownames(sexing) <- as.character(sexing$iid)
	return(sexing)
	
}
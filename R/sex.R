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

.predict.sex.xy <- function(gty, platform = "giga", clean = TRUE, ...) {
	
	## this ugly thing is the cluster model for GigaMUGA
	models <- structure(list(giga = structure(list(
		mu = structure(
			c(1.00179688827984, 0.579627096415587, 0.816608462699492, 1.00517965688652),
			.Dim = c(2L, 2L), .Dimnames = list(c("x", "y"), c("female", "male"))),
		sigma = structure(list(female = structure(c(0.000777449983304451, 0.000646632991729865, 0.000646632991729865, 0.00911353377861309),
												  .Dim = c(2L, 2L), .Dimnames = list(c("x", "y"), c("x", "y"))),
							   male = structure(c(0.000343819755693033, 0.00036077907454159, 0.00036077907454159, 0.0107555733491657),
							   				 .Dim = c(2L, 2L), .Dimnames = list(c("x", "y"), c("x","y")))), .Names = c("female", "male"))),
		.Names = c("mu", "sigma"))), .Names = "giga")
	
	
	.plat <- match.arg(platform)
	if (!(.plat %in% names(models)))
		stop("Unsupported platform.")
	
	## extract intensities for each sex chrom
	ychr <- subset(gty, grepl("Y", chr))
	if (clean)
		ychr <- subset(ychr, grepl("^JAX", marker))
	
	xchr <- subset(gty, grepl("X", chr))
	
	## matrix of samples x (X,Y)
	xy <- cbind( x = colMeans(.si(xchr), na.rm = TRUE),
				 y = colMeans(.si(ychr), na.rm = TRUE) )
	rownames(xy) <- colnames(gty)
	
	## evaluate probs
	themodel <- models[[.plat]]
	pmale <- mvtnorm::dmvnorm(xy, mean = themodel$mu[ ,"male" ], sigma = themodel$sigma$male, log = TRUE)
	pfemale <- mvtnorm::dmvnorm(xy, mean = themodel$mu[ ,"female" ], sigma = themodel$sigma$female, log = TRUE)
	
	score <- pmale-pfemale
	names(score) <- colnames(gty)
	sexes <- ifelse(score > 0, 1, 2)
	names(sexes) <- colnames(gty)
	
	return(list(sex = sexes, score = score))
	
}

#' Predict sample sexes based on genotype and intensity data
#' 
#' @param object a \code{genotypes} object
#' @param method how to go about making sex predictions (see Details)
#' @param clean logical; for *MUGA arrays, use only known-good Y chromosome markers
#' @param platform character; name of a specific array platform (used only for \code{"xy"} method)
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
#' 	Method \code{"xy"} calculates the sum-intensity on each sex chromsome and evaluates them against
#' 	a pre-constructed set of clusters in 2d space corresponding to each sex.  Only available for the
#' 	GigaMUGA array for mouse, at present.
#' 
#' @export predict.sex
predict.sex <- function(object, method = c("ycalls","xy"), clean = TRUE, platform = "giga", ...) {
	
	if (!(inherits(object, "genotypes") && .has.valid.map(object)))
		stop("Please supply an object of class 'genotypes' with valid marker map.")
	
	if (method == "ycalls") {
		message("Predicting sex using count of good calls on chrY...")
		sexes <- .predict.sex.ycalls(object, ...)
		prob <- setNames( rep(NA, ncol(object)), colnames(object) )
	}
	else if (method == "xy") {
		assigner <- .predict.sex.xy(object, clean = clean, platform = platform, ...)
		sexes <- assigner$sex
		prob <- assigner$score
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

#' Summary stats for the sex chromosomes
#' 
#' @param gty a \code{genotypes} object
#' @param clean logical; for *MUGA arrays, use only known-good Y chromosome markers
#' @param ... ignored
#' 
#' @return A dataframe of all pre-existing sample metadata with additional columns \code{Ygood} (count of good calls on the
#' 	Y chromosome), \code{Xhet} (count of heterozygous calls on the X chromosome), \code{ix} (mean intensity of X markers) and
#' 	\code{iy} (mean intensity of Y markers).  If no intensity data is availalble, the latter two columns are not included.
#' 	
#' @details If genotypes are not present for *both* sex chromosomes, the results are unpredictable.
#' 
#' @export
sexysum <- function(gty, clean = TRUE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	## strip all but sex chroms
	gty <- subset(gty, grepl("X", chr) | grepl("Y", chr))
	
	## count calls on sex chroms
	bychr <- genoapply(gty, 1, .(chr), summarize.calls)
	calls <- plyr::ldply(bychr)
	sexing <- plyr::ddply(calls, plyr::.(iid), plyr::summarise,
						  Ygood = A[ grep("X", chr)[1] ]+B[ grep("Y", chr)[1] ],
						  Xhet = H[ grep("X", chr)[1] ])
	sexing <- merge(sexing, attr(gty, "ped"))
	rownames(sexing) <- as.character(sexing$iid)
	
	## get sum-intensities, if available
	if (.has.valid.intensity(gty)) {
		ys <- with(markers(gty), grepl("Y", chr))
		if (clean)
			ys <- ys & with(markers(gty), grepl("^JAX", marker))
		siy <- colMeans( .si(subset(gty, ys)), na.rm = TRUE )
		six <- colMeans( .si(subset(gty, grepl("X", chr))), na.rm = TRUE )
		sexing$ix <- six[ rownames(sexing) ]
		sexing$iy <- siy[ rownames(sexing) ]
	}
	
	return(sexing)
	
}
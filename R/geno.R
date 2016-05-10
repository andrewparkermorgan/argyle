## geno.R
## utility functions for handling genotype matrices

#' Indexing into a \code{genotypes} object
#'
#' Allows indexing into a genotype matrix, trimming marker and sample metadata, filters,
#' and intensity matrices accordingly.
#'
#' @param x an object of class \code{genotypes}
#' @param i row index (logical, character, or numeric), recycled if needed
#' @param j column index (logical, character, or numeric), recycled if needed
#' @param ... ignored
#' @param drop allow reduction of the dimensions of \code{x} if \code{i,j} are of length 1?
#'
#' @return a \code{genotypes} object, preserving any attributes
#'
#' @examples
#' \dontrun{
#' geno[ 1:10, ]
#' geno[ ,c("sample1","sample2") ]
#'}
#' @export
`[.genotypes` <- function(x, i, j, ..., drop = FALSE) {
	
	if (missing(i))
		i <- rep(TRUE, nrow(x))
	if (missing(j))
		j <- rep(TRUE, ncol(x))
	
	r <- NextMethod("[", drop = drop)
	## handle "array indexing"
	#if (is.matrix(i)) {
	#	ii <- i
	#	i <- ii[,1, drop = TRUE]
	#	j <- ii[,2, drop = TRUE]
	#}
	i <- rownames(r)
	j <- colnames(r)
	if (!is.null(attr(x, "map"))) {
		attr(r, "map") <- attr(x, "map")[ i,,drop = FALSE ]
	}
	if (!is.null(attr(x, "ped"))) {
		attr(r, "ped") <- attr(x, "ped")[ j,,drop = FALSE ]
		rownames(attr(r, "ped")) <- as.character(attr(r, "ped")$iid)
	}
	if (.has.valid.intensity(x)) {
		x.new <- attr(x, "intensity")$x[ i,j, drop = FALSE ]
		y.new <- attr(x, "intensity")$y[ i,j, drop = FALSE ]
		if (any(rownames(x.new) != rownames(r)))
			stop("Rownames of genotypes and intensity matrices don't match.")
		attr(r, "intensity") <- list(x = x.new, y = y.new)
		if (!is.null(attr(x, "normalized")))
			attr(r, "normalized") <- attr(x, "normalized")
	}
	if (.has.valid.baflrr(x)) {
		baf.new <- attr(x, "baf")[ i,j, drop = FALSE ]
		lrr.new <- attr(x, "lrr")[ i,j, drop = FALSE ]
		attr(r, "baf") <- baf.new
		attr(r, "lrr") <- lrr.new
	}
	
	if (!is.null(attr(x, "alleles")))
		attr(r, "alleles") <- attr(x, "alleles")
	
	if (!is.null(attr(x, "filter.sites")))
		attr(r, "filter.sites") <- attr(x, "filter.sites")[i]
	if (!is.null(attr(x, "filter.samples")))
		attr(r, "filter.samples") <- attr(x, "filter.samples")[j]
	
	class(r) <- c("genotypes", class(r))
	return(r)
}

#' Access attributes of a \code{genotypes} object as if they were a list
#'
#' @param x an object of class \code{genotypes}
#' @param expr name of the attribute to extract
#' @param ... ignored
#' 
#' @export
`$.genotypes` <- function(x, expr, ...) {
	
	attributes(x)[[expr]]
	
}


#' @export
summary.genotypes <- function(object, ...) {
	
	nsamples <- ncol(object)
	nsites <- nrow(object)
	cat("---", deparse(substitute(object)), "---\nA genotypes object with", nsites, "sites x", nsamples, "samples\n")
	cat("Allele encoding:", attr(object, "alleles"), "\n")
	
	has.intens <- .has.valid.intensity(object)
	normed <- .null.false(attr(object, "normalized"))
	cat("Intensity data:", ifelse(has.intens, "yes","no"),
		ifelse(has.intens, ifelse(normed, "(normalized)", "(raw)"), ""), "\n")
	cat("Sample metadata:", ifelse(.has.valid.ped(object), "yes", "no"))
	if (.has.valid.ped(object) && "sex" %in% colnames(attr(object, "ped"))) {
		sex <- attr(object, "ped")$sex
		sexes <- tabulate(sex, nbins = 2)
		cat(" (", sexes[1], "male /", sexes[2], "female /", length(sex)-sum(sexes), "unknown )\n")
	}
	else {
		cat("\n")
	}
	
	filt <- colSums( summarize.filters(object) )
	cat("Filters set:", filt[1], "sites /", filt[2], "samples", "\n")
	
	if (!is.null(attr(object, "source")))
		cat("File source:", attr(object, "source"), "(on", format(attr(object, "timestamp")), ")\n")
	if (!is.null(attr(object, "md5")))
		cat("Checksum:", attr(object, "md5"), "\n")
	
	if (0) {
		if (.has.valid.map(object)) {
			cat("\n---\nSites by chromosome:\n")
			tbl.chr <- table(attr(object, "map")$chr)
			for (i in seq_along(tbl.chr)) {
				cat("\t", names(tbl.chr)[i], ":", tbl.chr[i], "\n")
			}
			cat("\n")
		}
	}
	
}


#' Extract vector of sample sexes
#'
#' @param gty an object of class \code{genotypes}
#' 
#' @return named numeric vector giving sample sexes: 0=unknown, 1=male, 2=female
#' 
#' @export
sex.genotypes <- function(gty) {
	if (.has.valid.ped(gty))
		rez <- setNames( attr(gty, "ped")$sex, colnames(gty) )
	else
		rez <- setNames( rep(0, ncol(gty)), colnames(gty) )
	class(rez) <- c("sexes", class(rez))
	return(rez)
}
sex <- function(gty) UseMethod("sex")

#' Print summary table of sample sexes
#' 
#' Counts the number of samples of each sex, converts numeric codes to human-friendly
#' labels, and returns the table.
#'
#' @param gty an object of class \code{genotypes}
#' @param ... ignored
#' 
#' @return table of sample counts by sex
#' 
#' @export
sextable <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply and object of class 'genotypes'.")
	
	s <- sex.genotypes(gty)
	x <- tabulate(s+1, nbins = 3)
	x <- x[ c(2,3,1) ]
	names(x) <- c("male","female","unknown")
	return(x)
}

#' @export
print.genotypes <- function(x, ...) {
	summary.genotypes(x)
	
	if (.has.valid.map(x)) {
		chrs <- factor(attr(x, "map")$chr)
		cat("\nCounts of markers by chromosome:\n")
		tbl <- table(chrs, dnn = NULL)
		print(tbl)
		cat("\n")
	}
	
}

#' Show the first few samples and markers for a \code{genotypes}
#'
#' @param x an object of class \code{genotypes}
#' @param n number of markers to show
#' @param nsamples number of samples to show
#' @param ... ignored
#' 
#' @details Sample names are abbreviated to keep output tidy.
#' 
#' @aliases head
#' @export
head.genotypes <- function(x, n = 10, nsamples = min(10,ncol(x)), ...) {
	
	if (!inherits(x, "genotypes"))
		stop("Please supply and object of class 'genotypes'.")
	
	w <- 5
	head.fmt <- paste0("%", w, "s")
	col.fmt <- paste0("%", w, "s")

	n <- min(nrow(x), n)
	nc <- min(ncol(x), nsamples)
	if (n >= 1 && nc >= 1) {
		
		x <- x[1:n,1:nc]
		G <- matrix( as.character(.copy.matrix.noattr(x)), nrow = nrow(x), ncol = ncol(x) )
		sm <- sprintf(head.fmt, substr(gsub("[\\_\\. \\:\\(\\)]", "", colnames(x)), 1, w))
		mk <- sprintf("%15s", rownames(x))
		
		cat("Genotypes matrix:\n")
		cat(sprintf("%15s", ""), sm, "\n")
		for (i in 1:nrow(G)) {
			cat(mk[i], paste(sprintf(col.fmt, G[i,])), "\n")
		}
		
		if (.has.valid.map(x)) {
			map <- attr(x, "map")[1:n,]
			cat("\nMarker map:\n")
			print(map[ ,1:min(6, ncol(map)) ], row.names = FALSE)
		}
		
		if (.has.valid.ped(x)) {
			cat("\n")
			ped <- attr(x, "ped")[1:nc,]
			cat("Sample info:\n")
			print(ped, row.names = FALSE)
		}
		
	}	
	else {
		## nothing much here...
		summary.genotypes(x)
	}

}

## set some S3 generics for useful accessor functions

#' Get marker map for a \code{genotypes} object
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' @return a dataframe of marker information (chr, cM, pos, ...) which should run parallel to the
#' 	rows of the genotypes matrix in the input.  If map is absent, returns \code{NULL.}
#' @seealso Other accessor functions: \code{\link{samples}} (sample metadata), \code{\link{filters}}
#' 	(site and sample filters), \code{\link{intensity}} (intensity matrices)
#' @aliases markers
#' @export
markers.genotypes <- function(gty, ...) {
	attr(gty, "map")
}
#' @export
markers <- function(gty, ...) UseMethod("markers")

#' Get sample metadata for a \code{genotypes} object
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' @return a dataframe of sample metadata (iid, ...) which should run parallel to the
#' 	columns of the genotypes matrix in the input.  If sample data is absent, returns \code{NULL.}
#' @seealso Other accessor functions: \code{\link{markers}} (marker map), \code{\link{filters}}
#' 	(site and sample filters), \code{\link{intensity}} (intensity matrices)
#' @aliases samples
#' @export
samples.genotypes <- function(gty, ...) {
	if (!is.null(attr(gty, "ped")))
		attr(gty, "ped")
	else
		colnames(gty)
}
#' @export
samples <- function(gty, ...) UseMethod("samples")

#' Get filters attached to a \code{genotypes} object
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' @return list of length 2 with elements \code{sites} and \code{samples} which should run parallel to the
#' 	row and columns, respectively, of the genotypes matrix in the input.  If fitlers are absent, returns vector
#' 	with appropriate dimensions in which all elements are empty.  These vectors have names matching the row
#' 	and column names of the genotypes matrix.
#' @seealso Other accessor functions: \code{\link{markers}} (marker map), \code{\link{samples}}
#' 	(sample metadata), \code{\link{intensity}} (intensity matrices)
#' @aliases filters
#' @export
filters.genotypes <- function(gty, ...) {
	get.filters(gty, ...)
}
#' @export
filters <- function(gty, ...) UseMethod("filters")

#' Get intensity matrices attached to a \code{genotypes} object
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' @return list of length 2 with elements \code{x} and \code{y}, respectively the x- and y-intensity matrices,
#' 	which both have dimensions and row/column names equal to those of the genotypes matrix in the input.
#' @seealso Other accessor functions: \code{\link{markers}} (marker map), \code{\link{samples}}
#' 	(sample metadata), \code{\link{filters}} (site and sample filters)
#' @aliases intensity
#' @export
intensity.genotypes <- function(gty, ...) {
	attr(gty, "intensity")
}
#' @export
intensity <- function(gty, ...) UseMethod("intensity")

#' Check if a \code{genotypes} object has intensity data attached
#'
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' 
#' @return logical scalar; TRUE is intensity data is present
#'
#' @export has.intensity
has.intensity <- function(gty, ...) {
	.has.valid.intensity(gty)
}

#' Extract hybridization intensities as a dataframe
#'
#' @param gty an object of class \code{genotypes}
#' @param markers indexing vector to choose which markers (rows) to extract
#' @param ... ignored
#' 
#' @return a dataframe containing sample metadata, marker metadata, and hybridization intensities
#' 
#' @export get.intensity
get.intensity <- function(gty, markers, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes' with intensity information attached.")
	
	rez <- reshape2:::melt.matrix(attr(gty, "intensity")$x[ markers,,drop = FALSE ])
	colnames(rez) <- c("marker","iid","x")
	rez <- cbind(rez, y = reshape2:::melt.matrix(attr(gty, "intensity")$y[ markers,,drop = FALSE ])[ ,3 ])
	rez$si <- with(rez, sqrt(x^2+y^2))
	
	if (.has.valid.map(gty))
		rez <- merge(rez, attr(gty, "map"))
	
	if (.has.valid.ped(gty))
		rez <- merge(rez, attr(gty, "ped"))
	
	return(rez)
	
}

#' Rename some or all samples
#' 
#' @param gty a \code{genotypes} object
#' @param nn a vector of new sample names
#' @param ... ignored
#' 
#' @details The values in \code{nn} are used to replace existing names using the following rules. First, if \code{nn}
#' 	is a named vector, it is assumed to map old names to new names.  Otherwise \code{nn} must have length equal to the
#' 	number of samples in \code{gty}, and each element should be a new name to replace an old one.  Missing values are
#' 	not permitted, for obvious reasons.
#' 
#' @export replace.names
replace.names <- function(gty, nn, ...) {

	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (length(nn) < 1 || any(is.na(nn)))
		stop("Names cannot have zero length or contain missing values.")
	
	## if 'nn' is named vector, assume it maps old names -> new names
	if (!is.null(names(nn))) {
		olds <- intersect(colnames(gty), names(nn))
		same <- setdiff(colnames(gty), olds)
		same <- setNames(same, same)
		nn <- c(nn, same)[ colnames(gty) ]
	}
	## otherwise assume 'nn' is parallel to columns of genotype matrix
	else if (length(nn) == ncol(gty)) {
		nn <- setNames(nn, colnames(gty))
	}
	else {
		stop("New sample names in 'nn' should be either named, or have same length as ",
			 "the number of samples.")
	}
	
	## check for duplicates
	on <- colnames(gty)
	rn <- nn[on]
	dups <- duplicated(rn)
	if (any(dups)) {
		warning("The following sample names are duplicates; they will be uniquified, but you ",
				"should check what's going on: ", paste(nn[dups], collapse = ","))
	}
	rn <- make.unique(rn)
	
	## now actually do the renaming
	colnames(gty) <- rn
	
	if (.has.valid.ped(gty)) {
		attr(gty, "ped")$iid <- as.character(attr(gty, "ped")$iid)
		attr(gty, "ped")$iid <- nn[on]
		rownames(attr(gty, "ped")) <- attr(gty, "ped")$iid
	}
	
	if (.has.valid.intensity(gty)) {
		colnames(attr(gty, "intensity")$x) <- rn
		colnames(attr(gty, "intensity")$y) <- rn
	}
	
	if (.has.valid.baflrr(gty)) {
		colnames(attr(gty, "baf")) <- rn
		colnames(attr(gty, "lrr")) <- rn
	}
	
	if (!is.null(attr(gty, "filter.samples")))
		names(attr(gty, "filter.samples")) <- rn
	
	return(gty)
	
}

#' Extract genotype calls as a dataframe
#'
#' @param gty an object of class \code{genotypes}
#' @param markers indexing vector to choose which markers (rows) to extract
#' @param ... ignored
#' 
#' @return a dataframe containing sample metadata, marker metadata, and genotype calls
#' 
#' @export get.call
get.call <- function(gty, markers, ...) {
	
	if (!(inherits(gty, "genotypes")))
		stop("Please supply an object of class 'genotypes' with intensity information attached.")
	
	rez <- reshape2:::melt.matrix(gty[ markers,,drop = FALSE ])
	colnames(rez) <- c("marker","iid","call")
	
	if (.has.valid.map(gty))
		rez <- merge(rez, attr(gty, "map"))
	
	if (.has.valid.ped(gty))
		rez <- merge(rez, attr(gty, "ped"))
	
	return(rez)
	
}

#' Extract BAF/LRR as a dataframe
#'
#' @param gty an object of class \code{genotypes}
#' @param markers indexing vector to choose which markers (rows) to extract
#' @param ... ignored
#' 
#' @return a dataframe containing sample metadata, marker metadata, and BAF/LRR values
#' 
#' @details BAF and LRR must have already been computed using the \code{\link{tQN}} function.
#' 
#' @export get.baf
get.baf <- function(gty, markers = NULL, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.baflrr(gty)))
		stop("Please supply an object of class 'genotypes' with BAF and LRR computed.")
	
	if (is.null(markers))
		markers <- TRUE
	
	rez <- reshape2:::melt.matrix(attr(gty, "baf")[ markers,,drop = FALSE ])
	colnames(rez) <- c("marker","iid","BAF")
	rez <- cbind(rez, LRR = reshape2:::melt.matrix(attr(gty, "lrr")[ markers,,drop = FALSE ])[ ,3 ])
	rez <- cbind(rez, call = reshape2:::melt.matrix(gty[ markers,,drop = FALSE ])[ ,3 ])
	
	if (.has.valid.map(gty))
		rez <- merge(rez, attr(gty, "map"))
	
	if (.has.valid.ped(gty))
		rez <- merge(rez, attr(gty, "ped"))
	
	return(rez)
	
}

#' Subset a \code{genotypes} object by markers or samples
#' 
#' @param x a \code{genotypes} object
#' @param expr logical expression indicating elements or rows to keep: missing values are taken as false
#' 	(just like \code{subset.data.frame()})
#' @param by should \code{expr} be evaluated in the context of markers (eg. by rows) or samples (eg. by columns)?
#' @param ... ignored
#' 
#' @return a \code{genotypes} object similar to \code{x}, but with only the selected rows (\code{by == "markers"})
#' 	or columns (\code{by == "samples"}).  Attributes of \code{x} including intensity matrices, filters,
#' 	marker map and sample metadata are all adjusted accordingly.
#' 	
#' @section Warning:
#' Since this function uses non-standard evaluation, be careful about expressions which include variables defined
#' in both the global environment and in the scope of the call (eg. variables named same as columns of the marker
#' map).  In that case results may not be what you expect.
#' 
#' @export
subset.genotypes <- function(x, expr, by = c("markers","samples"), ...) {
	
	if (!inherits(x, "genotypes"))
		stop("Only willing to subset an object of class 'genotypes.'")
	
	by <- match.arg(by)
	e <- substitute(expr)
	if (by == "markers") {
		if (!.has.valid.map(x))
			stop("Can't subset genotypes by marker without a marker map.")
		r <- eval(e, attr(x, "map"), parent.frame())
		r <- r & !is.na(r)
		return( x[ r, ] )
	}
	else if (by == "samples") {
		if (is.null(attr(x, "ped")))
			stop("Can't subset genotypes by samples without sample information.")
		r <- eval(e, attr(x, "ped"), parent.frame())
		r <- r & !is.na(r)
		return( x[ ,r ] )
	}
	else {
		stop()
	}
	
}

#' Shortcut for grabbing just the autosomes
#' 
#' @param gty a \code{genotypes} object
#' @param ... ignored
#'  
#' @return a copy of \code{gty} without chrX, chrY or chrM (or markers not assigned to a chromosome)
#' 
#' @export
autosomes <- function(gty, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.map(gty)))
		stop("Please supply an object of class 'genotypes' with valid marker map.")

	gty[ !grepl("[YXMPUu]", attr(gty, "map")$chr) & grepl("[0-9]+", attr(gty, "map")$chr), ]
	
}

#' Shortcut for grabbing just chrX
#' 
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' 
#' @return a copy of \code{gty} with only chrX
#'
#' @export
xchrom <- function(gty, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.map(gty)))
		stop("Please supply an object of class 'genotypes' with valid marker map.")
	
	gty[ grepl("[X]", attr(gty, "map")$chr), ]
	
}

#' Shortcut for grabbing just chrY
#' 
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' 
#' @return a copy of \code{gty} with only chrY
#'
#' @export
ychrom <- function(gty, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.map(gty)))
		stop("Please supply an object of class 'genotypes' with valid marker map.")
	
	gty[ grepl("[Y]", attr(gty, "map")$chr), ]
	
}

#' Concatenate two \code{genotypes} objects 'horizontally'
#'
#' @param a an object of class \code{genotypes}
#' @param b another object of class \code{genotypes}
#' @param ... ignored
#' 
#' @details The objects \code{a} and \code{b} must represent the same markers typed on non-overlapping sets
#' 	of samples in order for this operation to be meaningful. That means their row names must match.
#' 	Any sample or site filters are dropped in the result.
#' 
#' @export
cbind.genotypes <- function(a, b = NULL, ...) {
	
	if (is.null(b))
		return(a)
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to bind two objects of class 'genotypes.'")
	
	if (nrow(a) != nrow(b) | any(rownames(a) != rownames(b)))
		stop("Number and names of markers don't match.  Try merging genotypes instead of cbind-ing.")
	
	dups <- intersect(colnames(a), colnames(b))
	if (length(dups))
		stop( paste("Duplicate sample names encountered; make them unique before cbind-ing:", paste0(dups, collapse = ",")) )
	
	if (any(is.null(attr(a,"alleles")), is.null(attr(b,"alleles"))))
		stop("Ambiguous allele encoding in one or both of the objects to be merged.")
	if (attr(a,"alleles") != attr(b,"alleles"))
		stop("Allele encodings don't match -- nonsensical result.")
	
	message(paste0("Adding ",ncol(b)," individuals to the existing ",ncol(a),"."))
	
	rez <- cbind(unclass(a), unclass(b))
	class(rez) <- c("genotypes", class(rez))
	if (!is.null(attr(a, "map")))
		attr(rez, "map") <- attr(a, "map")
	if (!is.null(attr(a, "ped")))
		attr(rez, "ped") <- rbind( attr(a, "ped"), attr(b, "ped") )
	attr(rez, "alleles") <- attr(a, "alleles")
	
	if (.has.valid.intensity(a) && .has.valid.intensity(b)) {
		message("Updating intensity matrices...")
		x <- cbind( attr(a, "intensity")$x, attr(b, "intensity")$x )
		y <- cbind( attr(a, "intensity")$y, attr(b, "intensity")$y )
		attr(rez, "intensity") <- list(x = x, y = y)
	}
	
	return(rez)
	
}

#' Concatenate two \code{genotypes} objects 'vertically'
#'
#' @param a an object of class \code{genotypes}
#' @param b another object of class \code{genotypes}
#' @param ... ignored
#' 
#' @details The objects \code{a} and \code{b} must represent the same samples at non-overlapping sets
#' 	of markers in order for this operation to be meaningful. That means their column names must match.
#' 	Sample or site filters are dropped in the result.
#' 
#' @export
rbind.genotypes <- function(a, b, ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to bind two objects of class 'genotypes.'")
	
	dups <- intersect(rownames(a), rownames(b))
	if (length(dups))
		stop( paste("Duplicate marker names encountered; make them unique before rbind-ing:", paste0(dups, collapse = ",")) )
	
	if (any(is.null(attr(a,"alleles")), is.null(attr(b,"alleles"))))
		stop("Ambiguous allele encoding in one or both of the objects to be merged.")
	if (attr(a,"alleles") != attr(b,"alleles"))
		stop("Allele encodings don't match -- nonsensical result.")
	
	cols.a <- colnames(a)
	cols.b <- colnames(b)
	if (length(setdiff(cols.a,cols.b)) || length(setdiff(cols.b,cols.a)))
		stop("Number and names of samples don't match.  Try merging genotypes instead of rbind-ing.")
	
	message(paste0("Adding ",nrow(b)," markers to the existing ",nrow(a),"."))
	
	rez <- rbind(unclass(a)[ ,cols.a ], unclass(b)[ ,cols.a ])
	class(rez) <- c("genotypes", class(rez))
	if (!is.null(attr(a, "ped")))
		attr(rez, "ped") <- attr(a, "ped")
	if (!is.null(attr(a, "map")) & !is.null(attr(b, "map")))
		attr(rez, "map") <- rbind( attr(a, "map"), attr(b, "map") )
	attr(rez, "alleles") <- attr(a, "alleles")
	
	if (.has.valid.intensity(a) && .has.valid.intensity(b)) {
		message("Updating intensity matrices...")
		x <- rbind( attr(a, "intensity")$x, attr(b, "intensity")$x )
		y <- rbind( attr(a, "intensity")$y, attr(b, "intensity")$y )
		attr(rez, "intensity") <- list(x = x, y = y)
	}
	
	return(rez)
	
}

#' Merge two \code{genotypes} objects which share markers
#' 
#' @param x a \code{genotypes} object
#' @param y another \code{genotypes} object
#' @param join what sort of join to perform: at present, only \code{"inner"} (intersection of \code{x}'s and \code{y}'s
#' 	marker sets) is supported
#' @param check.alleles logical; if \code{TRUE}, attempt to verify compatibility of marker maps and fix allele or strand swaps
#' @param verbose logical; if \code{TRUE}, show lots of detail about the merge process, especially allele/strand swaps
#' @param ... ignored
#'  
#' @return a \code{genotypes} object containing all samples in \code{x} and \code{y}, at all markers shared
#' 	by \code{x} and \code{y}
#' 	
#' @details Both input objects should have the same allele coding (character vs. numeric, relative vs. absolute)
#' 	before merging.  Sharing of markers is established by the presence of overlapping rownames in the genotype
#' 	matrices of the input, but alleles are NOT checked.  Marker names must match exactly if they are to appear
#' 	in the output.  Intensity matrices, if present, are merged and carried over to the output. They will be marked
#' 	as normalized only if they are marked as such in both input objects.  Vectors of site and sample filters will
#' 	be present in the output only if they are present in both input objects, and have the expected names.  A
#' 	logical OR is applied to the site filters -- ie. a site will be marked as \code{TRUE} (filtered) if it is
#' 	filtered in either or both input datasets.
#' 	
#' @section Warnings:
#' Sample filters currently are *dropped* by the merge operation.
#' 
#' @export
merge.genotypes <- function(x, y, join = c("inner","left"), check.alleles = FALSE, verbose = TRUE, ...) {
	
	if (!inherits(y, "genotypes"))
		stop("Only willing to merge two objects of class 'genotypes.'")
	
	if (length(intersect(colnames(x), colnames(y))))
		stop("Some samples are shared: that merge isn't implemented yet.")
	
	if (mode(x) != mode(y))
		stop("The two objects appear to have different modes (character vs. numeric).")
	
	if (!is.null(attr(x, "alleles"))) {
		if (!is.null(attr(y, "alleles"))) {
			if (attr(y, "alleles") != attr(x, "alleles"))
				warning("Alleles are coded differently in the two datasets; consider fixing that before merging.")
		}
	}
	
	o <- setNames( 1:nrow(x), rownames(x) )
	keep <- intersect(rownames(x), rownames(y))
	if (!length(keep))
		stop("No shared markers; result would be empty.")
	
	message(paste0("Set A has ", nrow(x), " markers x ", ncol(x), " samples."))
	message(paste0("Set B has ", nrow(y), " markers x ", ncol(y), " samples."))
	
	join <- match.arg(join)
	new.o <- keep[ order(o[keep]) ]
	
	if (join == "inner") {
		
		## check for strand/allele swaps
		if (check.alleles) {
			
			if (!(attr(x, "alleles") == "01" && attr(y, "alleles") == "01"))
				stop("In order to perform allele check efficiently, both input datasets should be in",
					 "the '01' numeric encoding.")
				
			x <- x[ new.o, ]
			y <- y[ new.o, ]
			
			message("Checking for allele and strand swaps...")
			unswapped <- .fix.allele.swaps(y, x, verbose = verbose)
			message("Done fixing swaps.")
			x <- unswapped[[2]]
			y <- unswapped[[1]]
			
		}
		
		## keep intersection of marker sets
		rez <- cbind( unclass(x)[ new.o, ],
					  unclass(y)[ new.o, ] )
		if (.has.valid.intensity(x) && .has.valid.intensity(y)) {
			attr(rez, "intensity") <- list( x = cbind(attr(x, "intensity")$x[ new.o,, drop = FALSE ], attr(y, "intensity")$x[ new.o,, drop = FALSE ]),
											y = cbind(attr(x, "intensity")$y[ new.o,, drop = FALSE ], attr(y, "intensity")$y[ new.o,, drop = FALSE ]) )
			attr(rez, "normalized") <- .null.false(attr(x, "normalized")) && .null.false(attr(x, "normalized"))
		}
		if (.has.valid.baflrr(x) && .has.valid.baflrr(y)) {
			attr(rez, "baf") <- cbind( attr(x, "baf")[ new.o, drop = FALSE ],
									   attr(y, "baf")[ new.o, drop = FALSE ] )
			attr(rez, "lrr") <- cbind( attr(x, "lrr")[ new.o, drop = FALSE ],
									   attr(y, "lrr")[ new.o, drop = FALSE ] )
		}
	}
	else
		stop("Not yet implemented: merges other than 'inner join'.")
	
	message(paste0("Merged set has ", nrow(rez), " markers x ", ncol(rez), " samples."))
	
	## add class info
	class(rez) <- c("genotypes", class(rez))
	colnames(rez) <- c(colnames(x), colnames(y))
	attr(rez, "alleles") <- attr(x, "alleles")
	
	## merge markers and family information
	if (!is.null(attr(x, "map")))
		attr(rez, "map") <- attr(x, "map")[ new.o, ]
	if (!is.null(attr(x, "ped")) & !is.null(attr(y, "ped")))
		attr(rez, "ped") <- rbind(attr(x, "ped"), attr(y, "ped"))
	
	## merge filters
	if (!any(is.null(attr(y, "filter.sites")), is.null(attr(x, "filter.sites")))) {
		#fa <- attr(x, "filter.sites")
		#fb <- attr(y, "filter.sites")
		#attr(rez, "filter.sites") <- .merge.filters(fa, fb)
		attr(rez, "filter.sites") <- .init.filters(rownames(rez))
	}
	else {
		attr(rez, "filter.sites") <- .init.filters(rownames(rez))
	}
	if (!any(is.null(attr(y, "filter.samples")), is.null(attr(x, "filter.samples")))) {
		attr(rez, "filter.samples") <- c( attr(x, "filter.samples"), attr(y, "filter.samples") )[ colnames(rez) ]
	}
	
	return(rez)
	
}

## try to reconcile alleles in 'genotypes' objects with identical marker sets
.fix.allele.swaps <- function(a, b, verbose = FALSE, ...) {
	
	rc <- function(x) {
		bases <- c(A="T",C="G",T="A",G="C")
		if (is.na(bases[x]))
			return("N")
		else
			return(bases[x])
	}
	
	aa1 <- toupper(as.character(attr(a, "map")$A1))
	aa2 <- toupper(as.character(attr(a, "map")$A2))
	ba1 <- toupper(as.character(attr(b, "map")$A1))
	ba2 <- toupper(as.character(attr(b, "map")$A2))
	
	ii <- which(aa1 != ba1 | aa2 != ba2)
	message("\tinvestigating ", length(ii), " potential swaps...")
	
	for (i in ii) {
		#cat(rownames(a)[i], "\n")
		if (aa1[i] == ba1[i]) {
			if(aa2[i] == ba2[i]) {
				## matching alleles
				next
			}
			else {
				aa2[i] <- NA
			}
		}
		else if (rc(aa1[i]) == ba1[i]){
			if (rc(aa2[i]) == ba2[i]) {
				## matching allele order; opposite strand
				if (verbose)
					message("\tswapping strand at marker ", rownames(a)[i])
				aa1[i] <- rc(aa1[i])
				aa2[i] <- rc(aa2[i])
			}
			else {
				aa2[i] <- NA
			}
		}
		else if (aa1[i] == ba2[i]) {
			if (aa2[i] == ba1[i]) {
				## different allele order; same strand
				if (verbose)
					message("\tswapping order at marker ", rownames(a)[i])
				tmp <- aa1[i]
				aa1[i] <- aa2[i]
				aa2[i] <- tmp
				a[ i, ] <- abs(2 - a[ i, ]) # swap calls
			}
			else {
				aa2[i] <- NA
			}
		}
		else if (rc(aa1[i]) == ba2[i]) {
			if (rc(aa2[i]) == ba1[i]) {
				## different allele order; opposite strand
				if (verbose)
					message("\tswapping strand and order at marker ", rownames(a)[i])
				tmp <- rc(aa1[i])
				aa1[i] <- rc(aa2[i])
				aa2[i] <- tmp
				a[ i, ] <- abs(2 - a[ i, ]) # swap calls
			}
			else {
				aa2[i] <- NA
			}
		}
		else {
			aa2[i] <- NA
		}
	}

	attr(a, "map")$A1 <- aa1
	attr(a, "map")$A2 <- aa2
	attr(b, "map")$A1 <- ba1
	attr(b, "map")$A2 <- ba2
	
	return(list(a, b))
	
}

## internal helpers for validating the 'genotypes' data structure and its parts

.is.valid.map <- function(mm, ...) {
	
	pass <- is.data.frame(mm)
	pass <- pass && all(colnames(mm)[1:4] == c("chr","marker","cM","pos"))
	if ("marker" %in% colnames(mm))
		pass <- all(rownames(mm) == as.character(mm$marker))
	else
		pass <- FALSE
	
	return(pass)
	
}

.has.valid.map <- function(gty, ...) {
	
	map <- attr(gty, "map")
	pass <- FALSE
	
	if (!is.null(map))
		pass <- .is.valid.map(map) && (nrow(map) == nrow(gty))
	if (is.na(pass))
		pass <- FALSE
	
	return(pass)
	
}

.has.valid.ped <- function(gty, ...) {
	
	ped <- attr(gty, "ped")
	pass <- FALSE
	
	if (!is.null(ped))
		pass <- .is.valid.ped(ped) && (nrow(ped) == ncol(gty))
	if (is.na(pass))
		pass <- FALSE
	
	return(pass)
	
}

.is.valid.ped <- function(ped, ...) {
	
	pass <- is.data.frame(ped)
	pass <- pass && all(colnames(ped)[1:6] == c("fid","iid","mom","dad","sex","pheno"))
	if ("iid" %in% colnames(ped))
		pass <- all(rownames(ped) == as.character(ped$iid))
	else
		pass <- FALSE
	
	return(pass)
	
}

.has.valid.baflrr <- function(gty, ...) {
	
	pass <- FALSE
	if (!is.null(attr(gty, "baf")))
		if (all( dim(attr(gty, "baf")) == dim(attr(gty, "lrr")),
				 dim(gty) == dim(attr(gty, "baf"))))
			pass <- TRUE
	
	return(pass)
	
}

.null.false <- function(x) {
	
	if (is.null(x))
		FALSE
	else x
	
}

.has.valid.intensity <- function(gty, ...) {
	
	pass <- FALSE
	if (!is.null(attr(gty, "intensity")))
		if (is.list(attr(gty, "intensity")) && length(attr(gty, "intensity")) == 2)
			if ( all(dim(attr(gty, "intensity")$x) == dim(attr(gty, "intensity")$y)) &&
				 	all(dim(attr(gty, "intensity")$x) == dim(gty)) )
				pass <- TRUE
	
	return(pass)
	
}

.fix.sex <- function(x) {
	
	if (is.character(x) || is.factor(x))
		ifelse(grepl("^[fF]", x), 2, ifelse(grepl("^[mM]", x), 1, 0))
	else
		return(x)
	
}

#' Check the integrity of a \code{genotypes} object
#'
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' 
#' @return Logical scalar indicating whether object passes or fails integrity checks.
#'
#' @details A valid genotypes object must have the following: a genotypes matrix with row and column names;
#' 	a valid marker map (chr,marker,cM,pos) whose rownames match the rownames of the genotypes matrix;
#' 	filter flags for sites and samples, with appropriate dimensions; and an indicator for allele coding.
#' 	
#' 	If sample metadata is present, it must be a dataframe with at least a column 'iid' which is identical
#' 	to its rownames.  A column 'sex' is strongly encouraged but its absence will not cause validation to fail.
#' 	
#' 	If intensity data is present, the x- and y-intensity matrices will be checked for dimensions and names
#' 	in manner analogous to the checks on the genotypes matrix.  If LRR and BAF have been computed, they
#' 	will be checked also.
#'
#' @aliases validate
#' @export
validate.genotypes <- function(gty, ...) {
	
	## check genotypes matrix
	pass <- TRUE
	pass <- pass && is.matrix(gty)
	if (is.numeric(gty)) {
		if (!is.null(attr(gty, "alleles"))) {
			pass <- pass && (attr(gty, "alleles") %in% c("01","relative"))
		}
		else {
			pass <- FALSE
		}
	}
	else if (is.character(gty)) {
		if (!is.null(attr(gty, "alleles"))) {
			pass <- pass && (attr(gty, "alleles") %in% c("native"))
		}
		else {
			pass <- FALSE
		}
	}
	else {
		pass <- FALSE
	}
	if (!pass) {
		message("Genotype matrix is not a matrix, or allele coding is unrecognized.")
		return(pass)
	}
	
	## check marker map
	pass <- pass && .has.valid.map(gty)
	pass <- all(rownames(attr(gty, "map")) == rownames(gty))
	if (!pass) {
		message("Marker map is absent, malformed, or does not match genotypes matrix.")
		return(pass)
	}
	
	## check sample metadata
	if (.has.valid.ped(gty)) {
		pass <- pass && all( rownames(attr(gty, "ped")) == attr(gty, "ped")$iid,
							 rownames(attr(gty, "ped")) == colnames(gty) )
	}
	if (!pass) {
		message(paste0("Sample metadata does not match genotypes matrix.  Check that it has",
					   "a column 'iid' which matches its rownames, and that these rownames",
					   "match column names of genotypes matrix."))
		return(pass)
	}
	
	## check intensity
	if (!is.null(attr(gty, "intensity"))) {
		pass <- pass && .has.valid.intensity(gty)
		pass <- pass && all( rownames(attr(gty, "intenxity")$x) == rownames(attr(gty, "intenxity")$y),
							 colnames(attr(gty, "intenxity")$x) == colnames(attr(gty, "intenxity")$y),
							 colnames(attr(gty, "intenxity")$x) == colnames(gty),
							 rownames(attr(gty, "intenxity")$x) == rownames(gty) )
	}
	if (!pass) {
		message("Intensity matrices are malformed or do not match genotypes matrix.")
		return(pass)
	}
	
	## check BAF/LRR
	if (!is.null(attr(gty, "baf"))) {
		if (!is.null(attr(gty, "lrr"))) {
			pass <- pass && .has.valid.baflrr(gty)
			pass <- pass && all( rownames(attr(gty, "baf")) == rownames(attr(gty, "lrr")),
								 colnames(attr(gty, "baf")) == colnames(attr(gty, "baf")),
								 colnames(attr(gty, "baf")) == colnames(gty),
								 rownames(attr(gty, "baf")) == rownames(gty) )
		}
		else {
			pass <- FALSE
		}
	}
	if (!pass) {
		message(paste0("BAF and LRR matrices are malformed or do not match genotypes matrix.  Note: if BAF is present,",
					   "LRR must be also."))
		return(pass)
	}
	
	## check filters
	if (!is.null(attr(gty, "filter.sites"))) {
		pass <- pass && all( length(attr(gty, "filter.sites")) == nrow(gty),
							 !is.na(attr(gty, "filter.sites")) )
	}
	if (!is.null(attr(gty, "filter.samples"))) {
		pass <- pass && all( length(attr(gty, "filter.samples")) == ncol(gty),
							 !is.na(attr(gty, "filter.samples")) )
	}
	if (!pass) {
		message(paste("Filters are malformed.  Sites filter should have length equal to number of rows",
					  "of genotypes matrix, and samples filter, to number of columns.  Names are not required."))
		return(pass)
	}
	
	return(pass)
	
}
#' @export
validate <- function(gty, ...) UseMethod("validate")

#' Strip intensity matrices from a \code{genotypes} object
#' 
#' @param gty a \code{genotypes} object
#' @param ... ignored
#'
#' @return a copy of \code{gty} without intensity data (but with any QC summaries still present)
#' 
#' @export drop.intensity
drop.intensity <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	for (a in c("intensity","baf","lrr"))
		attr(gty, a) <- NULL
	
	return(gty)
	
}

#' Apply a function over a \code{genotypes} object, by sample or marker groups
#'
#' An homage to base \code{R}'s \code{apply} family of functions.
#'
#' @param gty an object of class \code{genotypes}
#' @param margin dimension of the array over which to apply \code{expr} (1=rows, 2=columns)
#' @param expr an expression defining subsets of \code{gty}, either by markers (\code{margin = 1})
#' 	or samples (\code{margin = 1})
#' @param fn the function to apply over the subsets defined by \code{expr}
#' @param strip logical; if \code{TRUE}, strip metadata before evaluating \code{fn()}
#' @param ... other arguments passed-through to \code{fn()}
#' 
#' @return a list with one element per value of the grouping expression (\code{expr}); the value of that
#' 	element is the return value of \code{fn()} applied to the subset of \code{gty} defined
#' 	by \code{expr}
#' 
#' @details The grouping expression \code{expr} is evaluated in the parent environment by default. Use the
#' 	dot function to get lazy evaluation like with \code{subset()} (see examples).
#' 
#' 	When \code{fn()} will involve heavy use of array slicing, use \code{strip} to avoid having to
#' 	slice sample/marker metadata repeatedly.
#' 
#' @examples
#' # lazy evaluation 1: apply over chromosomes
#' \dontrun{
#' genoapply(x, 1, .(chr), myfunction, extra.arg = TRUE)
#' }
# lazy evaluation 2: apply over families
#' \dontrun{
#' genoapply(x, 2, .(fid), myfunction, extra.arg = TRUE)
#' }
#' 
#' @export
genoapply <- function(gty, margin = c(1,2), expr = 1, fn = NULL, strip = FALSE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Only willing to subset an object of class 'genotypes.'")
	
	#e <- eval(expr)
	e <- expr
	nn <- as.quoted(e)
	if (margin == 2) {
		if (!is.null(attr(gty, "ped")))
			r <- lapply(e, eval, envir = attr(gty, "ped"), enclos = parent.frame())
		else
			r <- lapply(e, eval, envir = parent.frame())
	}
	else {
		if (!is.null(attr(gty, "map")))
			r <- lapply(e, eval, envir = attr(gty, "map"), enclos = parent.frame())
		else
			r <- lapply(e, eval, envir = parent.frame())
	}
	
	
	if (strip)
		gty <- bless(.copy.matrix.noattr(gty))
	
	r <- lapply(r, factor)
	splitter <- data.frame(expand.grid(lapply(r, levels)))
	colnames(splitter) <- names(nn)
	if (margin == 2) {
		vals <- split(seq_len(ncol(gty)), r)
		rez <- lapply(vals, function(v) {
			x <- gty[ ,v,drop = FALSE ]
			fn(x, ...)
		})
	}
	else {
		vals <- split(seq_len(nrow(gty)), r)
		rez <- lapply(vals, function(v) {
			x <- gty[ v,,drop = FALSE ]
			fn(x, ...)
		})
	}
	
	names(rez) <- names(vals)
	nulls <- sapply(rez, is.null)
	attr(rez, "split_labels") <- splitter
	return(rez)
	
}

#' Switch between character and numeric representations of genotype calls
#' 
#' @param gty a \code{genotypes} object
#' @param mode conversion mode (see Details)
#' @param allowed a vector of allowable alleles for character-encoded genotypes
#' @param alleles a 2-column character matrix containing the reference and alternate allele
#' 	at each marker in the input object; if null, obtained from marker map
#' @param ... ignored
#' 
#' @return a copy of \code{gty} with alleles converted to the requested encoding
#' 
#' @details Genotypes on Illumina Infinium arrays are assumed to correspond to biallelic SNPs.
#' 	Although the BeadStudio software reports calls in character form, the numeric representation
#' 	(coded 0/1/2 for homozygous reference / heterozygous / homozygous alternate) is computationally
#' 	much more convenient.  This function performs the character-to-numeric conversion or its inverse.
#' 	
#' 	When \code{mode == "pass"} (the default); the input is returned unchanged.
#' 	When \code{mode == "01"}, conversion from character to numeric is performed.  If \code{alleles}
#' 	is supplied or reference alleles are provided in the marker map, the coding will
#' 	be with respect to the reference (A1) or alternate (A2) alleles.  If reference alleles
#' 	are not provided, the recoding falls back to \code{mode == "relative"}.
#' 	When \code{mode == "relative"}, conversion from character to numeric is performed,
#' 	and the coding is with respect to the minor allele in the current dataset.  This can
#' 	be a problem if the dataset is small, and obviously hampers comparisons to other datasets.
#' 	Use this option with caution. (On the other hand, the concept of a minor allele may be
#' 	very useful in the context of true population samples.)
#' 	When \code{mode == "native"}, conversion from numeric back to character is attempted.
#' 	If the input had mode "relative", an error occurs -- that conversion is too error-prone.
#' 	Otherwise the converison uses the supplied reference alleles.
#' 
#' @aliases recode
#' @export
recode.genotypes <- function(gty, mode = c("pass","01","native","relative"),
							 allowed = c("A","C","G","T","H"), alleles = NULL, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	coding <- attr(gty, "alleles")
	
	## do the genotypes come with a map and reference alleles?
	alleles <- NULL
	if (.has.valid.map(gty))
		if (ncol(attr(gty, "map")) >= 6)
			alleles <- as.matrix(attr(gty, "map")[ ,5:6 ])
	
	## a suite of recoding functions which operate on (marker-wise) vectors of genotype calls
	top1 <- function(x, nbins) {
		tbl <- tabulate(x, nbins = nbins)
		a <- which.max(tbl)
		return(a)
	}
	# recode 0=major allele, 1=het, 2=minor allele
	.recode.numeric.by.freq <- function(calls, alleles = NULL) {
		
		calls.f <- factor( calls, levels = allowed[1:4] )
		maj <- allowed[ top1(calls.f, 4) ]
		new.calls <- rep(NA, length(calls))
		new.calls[ calls == "H" ] <- 1
		new.calls[ calls == maj ] <- 0
		new.calls[ is.na(new.calls) & !is.na(calls.f) ] <- 2
		return(new.calls)
		
	}
	
	.recode.numeric.by.freq.from.numeric <- function(calls, alleles = NULL) {
		
		maj <- which.max(tabulate(calls+1, nbins = 3)[ c(1,3) ])-1
		new.calls <- rep(0, length(calls))
		#new.calls[ calls == maj ] <- 0
		new.calls[ calls == 1 ] <- 1
		new.calls[ abs(calls-maj) == 2 ] <- 2
		new.calls[ is.na(calls) ] <- NA
		return(new.calls)
		
	}
	
	# recode 0=A1, 1=het, 2=A2
	.recode.numeric.by.ref <- function(calls, alleles = NULL) {
		
		new.calls <- rep(NA, length(calls))
		new.calls[ calls == alleles[1] ] <- 0
		new.calls[ calls == alleles[2] ] <- 2
		new.calls[ calls == "H" ] <- 1
		return(new.calls)
		
	}
	
	# recode same as above, only backwards
	.recode.character.by.ref <- function(calls, alleles = NULL) {
		
		new.calls <- rep("N", length(calls))
		new.calls[ calls == 0 ] <- as.character(alleles[1])
		new.calls[ calls == 2 ] <- as.character(alleles[2])
		new.calls[ calls == 1 ] <- "H"
		return(new.calls)
		
	}
	
	.dont.convert <- function(calls, alleles = NULL) identity(calls)
	
	mode <- match.arg(mode)
	recode.as <- FALSE
	if (mode == "pass") {
		recode.as <- coding
		converter <- .dont.convert
	}
	else if (mode == "01") {
		if ((is.null(coding) || coding == "01") && is.numeric(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			recode.as <- "01"
			converter <- .dont.convert
		}
		else if (!is.null(alleles)) {
			converter <- .recode.numeric.by.ref
			recode.as <- "01"
			message("Recoding to 0/1/2 using reference alleles.")
		}
		else {
			recode.as <- "relative"
			message("Recoding to 0/1/2 using empirical frequencies.")
		}
	}
	else if (mode == "relative") {
		if (is.numeric(gty) && (coding == "relative" || is.null(coding))) {
			message("Nothing to do; genotypes already in requested coding.")
			recode.as <- "relative"
			converter <- identity
		}
		else {
			if (is.character(gty))
				converter <- .recode.numeric.by.freq
			else
				converter <- .recode.numeric.by.freq.from.numeric
			recode.as <- "relative"
			message("Recoding to 0/1/2 using empirical frequencies.")
		}
	}
	else if (mode == "native") {
		if (is.character(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			converter <- identity
		}
		else if (!is.null(alleles) && coding != "relative") {
				message("Recoding to character using reference alleles.")
				converter <- .recode.character.by.ref
				recode.as <- "native"
		}
		else
			stop("Can only convert genotypes numeric->character given some reference alleles.")
	}
	
	.gty <- .copy.matrix.noattr(gty)
	rez <- matrix(NA, nrow = nrow(.gty), ncol = ncol(.gty))
	colnames(rez) <- colnames(.gty)
	rownames(rez) <- rownames(.gty)
	for (i in seq_len(nrow(.gty))) {
		rez[ i, ] <- converter(.gty[i,], alleles = alleles[i,])
	}
	
	if (.has.valid.map(gty))
		attr(rez, "map") <- attr(gty, "map")
	if (.has.valid.ped(gty))
		attr(rez, "ped") <- attr(gty, "ped")
	
	for (a in c("intensity","normalized","filter.sites","filter.samples","baf","lrr","qc")) {
		if (!is.null(attr(gty, a)))
			attr(rez, a) <- attr(gty, a)
	}
	class(rez) <- c("genotypes", class(rez))
	attr(rez, "alleles") <- recode.as
	return(rez)
	
}
#' @export
recode <- function(gty, ...) UseMethod("recode")

#' Recode genotypes against genotypes of a parent
#' 
#' @param gty a \code{genotypes} object
#' @param parent a numeric vector of parental genotypes; a scalar pointing to a reference sample in
#' 	the input (\code{gty}); or another \code{genotypes} with parental genotypes in the first column
#' @param ... ignored
#' 
#' @return a recoded \code{genotypes} object
#' 
#' @export recode.to.parent
recode.to.parent <- function(gty, parent, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (!inherits(parent, "genotypes")) {
		refs <- as.vector(gty[ ,parent ])
	}
	else if (inherits(parent, "genotypes")) {
		prn <- rownames(parent)
		if (!setequal(prn, rownames(gty))) {
			stop("Markers in target object and parental genotypes don't match.")
		}
		parent <- parent[ rownames(gty), ]
		refs <- as.vector(parent[,1])
	}
	else if (is.numeric(parent))
		refs <- as.vector(parent)
	
	if (!(is.numeric(gty) && is.numeric(refs)))
		stop("All genotypes must be in numeric encoding.")
	if (length(refs) != nrow(gty))
		stop("Dimensions of target object and parental genotypes don't match.")
	
	converter <- function(x) {
		2-abs(refs-x)
	}
	
	.gty <- .copy.matrix.noattr(gty)
	rez <- matrix(NA, nrow = nrow(.gty), ncol = ncol(.gty),
				  dimnames = list(rownames(.gty), colnames(.gty)))

	for (i in seq_len(ncol(.gty))) {
		rez[ ,i ] <- converter(.gty[,i])
	}
	
	if (.has.valid.map(gty))
		attr(rez, "map") <- attr(gty, "map")
	if (.has.valid.ped(gty))
		attr(rez, "ped") <- attr(gty, "ped")
	
	for (a in c("intensity","normalized","filter.sites","filter.samples","baf","lrr","qc")) {
		if (!is.null(attr(gty, a)))
			attr(rez, a) <- attr(gty, a)
	}
	class(rez) <- c("genotypes", class(rez))
	attr(rez, "alleles") <- "parent"
	return(rez)
	
}

#' Swap out the marker map for a different one
#' 
#' @param gty a \code{genotypes} object
#' @param newmap a valid marker map (especially: rownames same as marker names)
#' @param ... ignored
#' 		
#' @return a copy of \code{gty} with old marker map replaced by the new one
#' 
#' @details Replacement of the marker map relies on marker names: only markers with which are shared
#' 	by the old and new maps are retained in the output.  (As such, this function is also a backdoor
#' 	subsetting function.)  Intensity matrices, BAF/LRR matrices, and site filters are adjusted accordingly.
#' 	A validity check is performed before the function returns to flag problems.
#' 
#' @export replace.map
replace.map <- function(gty, newmap, ...) {
	
	if (!inherits(gty, "genotypes") || !.has.valid.map(gty))
		stop("Please supply a genotypes object with map as an attribute.")
	
	if (!.is.valid.map(newmap))
		stop(paste("New map is not a valid marker map.  It should have (at least) the following",
				   "columns, in order: chr, marker, cM, pos."))
	
	map <- attr(gty, "map")
	map$order <- 1:nrow(map)
	rownames(map) <- as.character(map$marker)
	message(paste("Starting with", nrow(map), "markers..."))
	
	if (!is.null(rownames(newmap))) {
		## extract only markers which are (1) present in current map; (2) accounted for in new map
		m <- match(as.character(map$marker), rownames(newmap), nomatch = 0)
		gty <- gty[ which(m > 0), ]
		rownames(gty) <- newmap[ m[m > 0],"marker" ]
		attr(gty, "map") <- newmap[ m[m > 0], ]
	}
	else {
		attr(gty, "map") <- newmap
	}
	
	## sort the result
	o <- order( attr(gty, "map")$chr, attr(gty, "map")$pos, attr(gty, "map")$marker )
	map <- attr(gty, "map")
	gty <- gty[ o, ]
	
	## sort intensities too, if present
	if (.has.valid.intensity(gty)) {
		x.new <- attr(gty, "intensity")$x[ which(m > 0), ]
		y.new <- attr(gty, "intensity")$y[ which(m > 0), ]
		x.new <- x.new[ o, ]
		y.new <- y.new[ o, ]
		attr(gty, "intensity") <- list(x = x.new, y = y.new)
	}
	if (.has.valid.baflrr(gty)) {
		baf.new <- attr(gty, "baf")[ which(m > 0), ]
		lrr.new <- attr(gty, "lrr")[ which(m > 0), ]
		attr(gty, "baf") <- baf.new[ o, ]
		attr(gty, "lrr") <- lrr.new[ o, ]
	}
	if (!is.null(attr(gty, "filter.sites"))) {
		attr(gty, "filter.sites") <- attr(gty, "filter.sites")[o]
	}
	
	cols <- colnames(map)
	required <- c("chr","marker","cM","pos")
	attr(gty, "map") <- map[ o,c(required, setdiff(cols, required)) ]
	
	message(paste("...and ending with", nrow(attr(gty, "map")), "markers."))
	if (!validate.genotypes(gty))
		stop("Oops -- replacing the map corrupted the genotypes object.")
	
	return(gty)
	
}

## convert between styles of chromosome names (I prefer UCSC)
convert.names <- function(chrs, to = c("ncbi","ensembl","ucsc","plink"), muga.style = FALSE, ...) {
	
	ucsc <- paste0("chr", c(1:19, "X","Y","M"))
	ncbi <- c(1:19,"X","Y","MT")
	
	to <- match.arg(to)
	converter <- function(x) { identity(x) }
	if (to == "ncbi" | to == "ensembl") {
		converter <- function(x) setNames(ncbi, ucsc)[x]
	}
	else if (to == "ucsc") {
		converter <- function(x) setNames(ucsc, ncbi)[x]
	}
	else if (to == "plink") {
		#message("converting to plink")
		converter <- function(x) {
			nn <- gsub("^chr","", x)
			nn <- gsub("^M$","MT", nn)
			return(nn)
		}
	}
	
	nn <- converter( as.character(chrs) )
	if (muga.style)
		nn <- gsub("^MT$","M", chrs)
	
	return(nn)
	
}

## get rid of unfriendly characters in sample or marker names
clean.names <- function(x, space = "", ...) {
	
	x <- gsub("^\\s+|\\s+$", "", x)
	x <- make.unique(x)
	x <- gsub(" ", space, x)
	return(x)
	
}


#' Fix alleles in marker map using those in observed genotypes
#' 
#' @param gty a \code{genotypes} object
#'
#' @return A new \code{genotypes} object with corrected alleles, whose marker map has two
#' 	additional columns: \code{swapped} indicates if alleles were swapped, and \code{unresolved}
#' 	indicates that observed and expected alleles don't match, but can't be reconciled with
#' 	a simple strand swap.
#' 	
#' @details Because this method relies on information in observed genotypes, it is advisable
#' 	to use it only with "many" samples, where "many" is enough to capture all alleles of interest.
#'  
#' @export
fix.strand.swaps <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes") || !.has.valid.map(gty))
		stop("Please supply a genotypes object with map as an attribute.")
	
	if (attr(gty, "alleles") != "native")
		stop("Genotypes should be in 'native' character encoding.")
	
	revcomp = c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
	
	geno <- unclass(gty)
	map <- attr(gty, "map")
	
	message("Checking for strand swaps...")
	failures <- rep(FALSE, nrow(map))
	swaps <- rep(FALSE, nrow(map))
	for (i in seq_len(nrow(geno))) {
		
		obs <- setdiff(unique(geno[ i, ]), c("N","H"))
		a1 <- map$A1[i]
		a2 <- map$A2[i]
		expect <- c(a1,a2)
		expect.comp <- revcomp[expect]
		if (!all(obs %in% expect)) {
			if (all(obs %in% expect.comp)) {
				map$A1[i] <- expect.comp[1]
				map$A2[i] <- expect.comp[2]
				swaps[i] <- TRUE
			}
			else {
				failures[i] <- TRUE
			}
		}
		else {
			failures[i] <- FALSE
		}
		
	}
	
	message("Detected ", sum(swaps)," strand swaps.")
	message("Couldn't resolve ", sum(failures)," conflicts.")
	
	map$swapped <- swaps
	map$uresolved <- failures
	
	attr(gty, "map") <- map
	return(gty)
	
}

## geno.R
## utility functions for handling genotype matrices

## The 'genotypes' class is a just a matrix (sites x samples) with row and column names.
## Attributes include:
##	* 'map' -- marker metadata in PLINK format (chr, marker, cM, pos, A1, A2, ...)
##	* 'ped' -- pedigree/sample metadata in PLINK format (individual ID, family ID, mom ID, dad ID, sex, phenotype, ...)
##	* 'inensity' -- list(x = [X-intensities], y = [y-intensities])
##	* 'normalized' -- have intensities been normalized?
##	* 'filter.sites' -- homage to the FILTER field in VCF format, a flag for suppresing sites (rows) in downstream analyses
##	* 'filter.samples' -- same as above, but along other dimension (columns)
##	* 'alleles' -- manner in which alleles are encoded: "native" (ACTGHN), "01" (allele dosage wrt ALT allele), "relative" (allele dosage wrt MINOR allele)
## NB: I use missing values (NAs/NaNs) for no-calls, in order to take advantage of R's behaviors on missing data.

## overload indexing operator to make genotypes object behave like normal matrix while preserving attributes
`[.genotypes` <- function(x, i = TRUE, j = TRUE, drop = FALSE, ...) {
	
	r <- NextMethod("[")
	if (!is.null(attr(x, "map"))) {
		attr(r, "map") <- attr(x, "map")[ i, ]
	}
	if (!is.null(attr(x, "ped"))) {
		attr(r, "ped") <- attr(x, "ped")[ j, ]
		rownames(attr(r, "ped")) <- as.character(attr(r, "ped")$iid)
	}
	if (.has.valid.intensity(x)) {
		x.new <- attr(x, "intensity")$x[ i,j, drop = FALSE ]
		y.new <- attr(x, "intensity")$y[ i,j, drop = FALSE ]
		attr(r, "intensity") <- list(x = x.new, y = y.new)
		if (!is.null(attr(x, "normalized")))
			attr(r, "normalized") <- attr(x, "normalized")
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

`$.genotypes` <- function(x, expr, ...) {
	
	attributes(x)[[expr]]
	
}


## set some S3 generics for useful accessor functions

markers <- function(x) UseMethod("markers")
markers.genotypes <- function(gty, ...) {
	attr(gty, "map")
}

samples <- function(x) UseMethod("samples")
samples.genotypes <- function(gty, ...) {
	if (!is.null(attr(gty, "ped")))
		attr(gty, "ped")
	else
		colnames(gty)
}

filters <- function(x) UseMethod("filters")
filters.genotypes <- function(gty, ...) {
	get.filters(gty, ...)
}

intensity <- function(x) UseMethod("intensity")
intensity.genotypes <- function(gty, ...) {
	attr(gty, "intensity")
}

## add a method to overload subset() on a genotypes object which
## preserves attributes (markers, samples, intensity...)
subset.genotypes <- function(gty, expr, by = c("markers","samples"), ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Only willing to subset an object of class 'genotypes.'")
	
	by <- match.arg(by)
	e <- substitute(expr)
	if (by == "markers") {
		if (!.has.valid.map(gty))
			stop("Can't subset genotypes by marker without a marker map.")
		r <- eval(e, attr(gty, "map"), parent.frame())
		r <- r & !is.na(r)
		return( gty[ r, ] )
	}
	else if (by == "samples") {
		if (is.null(attr(gty, "ped")))
			stop("Can't subset genotypes by samples without sample information.")
		r <- eval(e, attr(gty, "ped"), parent.frame())
		r <- r & !is.na(r)
		return( gty[ ,r ] )
	}
	else {
		stop()
	}
	
}

## overload cbind() to join not only genotype matrix but also sample metadata
## NB: intensities and filters are dropped at this point
cbind.genotypes <- function(a, b, ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to bind two objects of class 'genotypes.'")
	
	if (nrow(a) != nrow(b) | any(rownames(a) != rownames(b)))
		stop("Number and names of markers don't match.  Try merging genotypes instead of cbind-ing.")
	
	message(paste0("Adding ",ncol(b)," individuals to the existing ",ncol(a),"."))
	
	rez <- cbind(unclass(a), unclass(b))
	class(rez) <- c("genotypes", class(rez))
	if (!is.null(attr(a, "map")))
		attr(rez, "map") <- attr(a, "map")
	if (!is.null(attr(a, "ped")))
		attr(rez, "ped") <- rbind( attr(a, "ped"), attr(b, "ped") )
	
	return(rez)
	
}

## overload rbind() to join not only genotype matrix but also marker metadata
## NB: intensities and filters are dropped at this point
rbind.genotypes <- function(a, b, ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to bind two objects of class 'genotypes.'")
	
	cols.a <- colnames(a)
	cols.b <- colnames(b)
	if (length(setdiff(a,b)) || length(setdiff(b,a)))
		stop("Number and names of markers don't match.  Try merging genotypes instead of cbind-ing.")
	
	message(paste0("Adding ",nrow(b)," markers to the existing ",nrow(a),"."))
	
	rez <- rbind(unclass(a)[ ,cols.a ], unclass(b)[ ,cols.a ])
	class(rez) <- c("genotypes", class(rez))
	if (!is.null(attr(a, "ped")))
		attr(rez, "ped") <- attr(a, "ped")
	if (!is.null(attr(a, "map")) & !is.null(attr(b, "map")))
		attr(rez, "map") <- rbind( attr(a, "map"), attr(b, "map") )
	
	return(rez)
	
}

## add a method to overload merge(), again keeping attributes
merge.genotypes <- function(a, b, join = c("inner","left"), ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to merge two objects of class 'genotypes.'")
	
	if (length(intersect(colnames(a), colnames(b))))
		stop("Some samples are shared: that merge isn't implemented yet.")
	
	if (mode(a) != mode(b))
		stop("The two objects appear to have different modes (character vs. numeric).")
	
	if (!is.null(attr(a, "alleles"))) {
		if (!is.null(attr(b, "alleles"))) {
			if (attr(b, "alleles") != attr(a, "alleles"))
				warning("Alleles are coded differently in the two datasets; consider fixing that before merging.")
		}
	}
	
	o <- setNames( 1:nrow(a), rownames(a) )
	keep <- intersect(rownames(a), rownames(b))
	if (!length(keep))
		stop("No shared markers; result would be empty.")
	
	message(paste0("Set A has ", nrow(a), " markers x ", ncol(a), " samples."))
	message(paste0("Set B has ", nrow(b), " markers x ", ncol(b), " samples."))
	
	join <- match.arg(join)
	new.o <- keep[ order(o[keep]) ]
	if (join == "inner") {
		## keep intersection of marker sets
		rez <- cbind( unclass(a)[ new.o, ],
					  unclass(b)[ new.o, ] )
		if (.has.valid.intensity(a) && .has.valid.intensity(b)) {
			attr(rez, "intensity") <- list( x = cbind(attr(a, "intensity")$x[ new.o, ], attr(b, "intensity")$x[ new.o, ]),
											y = cbind(attr(a, "intensity")$y[ new.o, ], attr(b, "intensity")$y[ new.o, ]) )
			attr(rez, "normalized") <- .null.false(attr(a, "normalized")) && .null.false(attr(a, "normalized"))
		}
	}
	else
		stop("Not yet implemented: merges other than 'inner join'.")
	
	message(paste0("Merged set has ", nrow(rez), " markers x ", ncol(rez), " samples."))
	
	## add class info
	class(rez) <- c("genotypes", class(rez))
	
	## merge markers and family information
	if (!is.null(attr(a, "map")))
		attr(rez, "map") <- attr(a, "map")[ keep[ order(o[keep]) ], ]
	if (!is.null(attr(a, "ped")) & !is.null(attr(b, "ped")))
		attr(rez, "ped") <- rbind(attr(a, "ped"), attr(b, "ped"))
	
	return(rez)
	
}

## internal helpers for validating the 'genotypes' data structure and its parts

.is.valid.map <- function(map, ...) {
	return( all(colnames(map)[1:4] == c("chr","marker","cM","pos")) )
}

.has.valid.map <- function(gty, ...) {
	
	map <- attr(gty, "map")
	rez <- FALSE
	
	if (!is.null(map))
		rez <- .is.valid.map(map)
	if (is.na(rez))
		rez <- FALSE
	
	return(rez)
	
}

.has.valid.ped <- function(gty, ...) {
	
	## TODO
	return( !is.null(attr(gty, "ped")) )
	
}

.null.false <- function(x) {
	
	if (is.null(x))
		FALSE
	else x
	
}

.has.valid.intensity <- function(gty, ...) {
	
	rez <- FALSE
	if (!is.null(attr(gty, "intensity")))
		if (is.list(attr(gty, "intensity")) && length(attr(gty, "intensity")) == 2)
			if ( all(dim(attr(gty, "intensity")$x) == dim(attr(gty, "intensity")$y)) &&
				 	all(dim(attr(gty, "intensity")$x) == dim(gty)) )
				rez <- TRUE
	
	return(rez)
	
}

## apply a function over samples in a genotype matrix, by sample groups
gapply <- function(gty, expr, fn = NULL, unclass = FALSE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Only willing to subset an object of class 'genotypes.'")
	
	e <- substitute(expr)
	if (!is.null(attr(gty, "ped")))
		r <- eval(e, attr(gty, "ped"), parent.frame())
	else
		r <- eval(e)
	
	if (unclass)
		gty <- unclass(gty)
	
	vals <- unique(r)
	rez <- lapply(vals, function(v) {
		fn(gty[ ,(r == v), drop = FALSE ], ...)
	})
	names(rez) <- vals
	return(rez)
	
}

maf <- function(x, ...) {
	rowMeans(x, na.rm = TRUE)/2
}

consensus <- function(gty, nas.allowed = 0.0, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	rez <- apply(gty, 1, function(x) which.max(tabulate(x+1, nbins = 3))-1)
	nas <- apply(gty, 1, function(x) sum(is.na(x))/length(x) > nas.allowed)
	rez[nas] <- NA
	
	return(rez)
	#return(as.matrix(rez))
	
}

is.segregating <- function(x, ...) {
	
	if (is.character(x)) {
		x <- factor(x)
		if (any(!is.na(x)) | any(x != "N")) {
			flag <- sum(tabulate(x, nbins = 3) > 0) > 1
			return(flag | any(x == "H", na.rm = TRUE))
		}
		else
			return(FALSE)
	}
	else {
		if (any(!is.na(x))) {
			flag <- sum(tabulate(x+1, nbins = 3) > 0) > 1
			return(flag | any(x == 1, na.rm = TRUE))
		}
		else
			return(FALSE)
	}
	
}

segregating <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	keep <- apply(gty, 1, is.segregating)
	return(gty[ keep, ])
	
}

is.fixed.diff <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (!is.numeric(gty) || ncol(gty) != 2)
		stop("This function is for numeric genotypes for pairs of samples only.")
	
	gty <- unclass(gty)
	#print(head(gty))
	attr(gty, "map") <- NULL
	
	diffs <- (gty[ ,1 ] != gty[ ,2 ])
	ns <- is.na(rowSums(gty))
	hs <- (gty[ ,1 ] == 1) | (gty[ ,2 ] ==1)
	
	return(diffs & !ns & !hs)
	
}

heterozygosity <- function(x, het.char = "H", ...) {
	
	if (is.factor(x) || is.character(x))
		sum(x == het.char, na.rm = TRUE)/sum(!is.na(x))
	else
		sum(x == 1, na.rm = TRUE)/sum(!is.na(x))
	
}

## swap genotypes betewen different coding schemes
## "01": 0/1/2 copies of non-reference allele
## "relative": 0/1/2 copies of minor allele, wrt this dataset
## "native": character genotypes resperenting base calls (ACGT) or H for hets
recode.genotypes <- function(gty, mode = c("pass","01","native","relative"),
							 allowed = c("A","C","G","T","H"), alleles = NULL, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (is.null(attr(gty, "alleles"))) {
		attr(gty, "alleles") <- FALSE
	}
	
	## do the genotypes come with a map and reference alleles?
	alleles <- NULL
	if (.has.valid.map(gty))
		if (ncol(attr(gty, "map")) >= 6)
			alleles <- as.matrix(attr(gty, "map")[ ,5:6 ])
	
	## a suite of recoding functions which operate on (marker-wise) vectors of genotype calls
	
	# recode 0=major allele, 1=het, 2=minor allele
	.recode.numeric.by.freq <- function(calls, alleles = NULL) {
		
		calls <- factor( as.character(calls), levels = allowed )
		maj <- allowed[ top1(calls[ calls != "H" ]) ]
		new.calls <- rep(NA, length(calls))
		new.calls[ calls == "H" ] <- 1
		new.calls[ calls == maj ] <- 0
		new.calls[ is.na(new.calls) & !is.na(calls) ] <- 2
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
		
		new.calls <- rep(NA, length(calls))
		new.calls[ calls == 0 ] <- as.character(alleles[1])
		new.calls[ calls == 2 ] <- as.character(alleles[2])
		new.calls[ calls == 1 ] <- "H"
		return(new.calls)
		
	}
	
	mode <- match.arg(mode)
	coding <- attr(gty, "alleles")
	recode.as <- FALSE
	if (mode == "pass")
		converter <- identity
	else if (mode == "01") {
		if ((!coding || coding == "01") && is.numeric(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			recode.as <- "01"
			converter <- identity
		}
		else if (!is.null(alleles)) {
			converter <- .recode.numeric.by.ref
			recode.as <- "01"
			message("Recoding to 0/1/2 using reference alleles.")
		}
		else {
			converter <- .recode.numeric.by.freq
			recode.as <- "relative"
			message("Recoding to 0/1/2 using empirical frequencies.")
		}
	}
	else if (mode == "relative") {
		if ((!coding || coding == "relative") && is.numeric(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			recode.as <- "relative"
			converter <- identity
		}
		else {
			converter <- .recode.numeric.by.freq
			recode.as <- "relative"
			message("Recoding to 0/1/2 using empirical frequencies.")
		}
	}
	else if (mode == "native")
		if (is.character(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			converter <- identity
		}
	else
		if (!is.null(alleles)) {
			converter <- .recode.character.by.ref
			recode.as <- "native"
		}
	else
		stop("Can only convert genotypes numeric->character given some reference alleles.")
	
	.gty <- unclass(gty)
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
	
	for (a in c("intensity","normalized","filter.sites","filter.samples")) {
		if (!is.null(attr(gty, a)))
			attr(rez, a) <- attr(gty, a)
	}
	class(rez) <- c("genotypes", class(rez))
	attr(rez, "alleles") <- recode.as
	return(rez)
	
}

## replace an existing map with a new one, possibly with allele info
## if <newmap> has rownames, they should map old marker names to new ones
replace.map <- function(gty, newmap, ...) {
	
	if (!inherits(gty, "genotypes") | !.has.valid.map(gty))
		stop("Please supply a genotypes object with map as an attribute.")
	
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
		if (is.list(attr(gty, "intensity")) && length(attr(gty, "intensity")) == 2) {
			x.new <- attr(gty, "intensity")$x[ o, ]
			y.new <- attr(gty, "intensity")$y[ o, ]
			attr(gty, "intensity") <- list(x = x.new, y = y.new)
		}
	}
	
	cols <- colnames(map)
	required <- c("chr","marker","cM","pos")
	attr(gty, "map") <- map[ o,c(required, setdiff(cols, required)) ]
	
	message(paste("...and ending with", nrow(attr(gty, "map")), "markers."))
	if (nrow(attr(gty, "map")) != nrow(gty))
		stop("Map no longer matches genotypes object.")
	
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
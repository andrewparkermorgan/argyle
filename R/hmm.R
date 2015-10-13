## hmm.R
## *very* basic haplotype reconstruction from genotypes

#' Perform (simplistic) haplotype reconstruction, given parental genotypes
#' 
#' @param geno a \code{genotypes} object
#' @param parents a matrix of parental genotypes
#' @param scale scaling parameter for expected recombination fraciton (see Details)
#' @param err genotyping error rate
#' @param ... ignored
#' 
#' @return a dataframe of haplotype segments
#' 
#' @details This function attempts to reconstruct haplotypes of target individuals (in \code{geno})
#' 	in terms of the founder genotypes given in \code{parents}.  It expects founders to be inbred;
#' 	markers at which any founder has a heterozygous genotype are dropped.  (TODO: relax that
#' 	assumption.)
#' 	
#' 	Transition probabilities in the HMM are proportional to expected recombination fractions,
#' 	calculated with the Haldane map function from genetic positions of markers.  The \code{scale}
#' 	parameter is a multiplier for these fractions when the target individuals are >1 generation
#' 	removed from the founders.  Tuning the transition probabilities will be a matter of trial and
#' 	error, and optimal values will probably vary by dataset.
#' 	
#' 	Haplotype segments are pseudo-phased in the result: phase is set to greedily minimize the number of
#' 	haplotype transitions between adjacent segments.  (There is no look-ahead, so in some cases the
#' 	result may be less than optimal.)  Real phasing is beyond the scope of this package at present.
#' 
#'  The model is constrained to only the homozygous states when a sample should really be hemizygous:
#'  this is true for males on chrX and everybody on chrM.
#' 
#' @export
reconstruct.haps <- function(geno, parents, scale = 10, err = 0.001, ...) {
	
	if (!inherits(geno, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (any(sapply(dimnames(parents), is.null)))
		stop("Parental genotypes matrix should have row and column names.")
	
	## keep intersection of markers in target and ref sets
	keep <- intersect(rownames(geno), rownames(parents))
	geno <- geno[ keep,,drop = FALSE ]
	parents <- parents[ keep,,drop = FALSE ]
	
	## drop missing values, markers without genetic position, and parental hets
	drop <- is.na(rowSums(parents)) | is.na(rowSums(geno)) | is.na(attr(geno, "map")$cM)
	drop <- drop | apply(parents, 1, function(x) any(x == 1, na.rm = TRUE))
	message(paste("Dropping", sum(drop), "markers with missing values."))
	geno <- geno[ !drop, ]
	parents <- parents[ !drop, ]
	
	.reconstruct.haps.chr <- function(i) {
		
		g <- geno[ i,, drop = FALSE ]
		p <- parents[ i,,drop = FALSE ]
		
		chr <- attr(g, "map")$chr[1]
		
		message(paste0("Initializing HMM parameters on sequence ", chr, " (", nrow(g), " markers)..."))
		## set model params
		mod <- .init.model(p)
		mod$eprob <- .init.emit.probs(mod, err = err)
		mod$tprob <- .init.trans.probs(mod, markers(g)$cM, scale = scale)
		
		message(paste0("Running Viterbi algorithm on ", ncol(g), " samples..."))
		## initialize container for result
		rez <- matrix(NA, nrow = nrow(g), ncol = ncol(g))
		colnames(rez) <- colnames(g)
		## determine zygosity
		sex <- attr(g, "ped")[ ,"sex" ]
		hemi <- (sex == 1 & grepl("X", chr)) | grepl("M", chr)
		if (any(hemi))
			message("\t(of which ", sum(hemi), " expected to be hemizygous)")
		
		## finally actually run the HMM
		for (j in seq_len(ncol(g))) {
			rez[ ,j ] <- .run.hmm(g[,j], mod = mod, allow.het = !hemi[j])
		}
		
		attr(rez, "model") <- mod
		attr(rez, "map") <- attr(g, "map")
		class(rez) <- c("genohmm.result", class(rez))
		
		return(rez)
		
	}
	
	## split genotypes by chromosome for running HMM
	idx <- split(seq_len(nrow(geno)), droplevels(attr(geno, "map")$chr))
	## loop over chromosomes
	rez <- lapply(idx, .reconstruct.haps.chr)
	## convert haplotypes to blocks
	bl <- do.call(rbind, lapply(rez, as.intervals.genohmm.result))
	return(bl)
	
}

as.intervals.genohmm.result <- function(haps, phase = TRUE, ...) {
	
	rez <- vector("list", ncol(haps))
	for (i in seq_along(rez)) {
		h <- rle(haps[,i])
		bl <- .rle.to.blocks(h)
		tmp <- attr(haps, "model")$labels[ h$values ]
		tmp[ is.na(tmp) ] <- "NA::NA"
		diplo <- t(simplify2array(strsplit(tmp, "::")))
		df <- with(attr(haps, "map"), data.frame(chr = chr[ bl[,1] ],
										   start = pos[ bl[,1] ], end = pos[ bl[,2] ],
										   hap1 = diplo[,1], hap2 = diplo[,2]))
		if (phase)
			df <- pseudophase(df)
		rez[[i]] <- df
	}
	names(rez) <- colnames(haps)
	
	return(plyr::ldply(rez, as.data.frame))
	
}
as.intervals <- function(x, ...) UseMethod("as.intervals")

## convert a run-length encoding to tuples of (start, end)
.rle.to.blocks <- function(x, ...) {
	
	if (!inherits(x, "rle"))
		stop("Please supply an object of class 'rle' (base-R's run-length encoding.)")
	
	len <- x$lengths
	val <- x$values
	
	start <- c(1, (cumsum(len[ -length(len) ]))+1)
	end <- cumsum(len)
	
	cbind(start, end)
	
}

## fill DP matrix and run Viterbi decoding
.run.hmm <- function(obs, mod, ...) {
	
	V <- .fill.dpmat(mod, obs, ...)
	path <- .traceback(V)
	attr(path, "V") <- V
	return(path[ -1 ])
	
}

## given an HMM and an observation sequence, fill out the DP matrix
.fill.dpmat <- function(mod, obs, allow.het = TRUE, ...) {
	
	if (nrow(mod$eprob) != length(obs))
		stop("Emission sequence and emission probabilities must have matching dimensions.")
	obs <- obs+1 ## convert allele codes to 1-based indices
	
	V <- matrix(0, nrow = mod$nstates, ncol = nrow(mod$eprob)+1)
	V[,1] <- 1/nrow(V) # uniform prior
	for (i in 2:ncol(V)) {
		## which was best state last time
		lastmax <- which.max(V[ ,i-1 ])
		if (length(lastmax) == 0) {
			message("Stopping at marker ", i-1)
			print(V[ ,1:i ])
			print(mod$eprob[ i-1,, ])
			print(mod$tprob[ i-2,, ])
			stop("Bad emission or transmission probability.")
		}
		## prob of emission given state
		alpha <- mod$eprob[ i-1,,obs[i-1] ]
		## prob of getting to this state given last best
		if (i > 2)
			beta <- mod$tprob[ i-2, lastmax, ]
		else
			beta <- 1
		## joint prob of this state
		#print(cbind(alpha, beta, lastmax = lastmax))
		gamma <- log(alpha) + log(beta)
		## update
		V[,i] <- gamma + V[ lastmax,i-1 ]
		if (!allow.het)
			V[ mod$het.states,i ] <- -Inf
	}

	return(V)
	
}

## do the Viterbi algorithm on a filled-out DP matrix
.traceback <- function(V, ...) {
	
	path <- apply(V, 2, function(x) which(x == max(x)))
	final <- numeric(length(path))
	final[1] <- path[[1]][1]
	for (i in length(path):2) {
		if (length(path[[i]]) > 1) {
			#print(paste("resolving ambiguity at step", i))
			isect <- intersect(path[[i]], path[[i-1]])
			if (length(isect) > 0)
				final[i] <- isect[1]
			else
				final[i] <- path[[i]][1]
		}
		else
			final[i] <- path[[i]][1]
	}
	
	return(final)
	
}

## initialize an HMM for haplotype reconstruction
.init.model <- function(parents, allow.het = TRUE, ...) {
	
	mod <- list(parents = parents)
	mod$npar <- ncol(parents)
	
	mod$states <- which(upper.tri(diag(mod$npar), diag = TRUE), arr.ind = TRUE)
	mod$het.states <- mod$states[,1] != mod$states[,2]
	labels <- apply(mod$states, 1, function(x) paste(colnames(mod$parents)[x], collapse = "::"))
	rownames(mod$states) <- labels
	mod$labels <- labels
	mod$nstates <- nrow(mod$states)
	mod$nsymb <- 3
	
	class(mod) <- c("genohmm", class(mod))
	return(mod)
	
}

## construct matrix of emission probs at each marker
.init.emit.probs <- function(mod, err = 0.01, ...) {
	
	## set genotyping error rate, per wrong call
	noerr <- 1 - err
	
	## emission probs: sites x states x symbols
	eprob <- array(NA, dim = c(nrow(mod$parents), mod$nstates, 3))
	for (i in seq_len(nrow(mod$parents))) {
		for (j in seq_len(mod$nstates)) {
			expect <- as.integer(mean(mod$parents[ i,mod$states[j,] ]))+1
			if (!is.na(expect))
				eprob[ i,j,expect ] <- noerr
				eprob[ i,j,setdiff(1:3, expect) ]  <- err/2
		}
	}
	
	return(eprob)
	
}

## construct matrix of transition probs between markers
## no principled way to do this in arbitrary, samples; I'll use a scaled version
##	of the Haldane map function to incorporate information on genetic map
.init.trans.probs <- function(mod, cM, scale = 1, ...) {
	
	## get recombination fractions from Haldane formula
	## see <http://www.informatics.jax.org/silver/chapters/7-2.shtml>
	m <- diff(cM/100)
	m <- pmax(m, 1e-5)
	r <- scale * (1/2)*(1-exp(-2*m))
	r <- pmin(0.99999, r)
	
	## transition probs: sites x states x states
	tprob <- array(NA, dim = c(nrow(mod$parents), mod$nstates, mod$nstates))
	#tprob[1,,] <- diag(mod$nstates)
	for (i in seq_along(r)) {
		x <- matrix(0, ncol = mod$nstates, nrow = mod$nstates)
		pairs <- combn(1:mod$nstates, 2)
		for (j in 1:ncol(pairs)) {
			## number of recombinant chrs = number of parents not shared between states
			nr <- 2-length(intersect(mod$states[ pairs[1,j], ], mod$states[ pairs[2,j], ]))
			d <- (r[i]^nr)/ncol(pairs) ## account for above/below diagonal
			x[ pairs[1,j], pairs[2,j] ] <- d
			x[ pairs[2,j], pairs[1,j] ] <- d
		}
		## make sure rows/cols sum to 1
		diag(x) <- 1-colSums(x)
		x <- sweep(x, 1, rowSums(x), "/")
		## update
		tprob[i,,] <- x
	}
	
	return(tprob)
	
}

## pseudo-phasing to limit number of haplotype transitions
## just swaps (arbitrarily-defined) parental haplotype columns
## expect column names: hap1,hap2
pseudophase <- function(df, ...) {
	
	df$hap1 <- as.character(df$hap1)
	df$hap2 <- as.character(df$hap2)
	
	if (nrow(df) <= 1)
		return(df)
	
	for (i in 2:nrow(df)) {
		if ((df$hap1[i] == df$hap2[i-1] || df$hap2[i] == df$hap1[i-1])
			 && df$hap1[i] != df$hap2[i]) {
			x <- df$hap1[i]
			df$hap1[i] <- df$hap2[i]
			df$hap2[i] <- x
		}
	}
	
	return(df)
	
}
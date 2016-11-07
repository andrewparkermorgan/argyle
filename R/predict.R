#' Predict haplogroup (eg chrY or chrM), given labelled controls
#' 
#' @param gty a \code{genotypes} object with just the markers of interest
#' @param train a \code{genotypes} object with labelled training samples
#' 
#' @return a dataframe with predicted haplogroup membership for training samples;
#' 	the result from the LDA procedure (see Details) is attached as \code{attr(,"lda")}.
#' 
#' @details We use linear discriminant analysis via (\code{MASS::lda()}) to extract the most
#' 	informative dimensions of the input data, given the labelled training samples, and then
#' 	project the test samples onto that dimension.
#' 
predict.haplogroup <- function(gty, train, method = c("lda","kmeans"), what = c("intensity","genotypes"), noise = 0.1, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (!inherits(train, "genotypes"))
		stop("Training set should be a 'genotypes' object.")
	
	if (attr(train, "alleles") != "01" || attr(gty, "alleles") != "01")
		stop("Both training and testing sets should have encoding '01'.")
	
	## for filling in missing genotypes
	.impute.missing <- function(f) {
		f[ is.na(f) ] <- mean(f, na.rm = TRUE)
		if (all(is.na(f)))
			f[ is.na(f) ] <- 0
		return(f)
	}
	
	## keep only the markers in both training and testing set
	mk <- intersect(rownames(train), rownames(gty))
	gty <- gty[ mk, ]
	train <- train[ mk, ]
	
	message("Training set is ", ncol(train), " samples x ", nrow(train), " markers.")
	message("Unknown set is ", ncol(gty), " samples x ", nrow(gty), " markers.")
	
	what <- match.arg(what)
	method <- match.arg(method)
	if (method == "lda") {
		
		## step 1: do LDA on training samples
		
		if (what == "intensity") {
			## extract x- and y-intensity into matrix
			genomat <- cbind( t(.copy.matrix.noattr(intensity(train)$x)),
							  t(.copy.matrix.noattr(intensity(train)$y)) )
			## mask missing values to zero (only happens if post-tQN)
			genomat[ is.na(genomat) ] <- 0
		}
		else if (what == "genotypes") {
			## extract genotype matrix
			genomat <- t(.copy.matrix.noattr(train))
			genomat <- apply(genomat, 2, .impute.missing)
			## add some fuzz to allow LDA to work
			genomat <- genomat + rnorm(length(genomat), 0, noise)
		}
	
		## remove low-variance columns, by group
		group.labels <- factor(attr(train, "ped")$fid)
		bygroup <- lapply(split(seq_len(ncol(train)), group.labels),
						  function(i) genomat[i,,drop = FALSE])
		invar <- Reduce("|", lapply(bygroup, colVars, warn.missing = TRUE))
		invar <- (invar < 1e-4)
		genomat <- genomat[ ,!invar ]
		
		## run LDA to train the classifier
		message("Training LDA classifier ... ")
		message("\tretained ", sum(!invar), " feature(s)")
		message("\ton ", length(unique(group.labels)), " group(s)")
		suppressWarnings({ 
			mod <- MASS::lda(genomat, grouping = group.labels,
							 tol = 1e-6)
		})
		
		## step 2: apply classifier to test samples
		message("Applying classifier to unknown samples ...")
		
		if (what == "intensity") {
			testmat <- cbind( t(.copy.matrix.noattr(intensity(gty)$x)),
							  t(.copy.matrix.noattr(intensity(gty)$y)) )
			testmat[ is.na(testmat) ] <- 0
		}
		else if (what == "genotypes") {
			testmat <- t(.copy.matrix.noattr(gty))
			testmat <- apply(testmat, 2, .impute.missing)
			testmat <- testmat + rnorm(length(testmat), 0, noise)
		}
		testmat <- testmat[ ,!invar ]
		
		rez <- MASS:::predict.lda(mod, testmat)
		
		pidx <- cbind(seq_along(rez$class), as.numeric(rez$class))
		rez.df <- transform(as.data.frame(rez$x), iid = rownames(rez$x), cluster = rez$class,
							posterior = rez$posterior[ pidx ])
		rownames(rez.df) <- as.character(rez.df$iid)
		attr(rez.df, "lda") <- rez
		return(rez.df)
		
	}
	else if (method == "kmeans") {
		
		stop("Not implemented yet.")
		
		## step 1: do K-means on training samples
		genomat <- t(.copy.matrix.noattr(train))
		genomat <- apply(genomat, 2, .impute.missing)
		group.labels <- factor(attr(train, "ped")$fid)
		cl <- kmeans(genomat, nlevels(group.labels))
		
		## assign a group label to cluster centroids
		clbygrp <- table(cl$cluster, group.labels)
		print(clbygrp)
		renamer <- apply(clbygrp, 1, function(f) names(which.max(f)))
		print(renamer)
		
		testmat <- t(.copy.matrix.noattr(gty))
		testmat <- apply(testmat, 2, .impute.missing)
		
		
		invar <- rep(FALSE, ncol(genomat))
	}
	

}

recluster <- function(train, test = NULL, ...) {
	
	new1 <- matrix(NA, nrow = nrow(train), ncol = ncol(train),
				   dimnames = dimnames(train))
	if (!is.null(test)) {
		new2 <- matrix(NA, nrow = nrow(test), ncol = ncol(test),
					   dimnames = dimnames(test))
	}
	
	for (ii in seq_len(nrow(train))) {
		cl <- .recluster.marker( attr(train, "intensity")$x[ii,], attr(train, "intensity")$y[ii,], ... )
		new1[ii,] <- cl$cluster
		if (!is.null(test)) {
			xy <- cbind( as.vector(attr(test, "intensity")$x[ii,]),
						 as.vector(attr(test, "intensity")$y[ii,]) )
			pred <- apply(xy, 1, .closest, cl = cl)
			new2[ii,] <- pred
		}
	}
	rez <- list(train = genotypes(new1, map = markers(train), ped = samples(train), alleles = "01",
								  intensity = intensity(train), check = TRUE),
				test = NULL)
	
	if (!is.null(test)) {
		rez$test <- genotypes(new2, map = markers(test), ped = samples(test), alleles = "01",
							  intensity = intensity(test), check = TRUE)
	}
	
	return(rez)
	
}

.recluster.marker <- function(x, y, k = 2) {
	
	X <- cbind(as.vector(x), as.vector(y))
	X[ is.na(X) ] <- 0
	rez <- kmeans(X, k)
	return(rez)
	
}

.closest <- function(x, cl) {
	d <- apply(cl$centers, 1, function(y) sqrt(sum((x-y)^2)))
	return(which.min(d)[1])
}
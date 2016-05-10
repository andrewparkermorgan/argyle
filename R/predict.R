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
predict.haplogroup <- function(gty, train, what = c("intensity","genotypes"), ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (!inherits(train, "genotypes"))
		stop("Training set should be a 'genotypes' object.")
	
	if (attr(train, "alleles") != "01" || attr(gty, "alleles") != "01")
		stop("Both training and testing sets should have encoding '01'.")
	
	## keep only the markers in both training and testing set
	mk <- intersect(rownames(train), rownames(gty))
	gty <- gty[ mk, ]
	train <- train[ mk, ]
	
	message("Training set is ", ncol(train), " samples x ", nrow(train), " markers.")
	message("Unknown set is ", ncol(gty), " samples x ", nrow(gty), " markers.")
	
	what <- match.arg(what)
	if (what == "intensity") {
		
		## step 1: do LDA on training samples
		
		## extract x- and y-intensity into matrix
		genomat <- cbind( t(.copy.matrix.noattr(intensity(train)$x)),
						  t(.copy.matrix.noattr(intensity(train)$y)) )
						  #t(.copy.matrix.noattr(train)) )
		
		## mask missing values to zero (only happens if post-tQN)
		genomat[ is.na(genomat) ] <- 0
		#genomat <- genomat + rnorm(length(genomat), 0, 1e-3)
		
		## remove low-variance columns, by group
		group.labels <- factor(attr(train, "ped")$fid)
		bygroup <- lapply(split(seq_len(ncol(train)), group.labels),
						  function(i) genomat[i,])
		invar <- Reduce("|", lapply(bygroup, function(g) colVars(g) < 1e-4))
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
		
		## now make input object with test samples
		testmat <- cbind( t(.copy.matrix.noattr(intensity(gty)$x)),
						  t(.copy.matrix.noattr(intensity(gty)$y)) )
						  #t(.copy.matrix.noattr(gty)) )
		testmat[ is.na(testmat) ] <- 0

		testmat <- testmat[ ,!invar ]
		rez <- MASS:::predict.lda(mod, testmat)
		
		pidx <- cbind(seq_along(rez$class), as.numeric(rez$class))
		rez.df <- transform(as.data.frame(rez$x), iid = rownames(rez$x), cluster = rez$class,
							posterior = rez$posterior[ pidx ])
		rownames(rez.df) <- as.character(rez.df$iid)
		attr(rez.df, "lda") <- rez
		return(rez.df)
		
	}
	else {
		stop("Not implemented.")
	}
	
}
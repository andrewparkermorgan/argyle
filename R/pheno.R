## pheno.R
## functions for extracting genotype-phenotype associations at a small number of loci

#' Get one-way marker-phenotype association
#' 
#' @param gty a \code{genotypes} object
#' @param marker indexing vector into \code{gty}: can be logical, character or numeric
#' 
#' @return a dataframe with marker metadata, genotype calls, and sample metadata including phenotype
#' 
#' @details To facilitate plotting downstream, the function tries to guess if the phenotype is binary
#' 	(assuming the \code{PLINK} convention of 1=control/2=case), ordinal/categorical, or continuous.
#' 	Categorical phenotypes will be converted to a factor.
#' 
#' @export
oneway <- function(gty, marker, verbose = TRUE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	gty <- gty[ marker, ]
	gty <- recode(gty, "01")
	df <- get.call(gty)
	
	## make friendly-looking genotypes for plotting
	aa <- paste0(df$A1, df$A1)
	ab <- paste0(df$A1, df$A2)
	bb <- paste0(df$A2, df$A2)
	
	aas <- df$call == 0
	abs <- df$call == 1
	bbs <- df$call == 2
	nas <- is.na(df$call)
	
	df$call <- as.character(df$call)
	df$call[ aas & !is.na(aas) ] <- aa[ aas & !is.na(aas) ]
	df$call[ abs & !is.na(abs) ] <- ab[ abs & !is.na(abs) ]
	df$call[ bbs & !is.na(bbs) ] <- bb[ bbs & !is.na(bbs) ]
	df$call[nas] <- NA
	df$call <- factor(df$call, levels = c(unique(aa), unique(ab), unique(bb)))
	
	## guess at phenotype encoding
	.is.integer.like <- function(x) all(x == round(x), na.rm = TRUE)
	if (is.factor(df$pheno) || (.is.integer.like(df$pheno) && max(df$pheno, na.rm = TRUE) <= 2)) {
		## assume case-control
		df$pheno <- factor(df$pheno, levels = c(1,2))
		levels(df$pheno) <- c("control","case")
	}
	else if (.is.integer.like(df$pheno)) {
		## ordinal/categorical
		df$pheno[ df$pheno == -9 ] <- NA
		df$pheno <- factor(df$pheno)
	}
	
	if (verbose)
		print(xtabs(~ call + pheno, df))
	
	return(df[ ,c("chr","marker","cM","pos","call","fid","iid","pheno") ])
	
}

#' Make graphical representation of a 'one-way' genotype-phenotype contingency table
#' 
#' @param df a dataframe from \code{oneway()}
#' @param space numeric scalar; buffer to add between adjacent squares in the plot
#' 
#' @return a \code{ggplot} object equivalent to output from base-\code{R}'s \code{mosaicplot()},
#' 	plotting the (genotype x phenotype) contingency table for a single marker
#' 	
#' @details A two-way contingency table can be rendered as a tesselation of rectangles whose area
#' 	is proportional to the corresponding cell counts.  (This is sometimes called a "fluctiation plot".)
#' 	In this version the plot is divided into columns by genotype, and the width of those columns
#' 	is proportional to genotype frequency.
#' 
#' 	Due to low-level manipulation of axis breaks and labels, this function can only handle
#' 	a single marker at a time.  A warning will be issued if the input contains data from multiple
#' 	markerks, although a (probably bogus) plot will still be drawn.
#' 
#' @seealso \code{\link{oneway}}
#' 
#' @export
oneway.plot <- function(df, space = 1, ...) {
	
	if (!is.factor(df$pheno)) {
		warning("Forcing phenotype to be categorical.")
		df$pheno <- factor(df$pheno)
	}
	
	if (length(unique(df$marker)) > 1) {
		warning("More than a single marker in the input: result will likely be bogus.")
	}
	
	ttl <- paste0(df$marker[1], " @ ", df$chr[1], ": ", round(df$pos/1e6, 2), " Mbp\n")
	
	tbl <- xtabs(~ pheno + call, df)
	bins <- c(0, unname(cumsum(colSums(tbl, na.rm = TRUE))))/sum(tbl, na.rm = TRUE)
	breaks <- bins[ -length(bins) ] + diff(bins)/2
	
	mosaic <- .mosaic.layout(tbl)
	if (is.factor(df$pheno))
		mosaic$pheno <- factor(mosaic$pheno, levels = levels(df$pheno))
	ggplot2::ggplot(mosaic) +
		ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pheno),
				  colour = "white", size = space) +
		ggplot2::scale_y_continuous(label = scales::percent) +
		ggplot2::scale_x_continuous(breaks = breaks, labels = colnames(tbl)) +
		ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(colour = NULL))) +
		ggplot2::ggtitle(ttl) +
		theme_axesonly()
	
}

## do layout for a ggplot version of base-R mosaicplot()
.mosaic.layout <- function(tbl, ...) {
	
	rows <- rowSums(tbl, na.rm = TRUE)
	cols <- colSums(tbl, na.rm = TRUE)
	
	col.box <- cbind( c(0, cumsum(cols)[ -length(cols) ])/sum(cols),
					  cumsum(cols)/sum(cols) )
	row.box <- apply(tbl, 2, function(x) c(0,cumsum(x))/sum(x))
	
	i <- which(!is.na(as.matrix(tbl)), arr.ind = TRUE)
	xmin <- col.box[ i[,2], 1 ]
	xmax <- col.box[ i[,2], 2 ]
	ymin <- row.box[i]
	i[,1] <- i[,1] + 1 
	ymax <- row.box[i]
	
	data.frame(reshape2::melt(tbl), xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
	
}

#' Get two-way (marker x marker)-phenotype association
#' 
#' @param gty a \code{genotypes} object
#' @param markers indexing vector into \code{gty}: can be logical, character or numeric
#' 
#' @return a dataframe with marker metadata, genotype calls, and sample metadata including phenotype
#' 
#' @details A wrapper around \code{oneway()}, generalizing it to marker pairs.  The contingency table
#' 	built from the result of \code{twoway()} will be (genotype x genotype) instead of (genotype x phenotype),
#' 	so the phenotype must be summarizeable.  This function recodes binary phenotypes to 0=control/1=case
#' 	so that taking the mean gives a penetrance value.  Other ordinal or categorical phenotypes will be
#' 	converted from factor to the underlying integer representation, which may or may not be meaningful.
#' 
#' @export
twoway <- function(gty, markers, verbose = TRUE, ...) {
	
	gty <- gty[ markers, ]
	if (nrow(gty) != 2)
		stop("Indexing expression 'markers' must yield exactly two markers.")
	
	df <- oneway(gty, TRUE, verbose = FALSE)
	df.w <- reshape2::dcast(df, iid + pheno ~ marker, value.var = "call")
	## convert binary phenotypes to 0/1 encoding instead of 1/2, so mean=penetrance
	if (is.factor(df.w$pheno) && nlevels(df.w$pheno) == 2)
		df.w$pheno <- as.numeric(df.w$pheno)-1
	tbl <- tapply(df.w$pheno, list(df.w[,3], df.w[,4]), function(x) sum(!is.na(x)))
	
	if (verbose)
		print(tbl)
	
	return(df.w)
	
}

#' Make graphical representation of a 'one-way' genotype-phenotype contingency table
#' 
#' @param df a dataframe from \code{twoway()}
#' @param space numeric scalar; buffer to add between adjacent squares in the plot
#' 
#' @return a \code{ggplot} object equivalent to output from base-\code{R}'s \code{mosaicplot()},
#' 	plotting the (genotype x genotype) contingency table for a single marker, with fill color
#' 	indicating phenotypic mean in each cell
#' 	
#' @details A two-way contingency table can be rendered as a tesselation of rectangles whose area
#' 	is proportional to the corresponding cell counts.  (This is sometimes called a "fluctiation plot".)
#' 	In this version the rectangles are defined by two-locus genotypes, and fill color indicates
#' 	phenotypic mean (or penetrance, in the case of a binary phenotype.)
#' 
#' 	Due to low-level manipulation of axis breaks and labels, this function can only handle
#' 	a single marker at a time.  A warning will be issued if the input contains data from multiple
#' 	markerks, although a (probably bogus) plot will still be drawn.
#' 
#' @seealso \code{\link{twoway}}
#' 
#' @export
twoway.plot <- function(df, space = 1, ...) {
	
	markers <- colnames(df)[3:4]
	tbl <- tapply(df$pheno, list(df[,3], df[,4]), function(x) sum(!is.na(x)))
	x.bins <- c(0, unname(cumsum(colSums(tbl, na.rm = TRUE))))/sum(tbl, na.rm = TRUE)
	x.breaks <- x.bins[ -length(x.bins) ] + diff(x.bins)/2
	y.bins <- c(0, unname(cumsum(rowSums(tbl, na.rm = TRUE))))/sum(tbl, na.rm = TRUE)
	y.breaks <- y.bins[ -length(y.bins) ] + diff(y.bins)/2
	mosaic <- .mosaic.layout(tbl)
	
	df$pheno <- as.numeric(df$pheno)
	tbl2 <- reshape2::melt(tapply(df$pheno, list(df[,3], df[,4]), mean, na.rm = TRUE))
	mosaic <- merge(mosaic[,-3], tbl2)
	
	ggplot2::ggplot(mosaic) +
		ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = value),
						   colour = "white", size = space) +
		ggplot2::scale_y_continuous(breaks = y.breaks, labels = rownames(tbl)) +
		ggplot2::scale_x_continuous(breaks = x.breaks, labels = colnames(tbl)) +
		ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(colour = NULL))) +
		ggplot2::xlab(paste0("\n", markers[2])) +
		ggplot2::ylab(paste0(markers[1], "\n")) +
		theme_axesonly()
	
}
## plots.R
## various plotting functions, mostly for QC

plot.QC.result <- function(qc, show = c("point","label"), max.H = Inf, max.N = Inf, theme.fn = ggplot2::theme_bw, ...) {
	
	calls <- qc$calls
	intens <- qc$intensity
	if (is.null(calls$filter))
		calls$filter <- FALSE
	calls$filter <- with(calls, H > max.H | N > max.N)
	
	show <- match.arg(show)
	
	## plot 1: H vs N
	p1 <- ggplot2::ggplot(calls, ggplot2::aes(x = N, y = H, label = iid, colour = filter))
	if (show == "point")
		p1 <- p1 + ggplot2::geom_point()
	else if (show == "label")
		p1 <- p1 + ggplot2::geom_text()
	p1 <- p1 + ggplot2::scale_colour_manual(values = c("black",scales::muted("red")), na.value = "grey") +
		ggplot2::scale_x_continuous(label = function(x) sprintf("%.1f", x/1e3)) +
		ggplot2::scale_y_continuous(label = function(x) sprintf("%.1f", x/1e3)) +
		ggplot2::guides(colour = FALSE) +
		ggplot2::xlab("\ncount of N calls (x1000)") +
		ggplot2::ylab("count of H calls (x1000)\n") +
		theme.fn()
	
	## plot 2: intensity distribution
	if (is.null(intens))
		return(p1)

	p2 <- ggplot2::ggplot(intens) +
		ggplot2::geom_line(ggplot2::aes(x = iid, y = value, group = q, colour = q)) +
		ggplot2::scale_colour_gradientn("quantile", colours = c(scales::muted("red"), "grey", scales::muted("blue")),
										label = scales::percent) +
		#ggplot2::guides(colour = FALSE) +
		ggplot2::xlab("\nsamples (sorted by median intensity)") +
		ggplot2::ylab("\nintensity quantiles\n") +
		theme.fn() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
									legend.position = c(1,1), legend.justification = c(1,1))
	
	calls.m <- reshape2::melt(calls, id.vars = c("iid","filter"))
	colnames(calls.m) <- c("iid","filter","call","value")
	calls.m$iid <- factor(calls.m$iid, levels = levels(intens$iid))
	p3 <- ggplot2::ggplot(subset(calls.m, call %in% c("H","N"))) +
		ggplot2::geom_point(ggplot2::aes(x = iid, y = value, pch = call, colour = filter), fill = "white") +
		ggplot2::scale_shape_manual(values = c(H=21, N=19)) +
		ggplot2::scale_colour_manual(values = c("black", scales::muted("red")), na.value = "grey") +
		ggplot2::scale_y_continuous(label = function(x) sprintf("%.1f", x/1e3)) +
		ggplot2::guides(colour = FALSE) +
		ggplot2::ylab("number of calls (x1000)\n") +
		theme.fn() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
									panel.grid = ggplot2::element_blank(),
									legend.position = c(1,1), legend.justification = c(1,1))
		
	rez <- gtable:::rbind_gtable( ggplot2::ggplotGrob(p3), ggplot2::ggplotGrob(p2),
								  size = "first")
	panels <- rez$layout$t[grep("panel", rez$layout$name)]
	rez$heights[panels] <- lapply(c(1,2), grid::unit, "null")
	
	return(rez)
	
}

#' Produce a visual summary of QC measures
#' 
#' @param gty a \code{genotypes} object with intensity data attached
#' @param draw actually show the plot, in addition to returning it
#' 
#' @return a \code{gtable} (ineriting from \code{grid::grob}) containg the composite QC plot (see Details)
#' 
#' @details This function will plot any existing QC result attached to \code{gty}; if none is present,
#' 	\code{run.qc.checks(gty)} will be run to generate it.
#' 
#' 	The QC plot has two panels. The upper panel displays the count of heterozygous (H) and missing
#' 	(N, no-call) calls for each sample.  The lower panel displays a the "sum intensity" quantiles of each sample.
#' 	Samples are sorted from left to right in increasing order of median intensity, and the sort order is matched
#' 	between panels.  Samples for which the quality filter is set are shown in dark red in the upper panel.
#' 	
#' 	Although somewhat crude, the count of H and N calls relative to expectations is anecdotally a robust
#' 	measure of genotyping quality (see Didion et al. (2014)).  (But the expectations are important: an outbred
#' 	sample and an inbred sample should have very different numbers of Hs, but probably similar number of Ns.)
#' 	Higher variance in the hybridization intensities, as indicated by wider spread of the intensity quantiles in
#' 	the lower plot, suggests poor input DNA quality.  Contamination of one sample with another, if both had good
#' 	quality DNA, will increase the proportion of no-calls but probably will not shift the intensity quantiles much.
#' 	
#' 	Samples which are diverged from the reference genome used in the array design are expected to have a
#' 	higher proportion of no-calls due to off-target variation in or near the probe sequence.  See Didion
#' 	et al. (2012) for a fuller discussion.
#' 	
#' 	For discussion of the Illumina Infinium chemistry, see Steemers et al. (2006); for more on intensity QC,
#' 	see Staaf et al. (2008).
#' 	
#' @references
#' Steemers FJ et al. (2006) Whole-genome genotyping with the single-base extension assay. Nat Methods 3:31-33.
#' 	doi:10.1038/nmeth842.
#' 
#' Staaf J et al. (2008) BMC Bioinformatics. doi:10.1186/1471-2105-9-409.
#' 
#' Didion JP et al. (2012) Discovery of novel variants in genotyping arrays improves genotype retention and
#' 	reduces ascertainment bias. BMC Genomics 13:34. doi:10.1186/1471-2164-13-34.
#' 	
#' Didion JP et al. (2014) SNP array profiling of mouse cell lines identifies their strains of origin and reveals
#' 	cross-contamination and widespread aneuploidy. BMC Genomics 15:847. doi:10.1186/1471-2164-15-847.
#' 
#' @export
qcplot <- function(gty, draw = TRUE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (is.null(attr(gty, "qc")))
		gty <- run.qc.checks(gty, ...)
	
	p <- plot.QC.result(gty$qc, ...)
	if (draw)
		gtable:::plot.gtable(p)
	
	invisible(p)
	
}

#' Plot 2D hybridization intensities at a few markers
#' 
#' @param a \code{genotypes} object
#' @param markers indexing vector for selecting markers from \code{gty}
#' @param theme.fn a function specifying formatting options for the resulting \code{ggplot}
#' @param force logical; if \code{TRUE}, issues a stern warning before trying to make a huge
#' 	plot which might hang the \code{R} session
#' 	
#' @return 2D intensity plots, one panel per marker (via \code{ggplot})
#' 
#' @details Each point represents a genotype call for a single sample.  Points are coloured according
#' 	to numeric genotype call (0/1/2/N), with missing (N) calls in grey.  Samples for which the
#' 	quality filter is set are shown as open circles.
#' 	
#' 	Since the return value is a \code{ggplot}, the user can add modify it at will.  To recover the
#' 	underlying dataframe, use \code{...$data}.
#'
#' @seealso \code{\link{dotplot.genotypes}} for compact plotting of genotype calls only
#'
#' @export
plot.clusters <- function(gty, markers = NULL, theme.fn = ggplot2::theme_bw, force = FALSE, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes'.")
	
	if (is.null(markers))
		stop("No markers selected.")
	
	gty <- recode.genotypes(gty[ markers, ], "01")
	
	if (nrow(gty) > 10 && !force)
		stop("Resulting plot will have more than 10 panels. Use force = TRUE to proceed anyway.")
	
	df <- get.intensity(gty, markers)
	df <- merge(df, get.call(gty, markers))
	df$call <- as.character(df$call)
	df$call[ is.na(df$call) ] <- "N"
	df$call <- factor(df$call, levels = c(0,1,2,"N"))
	if (!nrow(df))
		stop("No data.")
	
	df$filter <- attr(gty, "filter.samples")[ as.character(df$iid) ]
	
	ggplot2::ggplot(df) +
		ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = call, shape = filter), fill = "white") +
		ggplot2::scale_shape_manual(values = c(19, 21)) +
		ggplot2::scale_colour_manual(values = c(RColorBrewer::brewer.pal(3, "Set1"), "grey")) +
		ggplot2::guides(shape = FALSE) +
		ggplot2::facet_wrap(~ marker) +
		ggplot2::coord_equal() +
		theme.fn()
	
}

#' Show a visual representation of a slice of a genotype matrix
#' 
#' @param a \code{genotypes} object
#' @param size numeric; size of points
#' @param shape if \code{"point"}, each genotype call is represented by a circle; if \code{"tile"},
#' 	the effect will be like a heatmap
#' @param meta dataframe of additional sample metadata to be joined to the input before plotting
#' @param force logical; if \code{TRUE}, issues a stern warning before trying to make a huge
#' 	plot which might hang the \code{R} session
#' 	
#' @return a grid of genotype calls, laid out horizontally in genomic coordinates
#' 
#' @details Each point represents a genotype call for a single sample.  Points are coloured according
#' 	to numeric genotype call (0/1/2/N), with missing (N) calls hidden.  Points are spaced evenly along
#' 	the x-axis to facilitate visual inspection; marker spacing in genomic coordinates is indicated in
#' 	a track along the bottom of the plot which connects each column to its relative genomic position.
#' 	
#' 	Since the return value is a \code{ggplot}, the user can add modify it at will.  To recover the
#' 	underlying dataframe, use \code{...$data}.
#'
#' @seealso \code{\link{plot.clusters}} for intensity plots by marker
#'
#' @export
dotplot.genotypes <- function(gty, size = 2, meta = NULL, shape = c("point","tile"), force = FALSE, ...) {
	
	## check that input is right sort
	if (!(inherits(gty, "genotypes") && .has.valid.map(gty)))
		stop("Please supply an object of class 'genotypes' with valid marker map.")
	
	if (prod(dim(gty)) > 10e3 && !force)
		stop("Resulting plot will have >10k points. Use force = TRUE to proceed anyway.")
	
	## set up horizontal positions of markers
	map <- attr(gty, "map")
	map$i <- as.numeric(factor(map$pos))
	map$ipos <- map$pos - min(map$pos)
	map$relpos <- map$ipos/diff(range(map$pos)) * max(map$i)
	map <- map[ with(map, order(i)), ] # important
	
	## convert genotypes to long-form for ggplot
	df <- reshape2::melt(gty)
	colnames(df) <- c("marker","id","ibs")
	if (is.numeric(gty))
		df$ibs.code <- factor(df$ibs, levels = c(0,1,2))
	else
		df$ibs.code <- factor(df$ibs, levels = c("A","C","G","T","H"))
	df <- merge(df, map)
	
	## add family info, if present
	if (!is.null(attr(gty, "ped")))
		df <- merge(df, attr(gty, "ped"), by.x = "id", by.y = "iid", all.x = TRUE)
	
	## add metadata, if provided
	if (!is.null(meta))
		df <- merge(df, meta, all.x = TRUE)
	
	## check that dimensions still match
	stopifnot(nrow(df) == prod(dim(gty)))
	
	shape <- match.arg(shape)
	if (shape == "tile")
		geom.fn <- ggplot2::geom_tile
	else
		geom.fn <- ggplot2::geom_point
	
	p <- ggplot2::ggplot(subset(df, !is.na(ibs))) +
		geom.fn(ggplot2::aes(x = i, y = id, fill = ibs.code), pch = 21, size = size) +
		ggplot2::geom_segment(data = map,
							  ggplot2::aes(x = i, xend = relpos, y = 0.8, yend = 0), colour = "grey70") +
		ggplot2::geom_hline(yintercept = 0) +
		ggplot2::scale_x_continuous(breaks = warped.breaks(map$relpos, map$relpos, map$pos),
						   labels = round(warped.breaks(map$relpos, map$relpos, map$pos, scale = "transformed")/1e6, 2)) +
		ggplot2::theme_minimal() + ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
		ggplot2::xlab("\nposition (Mbp)")
	
	if (is.numeric(gty))
		p <- p + ggplot2::scale_fill_manual("call", values = c("white","grey","black"), drop = FALSE)
	else
		p <- p + ggplot2::scale_fill_manual("call", values = RColorBrewer::brewer.pal(9, "Set1")[ c(1:4,9) ], drop = FALSE)
	
	attr(p, "map") <- map
	
	return(p)
	
}
dotplot <- function(x, ...) UseMethod("dotplot")
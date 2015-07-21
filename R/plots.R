## plots.R
## various plotting functions, mostly for QC

plot.QC.result <- function(qc, show = c("point","label"), theme.fn = ggplot2::theme_bw, ...) {
	
	calls <- qc$calls
	intens <- qc$intensity
	if (is.null(calls$filter))
		calls$filter <- FALSE
	
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
#' 	\code{run.sample.qc(gty)} will be run to generate it.
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
		gty <- run.sample.qc(gty, ...)
	
	p <- plot.QC.result(gty$qc, ...)
	if (draw) {
		if (inherits(p, "ggplot"))
			plot(p)
		else if (inherits(p, "gtable"))
			gtable:::plot.gtable(p)
	}
	
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
	
	df$filter <- is.filtered(gty)$samples[ as.character(df$iid) ]
	
	ggplot2::ggplot(df) +
		ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = call, shape = filter), fill = "white") +
		ggplot2::scale_shape_manual(values = c(19, 21)) +
		ggplot2::scale_colour_manual(values = c(RColorBrewer::brewer.pal(3, "Set1"), "grey")) +
		ggplot2::guides(shape = FALSE) +
		ggplot2::facet_wrap(~ marker) +
		ggplot2::coord_equal() +
		theme.fn()
	
}

#' @export
heatmap <- function(gty, ...) {
	
	if (!(inherits(gty, "genotypes")))
		stop("Please supply an object of class 'genotypes'.")
	
	## recode genotypes to numeric
	gty <- recode(gty, "01")
	
	## get distance matrix
	message("Computing distance matrix from genotypes...")
	d <- dist(t(gty), "manhattan")
	cl <- hclust(d)
	k <- as.matrix(d)/(2*nrow(gty))
	k <- k[ cl$order,cl$order ]
	
	ggplot2::ggplot(reshape2::melt(1-k)) +
		ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
		scale_fill_heatmap("P(IBS)") +
		ggplot2::coord_equal() +
		theme_heatmap()
		
	
}

#' Plot histogram of (sum-)intensities by sample
#' @export
intensityhist <- function(gty, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes' with valid intensity matrices attached.")
	
	df <- get.intensity(gty, TRUE)
	df$si <- with(df, sqrt(x^2 + y^2))
	ggplot2::ggplot(df) +
		ggplot2::geom_histogram(ggplot2::aes(x = si), binwidth = diff(range(df$si, na.rm = TRUE))/50) +
		ggplot2::facet_wrap(~ iid, ncol = 3) +
		ggplot2::ylab("count of markers\n") + 
		ggplot2::xlab(expression(atop("", paste("sum-intensity = ", sqrt(x^2 + y^2)))))
	
}

#' Plot B-allele frequency (BAF) and log2-intensity ratio (LRR) for a sample
#'
#' @param gty a \code{genotypes} object with BAF and LRR pre-computed via \code{tQN()}
#' @param sm indexing vector which should extract exactly one sample from \code{gty}
#' @param ... additional parameters passed on to \code{supsmu()} to control smoothing
#'
#' @return a two-panel plot, BAF in upper panel and LRR in lower panel
#' 
#' @details For a detailed description of the BAF and LRR metrics, their calculation,
#' 	and literature references, see \code{?tQN}.
#' 
#' @seealso \code{\link{tQN}}
#'
#' @export
bafplot <- function(gty, sm = TRUE, smooth = "cv", draw = TRUE, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.baflrr(gty)))
		stop("Please supply an object of class 'genotypes' with valid BAF and LRR pre-computed.")
	
	gty <- gty[,sm]
	if (ncol(gty) != 1)
		stop("BAF+LRR plot only implemented for one sample at a time.")
	
	if (!is.numeric(gty))
		gty <- recode.genotypes(gty, "01")
	
	baf <- get.baf(gty, TRUE)
	
	## apply smoothing to BAF and LRR
	baf <- plyr::ddply(baf, .(chr), function(d) {
		d$BAF.smooth <- .smooth.me(d$pos, d$BAF, ...)
		d$LRR.smooth <- .smooth.me(d$pos, d$LRR, ...)
		return(d)
	})
	
	## plot BAF, coloured by call
	baf$.call <- factor( ifelse(is.na(baf$call), NA, ifelse(baf$call == 1, 1, ifelse(as.integer(baf$BAF > 0.5), 2, 0))),
						 levels = c(0,1,2) )
	baf$.col <- factor(baf$chr):factor(baf$.call)
	call.cols <- rep_len( c("lightblue","mediumpurple1","pink","darkblue","purple4","red"),
						  length.out = nlevels(baf$.col) )
	
	p1 <- ggmanhattan(baf) +
		ggplot2::geom_point(ggplot2::aes(y = BAF, colour = .col)) +
		ggplot2::geom_point(ggplot2::aes(y = BAF.smooth), colour = "red") +
		ggplot2::scale_alpha_discrete(range = c(0.5,1)) +
		ggplot2::scale_colour_manual(values = call.cols, na.value = "grey") +
		ggplot2::guides(colour = FALSE, alpha = FALSE) +
		ggplot2::theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
											 axis.text.x = ggplot2::element_blank(),
											 plot.margin = grid::unit(c(1, 1, 0, 0.5), "lines")) +
		ggplot2::ggtitle(paste0(colnames(gty)[1],"\n"))
	
	## plot LRR with smoothed fit
	lrr.cols <- brewer.interpolate("Spectral")(6)
	p2 <- ggmanhattan(baf) +
		ggplot2::geom_point(ggplot2::aes(y = LRR), colour = "grey") +
		ggplot2::geom_point(ggplot2::aes(y = LRR.smooth), colour = "red") +
		#ggplot2::scale_colour_gradientn(colours = lrr.cols, breaks = c(-1, 0, 1)) +
		ggplot2::scale_y_continuous(limits = c(-2,2)) +
		ggplot2::guides(colour = FALSE) +
		#scale_x_genome() +
		ggplot2::theme_bw() +ggplot2::theme(axis.title.x = ggplot2::element_blank(),
											plot.margin = grid::unit(c(0.2, 1, 0.5, 0.5), "lines"))
	
	## combine plots
	rez <- gtable:::rbind_gtable(ggplot2::ggplotGrob(p1),
								 ggplot2::ggplotGrob(p2),
								 size = "first")
	
	if (draw) {
		grid::grid.newpage()
		grid::grid.draw(rez)
	}
	
	invisible(rez)
	
}

## smoothing helper
.smooth.me <- function(x, y, smooth = "cv", ...) {
	smth <- do.call(cbind, supsmu(x, y, span = smooth, ...))
	z <- rep(NA, length(y))
	i <- match(x, smth[ ,1 ], nomatch = 0)
	z[ i > 0] <- smth[ i,2 ]
	return(z)
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
#' 	At present, the function will only do a single chromosome at a time.
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
	if (length(unique(map$chr)) > 1)
		stop("Dotplot works only for one chromosome at a time.")
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


#' Create skeleton of a 'Manhattan plot' (concatenated chromosomes) with \code{ggplot2}
#' 
#' @param df a dataframe with columns (at least) "chr" and "pos"
#' @param chroms chromosome names (in order); will try to guess them if not supplied
#' @param scale one of \code{"Mbp"} (millions of basepairs; default) or \code{"cM"} (centimorgans)
#' @param space padding factor to add in between chromosomes
#' @param cols vector of length 2, giving alternating colours for alternating chromosomes
#' 
#' @return a \code{ggplot2} object whose x-axis is defined by chromosome and position (with
#' 	chromosomes concatenated in karyotype order), to which data in \code{df} can be added as
#' 	additional layers
#' 
#' @export
ggmanhattan <- function(df, chroms = NULL, scale = c("Mbp", "cM"), space = NULL, cols = c("grey30","grey60"), ...) {
	
	if (all(c("chr","pos") %in% colnames(df)))
		df <- df[ ,c("chr","pos", setdiff(colnames(df), c("chr","pos"))) ]
	else
		stop("We require a dataframe with at least columns 'chr' and 'pos'.")
	
	scale <- match.arg(scale)
	if (scale == "Mbp") {
		denom <- 1e6
		if (is.null(space))
			space <- 10
	}
	else {
		denom <- 1
		if (is.null(space))
			space <- 5
	}
		
	if (!is.factor(df[,1]))
		if (!is.null(chroms))
			df[,1] <- factor(df[,1], levels = chroms)
	else
		df[,1] <- factor(df[,1])
	else
		df[,1] <- factor(df[,1])
	
	chrlen <- tapply(df[,2], df[,1], max, na.rm = TRUE)/denom + space
	adj <- c(0, cumsum(chrlen))
	ticks <- adj[ -length(adj) ] + diff(adj)/2
	#names(adj) <- c(names(chrlen),"z")
	#print( cbind(chrlen, adj[-1], ticks) )
	#print(length(chrlen))
	#print(length(ticks))
	df$.adj <- adj[ as.numeric(df[,1]) ]
	df$.x <- df$.adj + df[,2]/denom
	df$.chr <- df[,1]
	colmap <- setNames( rep_len(c(0,1), nlevels(df$.chr)), levels(df$.chr)  )
	df$.colour <- factor(colmap[ df$.chr ])
	#print(tapply(df$.x, df$.chr, max))
	
	rez <- ggplot2::ggplot(df, ggplot2::aes(x = .x, colour = .colour))
	rez <- rez +
		ggplot2::scale_x_continuous(breaks = ticks, minor_breaks = adj,
						   labels = gsub("^chr", "", names(chrlen))) +
		ggplot2::scale_colour_manual(values = cols) +
		ggplot2::guides(colour = FALSE) +
		ggplot2::theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
						   panel.grid.minor = ggplot2::element_line(colour = "grey90"),
						   panel.grid.minor = ggplot2::element_blank())
	return(rez)
	
}

#' Auto-plotting of a PCA result
#' 
#' @param pc a \code{pca.result} object
#' @param K numeric vector of length 2 specifing which PCs to plot against each other; first is on x-axis and second on y-axis
#' @param screeplot logical; if \code{TRUE}, show both the plot of two PCs against each other and the variances explained of
#' 	all available PCs
#' @param show length-1 character vector; show just points, or sample IDs.  Use \code{"nothing"} to have all the aesthetics
#' 	and axes set, but not actually draw anything.
#' @param theme.fn a \code{ggplot2}-compatible function to specify formatting
#' 
#' @return if \code{screeplot = FALSE}, a \code{ggplot}; if \code{screeplot = TRUE}, a \code{gtable::gtable} object which can be
#' 	modified or re-rendered with \code{grid::grid.draw()}.
#' 
#' @export
plot.pca.result <- function(pc, K = c(1,2), screeplot = FALSE, show = c("points","labels","nothing"), theme.fn = ggplot2::theme_bw, ...) {
	
	if (!inherits(pc, "pca.result"))
		stop("Please supply an object of class 'pca.result'.")
	
	if (length(K) != 2)
		stop("Can only plot pairs of PCs.")
	
	show <- match.arg(show)
	if (show == "points") {
		geom.fn <- ggplot2::geom_point
		size <- 2
	}
	else if (show == "labels") {
		geom.fn <- ggplot2::geom_text
		size <- 3.5
	}
	else if (show == "nothing") {
		geom.fn <- ggplot2::geom_blank
		size <- NULL
	}
		
	xll <- paste0("\nPC", K[1], " (", sprintf("%.1f", 100*attr(pc, "explained")[ K[1] ]), "%)")
	yll <- paste0("PC", K[2], " (", sprintf("%.1f", 100*attr(pc, "explained")[ K[2] ]), "%)\n")
	
	## plot PCs against each other
	p1 <- ggplot2::ggplot(pc, ggplot2::aes_string(paste0("PC", K[1]), paste0("PC", K[2]),
												  label = "iid"), size = size) +
		geom.fn() +
		ggplot2::xlab(xll) + ggplot2::ylab(yll) +
		ggplot2::coord_equal() +
		theme.fn()
	
	if (screeplot) {
		
		## draw "screeplot" (PCs by variance explained)
		evdf <- data.frame(PC = as.integer(gsub("PC","",grep("^PC", colnames(pc), value = TRUE))),
						   variance = attr(pc, "explained"))
		
		p2 <- ggplot2::ggplot(evdf, ggplot2::aes(x = PC, y = variance)) +
			ggplot2::geom_line(colour = "darkblue") +
			ggplot2::geom_point(colour = "darkblue") +
			ggplot2::scale_y_continuous(label = scales::percent) +
			ggplot2::scale_x_continuous(breaks = 1:max(evdf$PC, na.rm = TRUE),
										labels = function(x) paste0("PC",round(x)),
										expand = c(0,0.5)) +
			ggplot2::ylab("% variance explained\n") +
			theme.fn() + ggplot2::theme(axis.title.x = ggplot2::element_blank())
		
		p1 <- ggplot2::ggplotGrob(p1)
		p2 <- ggplot2::ggplotGrob(p2)
		rez <- gtable:::cbind_gtable( p1, p2,
									  size = "first")
		panels <- rez$layout$t[grep("panel", rez$layout$name)]
		#rez$widths[panels] <- lapply(c(1,1), grid::unit, "null")
		
		grid::grid.newpage()
		grid::grid.draw(rez)
		invisible(rez)
		
	}
	else {
		return(p1)
	}
	
}

#' Plot frequencies of missing, heterozygous, and minor-allele calls across genome
#' 
#' @param gty a \code{genotypes} object
#' @param max.H show markers with heterozygosity greater than this threshold
#' @param max.N show markers with missingness greater than this threshold
#' @param max.H show markers with minor-allele frequency less than this threshold
#' 
#' @return a \code{ggplot} object with a Manhattan-style plot of the above values
#' 
#' @seealso \code{\link{qcplot}}, \code{\link{bafplot}}, \code{\link{summarize.calls}}
#' 
#' @export
freqplot <- function(gty, max.H = -1, max.N = -1, min.maf = Inf, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	gty <- recode.genotypes(gty, "relative")
	calls <- summarize.calls(gty, "markers", counts = FALSE)

	calls <- merge(calls, attr(gty, "map"))
	calls <- transform(calls, maf = B+0.5*H)
	calls <- subset(calls[ ,c("chr","pos","H","N","maf") ], H > max.H | N > max.N | maf <= min.maf)
	calls.m <- reshape2::melt(calls, id.vars = c("chr","pos"))
	
	call.labs <- c("minor allele","H","N")
	call.cols <- c("black","grey",scales::muted("red"))
	calls.m$variable <- factor(calls.m$variable, levels = c("maf","H","N"),
							   labels = call.labs)

	message("Markers failing by")
	message("\t", sprintf("%13s", "no-call rate:"), sprintf("%7d", sum(calls$N > max.N)))
	message("\t", sprintf("%13s", "het rate:"), sprintf("%7d", sum(calls$H > max.H)))
	message("\t", sprintf("%13s", "MAF:"), sprintf("%7d", sum(calls$maf <= min.maf)))
	message("\t", sprintf("%13s", "Total"), sprintf("%7d", nrow(calls)))
	
	ggmanhattan(calls.m) +
		ggplot2::geom_point(ggplot2::aes(y = value, colour = variable)) +
		ggplot2::scale_colour_manual("frequency of", values = call.cols, labels = call.labs) +
		ggplot2::facet_grid(variable ~ .) +
		ggplot2::ylab("relative frequency\n")
	
}

#' Plot the results of haplotype reconstruction
#' 
#' @param haps a dataframe containing the results of \code{reconstruct.haps()}
#' @param space numeric scalar; how much buffer to add between adjacent chromosomes
#' 
#' @return A \code{ggplot} object with haplotype mosaics arranged in a grid layout
#' 
#' @details The input dataframe is assumed to have the following columns: \code{.id}, \code{chr},
#' \code{start}, \code{end}, \code{hap1} and \code{hap2}.  Since the result is a \code{ggplot},
#' it can be futher decorated/modified at will.
#' 
#' @export
plot.haplotypes <- function(haps, space = 0.1, ...) {
	
	haps <- droplevels(haps)
	if (!is.factor(haps$chr))
		haps$chr <- factor(haps$chr)
	haps$chr <- factor(haps$chr, levels = rev(levels(haps$chr)))
	haps$space <- space
	
	ggplot2::ggplot(haps, ggplot2::aes(xmin = start, xmax = end)) +
		# 'maternal' haplotypes
		ggplot2::geom_rect(ggplot2::aes(ymin = as.numeric(chr)-(0.5-space), ymax = as.numeric(chr)-space/2, fill = hap1)) +
		# 'paternal' haplotypes
		ggplot2::geom_rect(ggplot2::aes(ymin = as.numeric(chr)+space/2, ymax = as.numeric(chr)+(0.5-space), fill = hap2)) +
		scale_x_genome() +
		ggplot2::scale_y_continuous(breaks = seq_len(nlevels(haps$chr)), labels = gsub("^chr", "", levels(haps$chr))) +
		ggplot2::scale_fill_discrete("founder") +
		ggplot2::facet_wrap(~ .id, ncol = 4) +
		ggplot2::theme_bw()
	
}
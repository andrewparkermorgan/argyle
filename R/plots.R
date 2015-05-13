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
	grid::grid.draw(rez)
	
}

qcplot <- function(gty, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	if (is.null(attr(gty, "qc")))
		gty <- run.qc.checks(gty, ...)
	
	plot.QC.result(gty$qc, ...)
	
}

plot.clusters <- function(gty, markers = NULL, theme.fn = ggplot2::theme_bw, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("Please supply an object of class 'genotypes'.")
	
	if (is.null(markers))
		stop("No markers selected.")
	
	gty <- recode.genotypes(gty[ markers, ], "01")
	
	if (nrow(gty) > 10)
		warning("Grabbing a lot of markers -- this might be a bad idea.")
	
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
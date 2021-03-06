## themes.R
## auxiliary plotting functions, graphical themes, etc.

## map breaks on 'factor' scale (evenly-spaced) to breaks on continus scale (unevenly-spaced)
warped.breaks <- function(x, y1, y2, scale = c("native","transformed"), ...) {
	
	stopifnot(length(y1) == length(y2))
	i <- findInterval(x, y1)
	b <- pretty(y2[i])
	b <- b[ b > min(y2) & b < max(y2) ]
	
	scale <- match.arg(scale)
	if (scale == "native")
		return(y1[ findInterval(b, y2) ])
	else
		return(b)
}

## scales & transformations
trans_genome <- function(scale = 1e6, ...) {
	scales::trans_new("genome coordinates", transform = function(x) x/scale, 
			  inverse = function(X) X * scale)
}

#' Custom x-axis scale for genome coordinates.
#' 
#' @param scale scale factor, in base pairs
#' @param unit abbreviation for units corresponding to scale factor (see Details)
#' @param ... additional arguments to \code{ggplot2::scale_x_continuous()}
#' 
#' @return a \code{ggplot2} scale function
#' 
#' @details For example, if \code{scale = 1000}, we might set \code{unit = "kbp"}.
#' 
#' @seealso \code{\link[ggplot2]{scale_x_continuous}}
#' 
#' @export
scale_x_genome <- function(scale = 1e6, unit = NULL, ...) {
	
	if (is.null(unit)) {
		if (scale == 1e3)
			unit <- "kbp"
		else if (scale == 1e6)
			unit <- "Mbp"
		else if (scale == 1e9)
			unit <- "Gbp"
	}
	
	gtrans <- trans_genome(scale, unit)
	return(list( ggplot2::scale_x_continuous(trans = gtrans, labels = function(x) gtrans$transform(x), ...),
				 ggplot2::xlab(paste0("\nposition (", unit, ")")) ))
	
}

neglog_trans <- function (base = exp(1)) 
{
	trans <- function(x) -1*log(x, base)
	inv <- function(x) base^(-1*x)
	scales::trans_new(paste0("-log-", format(base)), trans, inv, scales::log_breaks(base = base), 
			  domain = c(1e-100, Inf))
}

#' Custom y-axis scale for -log10(p-values)
#' 
#' @param ... additional arguments to \code{ggplot2::scale_y_continuous()}
#' 
#' @return a \code{ggplot2} scale function
#' 
#' @seealso \code{\link[ggplot2]{scale_x_continuous}}
#' 
#' @export
scale_y_logp <- function(...) {
	
	ggplot2::scale_y_continuous(..., trans = neglog_trans(10))
	
}

## themes

#' A very minimal graphics theme
#' 
#' @param base_size baseline font size in points
#' @param base_family baseline font family
#' 
#' @return a \code{ggplot2} theme function
#' 
#' @seealso \code{\link[ggplot2]{theme}}
#' 
#' @export
theme_nothing <- function(base_size = 12, base_family = "Helvetica")
{
	ggplot2::`%+replace%`(
		ggplot2::theme_bw(base_size = base_size, base_family = base_family),
		ggplot2::theme(
			rect             = ggplot2::element_blank(),
			line             = ggplot2::element_blank(),
			text             = ggplot2::element_blank(),
			axis.ticks.margin = grid::unit(0, "lines")
		))
}

#' A genome-browser like graphics theme
#' 
#' @param base_size baseline font size in points
#' @param base_family baseline font family
#' 
#' @return a \code{ggplot2} theme function
#' 
#' @export
theme_gbrowse <- function(base_size = 12, base_family = "Helvetica") {
	ggplot2::`%+replace%`(
		ggplot2::theme_bw(base_size = base_size, base_family = base_family),
		ggplot2::theme(
			panel.grid = ggplot2::element_blank(),
			panel.border = ggplot2::element_blank(),
			axis.line.y = ggplot2::element_blank(),
			axis.ticks.y = ggplot2::element_blank(),
			axis.text.y = ggplot2::element_blank(),
			axis.title.y = ggplot2::element_blank(),
			axis.ticks.margin = grid::unit(0.25, "lines")
		))
}

#' A graphics theme with only an x-axis
#' 
#' @param base_size baseline font size in points
#' @param base_family baseline font family
#' 
#' @return a \code{ggplot2} theme function
#' 
#' @seealso \code{\link[ggplot2]{theme}}
#' 
#' @export
theme_axesonly <- function(base_size = 12, base_family = "Helvetica") {
	ggplot2::`%+replace%`(
		ggplot2::theme_bw(base_size = base_size, base_family = base_family),
		ggplot2::theme(
			panel.grid = ggplot2::element_blank(),
			panel.border = ggplot2::element_blank(),
			axis.line.y = ggplot2::element_blank(),
			#axis.ticks.y = ggplot2::element_blank(),
			#axis.text.y = ggplot2::element_blank(),
			#axis.title.y = ggplot2::element_blank(),
			axis.ticks.margin = grid::unit(0.25, "lines")
		))
}

#' A graphics theme suitable for heatmaps
#' 
#' @param ... ignored
#' 
#' @return a \code{ggplot2} theme function
#' 
#' @seealso \code{\link[ggplot2]{theme}}
#' 
#' @export
theme_heatmap <- function(...) {
	
	ggplot2::`%+replace%`(
		ggplot2::theme_classic(),
		ggplot2::theme(axis.text = ggplot2::element_text(size = 8),
					   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
					   axis.line = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
	)
		
}

#' A red-to-orange color scale for heatmaps
#' 
#' @param ... additional arguments passed to \code{ggplot2::scale_fill_gradientn()}
#' 
#' @return a \code{ggplot2} scale function
#' 
#' @export
scale_fill_heatmap <- function(...) {
		
	heat.cols <- colorRampPalette(RColorBrewer::brewer.pal(4, "Spectral")[2:1])(4)
	ggplot2::scale_fill_gradientn(..., colours = heat.cols, na.value = "grey70")
	
}

#' Make a blank grob for plotting
#' 
#' @param ... ignored
#' 
#' @return a \code{grid::grob} object
#' 
#' @seealso \code{\link[grid]{grid.rect}}
#' 
#' @export
blank.grob <- function(...) {
	grid::grid.rect(gp = grid::gpar(col = NA))
}

## get an interpolator function for a RColorBrewer palette
brewer.interpolate <- function(palette, ...) {
	
	maxcol <- RColorBrewer::brewer.pal.info[ palette,"maxcolors" ]
	pal <- colorRampPalette(RColorBrewer::brewer.pal(maxcol, palette))
	return(pal)
	
}
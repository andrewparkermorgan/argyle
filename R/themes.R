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

#' @export
scale_y_logp <- function(...) {
	
	ggplot2::scale_y_continuous(..., trans = neglog_trans(10))
	
}

## themes
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
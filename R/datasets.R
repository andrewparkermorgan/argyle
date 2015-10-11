#' Genotype calls and intensities from the GigaMUGA array for mouse
#'
#' A \code{genotypes} object containing genotype calls and hybridization intensities for
#' 116 control samples (67 male, 49 female) used to calibrate the GigaMUGA array.  To
#' keep the object at reasonable size, only calls from chr10 are included (6811 markers).
#'
#' The samples span 8 inbred mouse strains (n = 58) denoted with 2-letter codes:
#' \itemize{
#' \item AA -- A/J
#' \item BB -- C57BL/6J
#' \item CC -- 129S1/SvImJ
#' \item DD -- NOD/ShiLtJ
#' \item EE -- NZO/HlLtJ
#' \item FF -- CAST/EiJ
#' \item GG -- PWK/PhJ
#' \item HH -- WSB/EiJ
#' }
#' The remaining 58 samples are F1s between these 8 strains, also denoted with 2-letter codes
#' corresponding to their two parental strains.
#' 
#' The dataset \code{inbreds} contains just the 58 inbred samples.
#' 
#' @aliases inbreds
#' 
#' @usage
#' 	data(ex); data(inbreds)
#' @source Morgan AP et al. (2015) The Mouse Universal Genotying Array: from substrains
#' 	to subspecies. Submitted to G3.
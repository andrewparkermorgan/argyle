#' argyle: An \code{R} package for import and QC of genotypes from Illumina Infinium arrays
#' 
#' @section The \code{genotypes} object:
#' The \code{genotypes} class is just a matrix (sites x samples) with row and column names, and a dataframe
#' (in \code{attr(,"map")}) describing marker positions.  Nearly all the functions in this package expect
#' a \code{genotypes} object as input.  For details, see \code{\link{genotypes}}.
#'
#' @section Accessing metadata:
#' The \code{`$`} operator is overloaded for the \code{genotypes} class, so that writing \code{attr(g, "normalized")}
#' is equivalent to writing \code{g$normalized}.  But this works only in one direction: \code{g$normalized <- FALSE}
#' fails.  Nefarious users can modify attributes directly using the standard and somewhat convoluted syntax
#' \code{attr(,"x") <- y} but do so at their own risk.  For safety, always check that the resulting object remains
#' valid (all internal parts having matching dimensions and names) with a call to \code{validate(g)}.
#'
#' Accessor functions are provided for the marker map (\code{markers(g)}), sample metadata (\code{samples(g)}),
#' and intensity matrices (\code{intensity(g)}).
#' 
#' @section Notes on allele encoding:
#' For the purposes of this package, all markers on an array are treated as biallelic SNPs, and all samples are assumed
#' to be diploid for the autosomes. Genotype calls are reported by Illumina BeadStudio as a two-character vector of
#' nucleotides: eg. \code{AA}, \code{AG}, \code{GG} for an [A/G] SNP.  The \code{-} character indicates a missing call
#' ("no-call").  On import, these calls are summarized to a single character, one of \code{ACGTHN}
#' (\code{H} = heterozygous, \code{N} = no-call).
#' 
#' For most analyses a numeric representation of genotypes is desirable.  The function \code{recode.genotypes()}
#' performs this conversion.  When reference alleles are provided in columns "A1" (REF) and "A2" (ALT) in the marker map,
#' genotypes are recoded 0 (homozygous REF), 1 (heterozygous), 2 (homozygous ALT) or \code{NA} (missing).
#' 
#' Recoding can also be "relative": that is, performed with respect to the major and minor allele as defined by
#' the dataset itself.  In this case the recoded 0 (homozygous major allele), 1 (heterozygous), 2 (homozygous
#' minor allele) or \code{NA} (missing).
#'
#' The attribute \code{alleles} tracks the current allele encoding.
#'
#' @docType package
#' @name argyle
#' @importFrom Rcpp sourceCpp
#' @useDynLib argyle, .registration=TRUE
NULL

#' Constructor for a \code{genotypes} object
#'
#' @param G a genotype matrix with markers in rows and samples in columns, with both row and column names
#' @param map a valid marker map (see Deatils) corresponding to \code{G}, with row names
#' @param ped a valid "pedigree" (dataframe containing sample metadata)
#' @param alleles character vector describing allele encoding (see \code{\link{argyle}} for details); \code{"auto"}
#' 	lets the package try to guess the encoding
#' @param intensity a list with elements \code{x} and \code{y} containing hybridization intensities; each
#' 	is a matrix with same dimensions and same row/column names as \code{G}
#' @param normalized logical; have intensities been normalized?
#' @param filter.sites character vector of filters attached to markers
#' @param filter.samples character vector of filters attached to samples
#' @param check logical; if \code{TRUE}, do sanity checks on input
#' @param ... ignored
#'
#' @return a new \code{genotypes} object
#'
#' @details The input matrix \code{G} *must* have row and column names to help the package keep the marker
#' 	map, sample metadata, and genotypes themselves in sync.
#' 
#' @section The \code{genotypes} class:
#' 	This class is designed to be a lightweight container for genotype data on a set of samples typed for a
#' 	panel of biallelic SNP markers on a microarray.  The object inherits from base-\code{R}'s class \code{matrix},
#' 	so any code which accepts a matrix (including the \code{apply} family) will work on a \code{genotypes} object.
#' 	
#' 	Attributes of \code{genotypes} objects include:
#'	\itemize{
#'	\item \code{map} -- marker metadata in PLINK format (chr, marker, cM, pos, A1, A2, ...)
#'	\item \code{ped} -- pedigree/sample metadata in PLINK format (individual ID, family ID,
#'		mom ID, dad ID, sex, phenotype, ...)
#'	\item \code{intensity} -- \code{list}(\code{x} = [X-intensities], \code{y} = [y-intensities])
#'	\item \code{normalized} -- have intensities been normalized?
#'	\item \code{baf} -- matrix of B-allele frequencies (BAFs; see \code{\link{tQN}})
#'	\item \code{lrr} -- matrix of log2 intensity rations (LRRs; see \code{\link{tQN}})
#'	\item \code{filter.sites} -- homage to the FILTER field in VCF format, a flag for suppresing
#'		sites (rows) in downstream analyses
#'	\item \code{filter.samples} -- same as above, but along other dimension (columns)
#'	\item \code{alleles} -- manner in which alleles are encoded: "native" (ACTGHN),
#'	"01" (allele dosage wrt ALT allele), "relative" (allele dosage wrt MINOR allele)
#' 	}
#'	All attributes are maintained "parallel" to the genotypes matrix itself, and additionally have names
#'	to avoid ambiguity.
#' 
#'	Note that missing values (NAs/NaNs) are used for no-calls, in order to take advantage of R's behaviors on missing data.
#'
#' @section The marker map:
#' A valid marker map is a required attribute of a \code{genotypes} object.  It is a dataframe with (at least)
#' the following columns, in the following order.  Columns followed by an asterisk (*) are optional but may be
#' required for some downstream operations.
#' \itemize{
#' \item \code{chr} -- (character, factor) chromosome identifier; use \code{NA} for missing
#' \item \code{marker} -- (character, factor) *globally-unique* marker name, cannot be missing
#' \item \code{cM} -- (numeric) genetic position of this marker in cM; use zero for missing
#' \item \code{pos} -- (integer) position of this marker in basepairs; use zero for missing
#' \item \code{A1}* -- (character, factor) REFERENCE allele, case-insensitive, cannot be missing
#' \item \code{A2}* -- (character, factor) ALTERNATE allele, case-insensitive, cannot be missing
#' }
#' Rownames must be present and must match the contents of column "marker".
#' 
#' @section The "pedigree":
#' Although "pedigree" is used in homage to the nomenclature of the PLINK package, this attribute simply contains
#' sample metadata even if true pedigrees are unknown.  It is a dataframe with (at least) the following columns,
#' the first 6 of which are for PLINK compatibility, in the following order.
#'  \itemize{
#' \item \code{fid} -- (character, factor) "family" ID (aka group ID); can indicate family, population, batch...
#' \item \code{iid} -- (character, factor) *globally-unique* individual ID
#' \item \code{mom} -- (character, factor) individual ID of this sample's mother; use zero for missing
#' \item \code{dad} -- (character, factor) individual ID of this sample's father; use zero for missing
#' \item \code{sex} -- (integer) 1=male, 2=female, 0=unknown/missing
#' \item \code{pheno} -- (numeric) phenotype; 0/-9=missing, 1=control, 2=case, any other values allowed
#' 	are taken to be a quantitative trait
#' }
#' Rownames must be present and must match the contents of column "iid".  The pedigree is auto-generated
#' when missing, and in that case every sample is assigned an "fid" identical to its "iid".
#'
#' @section Filters:
#' The \code{filter.*} fields are character vectors describing the filter(s), if any, with which to mark markers
#' or samples.  An empy string (\code{""}) indicates a "passing" marker or sample.  Filters are appended to the
#' filter string as single characters: \code{H} for excess heterozygosity; \code{N} for excess no-call rate;
#' \code{I} (for sampes only) for abnormal intensity pattern; \code{F} (for markers only) for abberrant allele frequency.
#'
#' @export
genotypes <- function(G, map, ped = NULL, alleles = c("auto","native","01","relative"),
					  intensity = NULL, normalized = FALSE, filter.sites = NULL, filter.samples = NULL,
					  check = TRUE, ...) {
	
	## --- genotypes matrix --- ##
	## check genotypes matrix for valid size and names
	if (!is.matrix(G))
		stop("Input is not a matrix.")
	if (is.null(dimnames(G)) || length(rownames(G)) != nrow(G) || length(colnames(G)) != ncol(G))
		stop("Input matrix is missing row and/or column names.")
	if (length(unique(colnames(G))) != ncol(G))
		stop("Column names of input matrix (eg. sample IDs) are not unique.")
	if (length(unique(rownames(G))) != nrow(G))
		stop("Row names of input matrix (eg. marker IDs) are not unique.")
	
	## --- marker map --- ##
	## quick validity check
	if (!.is.valid.map(map))
		stop("Marker map lacks the four required columns: chr, marker, cM, pos.")
	## force datatypes
	map$marker <- as.character(map$marker)
	map$cM <- as.numeric(map$cM)
	map$pos <- as.integer(map$pos)
	## check alleles
	if (ncol(map) >= 6) {
		if (all(colnames(map)[5:6] == c("A1","A2"))) {
			map$A1 <- as.character(toupper(map$A1))
			map$A2 <- as.character(toupper(map$A2))
			monomorphic <- (map$A1 == map$A2)
			if (any(monomorphic)) {
				warning( paste("REF and ALT alleles shouldn't match, but they do at these markers:\n",
							   paste(rownames(map)[ which(monomorphic) ], collapse = ",")) )
			}
		}
	}
	mk <- rownames(G)
	if (length(intersect(rownames(map), mk)) != nrow(G))
		stop("All markers in genotypes matrix should be present in marker map.")
	
	## check positions
	if (any(map$cM < 0, na.rm = TRUE))
		stop("Some genetic positions (cM) are negative; that shouldn't be possible.")
	
	## --- allele encoding --- ##
	alleles <- match.arg(alleles)
	if (check) {
		alleles.final <- "native"
		pass <- .check.encoding(G, alleles, map)
		if (!pass$result)
			stop(pass$message)
		else
			alleles.final <- pass$final
	}
	else {
		alleles.final <- alleles
	}
	
	## --- pedigree --- ##
	colnames(G) <- gsub(" ","",colnames(G))
	sm <- colnames(G)
	if (!is.null(ped)) {
		if (!.is.valid.ped(ped))
			stop(paste("The sample metadata provided is invalid; it should be a dataframe with",
					   "these columns: fid, iid, mom, dad, sex, pheno.  See ?genotypes."))
		if (length(intersect(rownames(ped), sm)) != ncol(G))
			stop("All samples in genotypes matrix should be present in marker map.")
	}
	else {
		ped <- make.fam(colnames(G))
	}
	
	## --- filters --- ##
	if (!is.null(filter.sites)) {
		if (!is.list(filter.sites) || length(filter.sites) != nrow(G))
			stop(paste("Site filters are invalid: this should be list of zero or more filters per sample",
					   "with length equal to number of rows in genotypes matrix."))
	}
	else
		filter.sites <- .init.filters(rownames(G))
	
	if (!is.null(filter.samples)) {
		if (!is.list(filter.samples) || length(filter.samples) != ncol(G))
			stop(paste("Sample filters are invalid: this should be list of zero or more filters per sample",
					   "with length equal to number of columns in genotypes matrix."))
	}
	else
		filter.samples <- .init.filters(colnames(G))
	
	## construct the new 'genotypes' object
	rez <- structure(G, map = map[ mk, ], ped = ped[ sm, ], alleles = alleles,
					 filter.sites = filter.sites[mk], filter.samples = filter.samples[sm])
	class(rez) <- c("genotypes", class(rez))
	
	## --- intensity matrices --- ##
	if (!is.null(intensity)) {
		if (is.list(intensity) && length(intensity) == 2) {
			pass <- sapply(intensity, function(x) {
				colnames(x) <- gsub(" ","", colnames(x))
				pass <- all(ncol(x) == ncol(G), nrow(x) == nrow(G))
				pass <- pass && !any(sapply(dimnames(x), is.null))
				pass <- pass && all(sm %in% colnames(x)) && all(mk %in% rownames(x))
				return(pass)
			})
			if (!all(pass))
				stop(paste("Intensity matrices are invalid: they should have same dimensions as the",
						   "genotypes matrix, and matching row and column names."))
			else
				intensity <- lapply(intensity, function(x) {
						colnames(x) <- gsub(" ","", colnames(x))
						x[mk,sm, drop = FALSE]
					})
		}
		else {
			stop("Intensity matrices should be supplied as list of length 2.")
		}
		attr(rez, "intensity") <- intensity
		attr(rez, "normalized") <- as.logical(normalized)[1]
	}
	
	## do a final check, if requested
	is.valid <- TRUE
	if (check) {
		is.valid <- validate.genotypes(rez)
	}
	
	## return
	if (is.valid)
		return(rez)
	else
		stop("The new genotypes object failed validation.")
	
	
}

## helper function to check allele encodings
.check.encoding <- function(G, alleles = c("auto","native","01","relative"), map, ...) {
	
	.reply <- function(msg, rez, ...) {
		return(list(message = msg, result = rez, ...))
	}
	
	if (!is.matrix(G))
		.reply("Input not a matrix.", FALSE)
	
	final <- NA	
	if (alleles == "auto") {
		if (is.character(G)) {
			final <- "native"
		}
		else if (is.numeric(G)) {
			if (all(c("A1","A2") %in% colnames(map)))
				final <- "01"
			else
				final <- "relative"
		}
		else {
			.reply("Can't auto-detect allele encoding.", FALSE)
		}	
	}
	else if (alleles == "01" || alleles == "relative") {
		if (!is.numeric(G))
			.reply("Allele encoding claimed doesn't match data type of genotypes matrix.", FALSE)
	}
	else if (alleles == "native") {
		if (is.character(G))
			.reply("Allele encoding claimed doesn't match data type of genotypes matrix.", FALSE)
	}
	else {
		.reply("Can't determine allele encoding", FALSE)
	}
	
	## now sanity-check the numeric encodings
	if (is.numeric(G)) {
		
		## check that all are 0/1/2/NA
		g <- G
		g[ is.na(G) ] <- 3
		if (sum(tabulate(g+1, nbins = 4)) != prod(dim(g)))
			.reply("Numeric genotypes don't seem to be encoded 0/1/2/NA.", FALSE)
	
		## check that 'relative' encoding truly relative
		if (alleles == "relative") {
			if (rowMeans(G, na.rm = TRUE) > 1) {
				warning("Minor allele isn't really minor; assuming encoding is actually REF/ALT.")
				final <- "01"
			}
		}
		
	}
	
	.reply("PASS", TRUE, final = final)
	
}
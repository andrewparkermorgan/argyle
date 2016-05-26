## io.R
## functions for import and export of data from 'genotypes' objects

#' Read genotype calls and hybridization from Illumina BeadStudio output.
#'
#' @param prefix filename prefix, without working directory: the \code{*} in \code{*_FinalReport.zip}
#' @param snps dataframe containing marker map for this array, in PLINK's \code{*.bim} format
#' 	(chromosome, marker name, cM position, bp position); rownames should be set to marker names,
#' 	and those names should match those in the BeadStudio output.
#' @param in.path directory in which to search for input files
#' @param keep.intensity should hybridization intensities be kept in addition to genotype calls?
#' @param colmap named character vector mapping column names in \code{*FinalReport} to required columns
#' 	for \code{argyle} (see Details)
#' @param verify logical; if \code{TRUE}, check that \code{FinalReport} file is of expected size
#' @param checksum logical; if \code{TRUE}, generate an md5 checksum for the result
#' @param ... ignored
#'
#' @return A \code{genotypes} object with genotype calls, marker map, sample metadata and (as requested)
#' 	intensity data.
#'
#' @details This function initializes a \code{genotypes} object from Illumina BeadStudio output. (For an
#' 	example of the format, see the files in this package's \code{data/} directory.)  The two relevant
#' 	files are \code{Sample_Map.zip} and \code{*FinalReport.zip}, which contain the sample manifest
#' 	and genotype/intensity data, respectively.  On platforms with \code{unzip} available on the
#' 	command line, files will be unzipped on the fly.  Otherwise \code{FinalReport.zip} (but not
#' 	\code{Sample_Map.zip}) must be unzipped first.  This is due to the use of \code{data.table} to
#' 	handle the usually very large genotypes file.
#' 	
#' 	Use the \code{colmap} vector to assign column names in the \code{*FinalReport} file to the required
#' 	columns for argyle.  The required columns are \code{iid} (individual ID), \code{marker} (SNP/marker name),
#' 	\code{call1} (allele 1, in the same strand as in the marker map), \code{call2} (allele 2, in the
#' 	same strand as in the marker map), \code{x} (hybridization x-intensity) and \code{y} (hybridization
#' 	y-intensity).  The default column mapping is:
#'  \itemize{
#' 	\item \code{SNP Name} = \code{marker}
#' 	\item \code{Sample ID} = \code{iid}
#' 	\item \code{Allele1 - Forward} = \code{call1}
#' 	\item \code{Allele2 - Forward} = \code{call2}
#' 	\item \code{X} = \code{x}
#' 	\item \code{Y} = \code{y}
#' 	}
#' 	Note that \code{colmap} must be a named character vector, with old column headers in the \code{names()}
#' 	and new column names in the vector itself: eg. write \code{colmap = setNames( new, old )}.  An error
#' 	will be thrown if the column mapping does not provide enough information to read the input properly.
#' 	Particular attention should be paid to the encoding of the alleles in the \code{snps} object, which
#' 	will be platform-specific.  For users of the Mouse Universal Genotyping Array series from Neogen Inc,
#' 	alleles \code{A1,A2} in \code{snps} will be on the forward strand, so columns \code{Allele * - Forward}
#' 	(not \code{Allele * - Top} or \code{Allele * - AB}) are the ones to use.
#' 	
#' 	The behavior of this function with respect to missing data in the genotypes versus the contents
#' 	of \code{snps} is asymmetric.  Markers in \code{snps} which are absent in the input files will
#' 	be present in the output, but with missing calls and intensities.  Markers in the input files
#' 	which are missing from \code{snps} will simply be dropped.  If that occurs, check that the marker
#' 	names in \code{snps} match exactly those in the input file.
#' 	
#' 	Provenance of the resulting object can be traced by checking \code{attr(,"source")}.  For the paranoid,
#' 	a timestamp and checksum are provided in \code{attr(,"timestamp")} and \code{attr(,"md5")}.
#'
#' @references
#' 	Inspiration from Dan Gatti's DOQTL package: <https://github.com/dmgatti/DOQTL/blob/master/R/extract.raw.data.R>
#'
#' @export
read.beadstudio <- function(prefix, snps, in.path = ".", keep.intensity = TRUE, colmap = NULL, verify = TRUE, checksum = TRUE, ...) {

	## stop here if marker map is not well-formed
	if (!.is.valid.map(snps)) {
		if (!all(rownames(snps) == snps$marker))
			stop(paste("Marker manifest is not well-formed.  It should follow the format of a PLINK",
					   "*.bim file: a dataframe with columns <chr,marker,cM,pos> with rownames",
					   "same as 'marker' column.  If genetic positions are unknown, set them to zero."))
	}
	
	## read files from Illumina BeadStudio
	data <- .read.illumina.raw(prefix, in.path, colmap, check.dims = verify)
	rownames(data$samples) <- gsub(" ","", rownames(data$samples))
	
	## convert to matrices using data.table's optimized code
	message("Constructing genotype matrix...")
	calls <- .raw.to.matrix(data$intens, snps, keep.map = TRUE, value.col = "call")
	if (keep.intensity) {
		message("Constructing intensity matrices...")
		x <- .raw.to.matrix(data$intens, snps, value.col = "x", keep.map = FALSE)
		y <- .raw.to.matrix(data$intens, snps, value.col = "y", keep.map = FALSE)
		## verify that shapes match
		all(dim(calls) == dim(x), dim(calls) == dim(y))
		## verify that sample names are in sync
		all(colnames(calls) == colnames(x), colnames(calls) == colnames(y))
		## verify that marker names are in sync
		all(rownames(calls) == rownames(x), rownames(calls) == rownames(y))
	}
	
	message(paste("\t", nrow(calls), "sites x", ncol(calls), "samples"))
	samples.kept <- colnames(calls)
	fam <- make.fam(samples.kept, sex = data$samples[ samples.kept,"Gender" ])
	
	## construct the return 'genotypes' object
	if (keep.intensity) {
		calls <- genotypes(.copy.matrix.noattr(calls),
						   map = attr(calls, "map"), ped = fam,
						   alleles = "native",
						   intensity = list(x = x, y = y), normalized = FALSE,
						   check = FALSE)
	}
	else {
		calls <- genotypes(.copy.matrix.noattr(calls),
						   map = attr(calls, "map"), ped = fam,
						   alleles = "native",
						   check = FALSE)
	}
	
	## make a checksum, then record file source and timestamp (which would mess up checksum comparisons)
	if (checksum)
		attr(calls, "md5") <- digest::digest(calls, algo = "md5")
	
	attr(calls, "source") <- normalizePath(file.path(in.path))
	attr(calls, "timestamp") <- Sys.time()
	
	## check that all pieces of result have matching dimensions, names, ...
	if (!validate.genotypes(calls)) {
		warning("The assembled genotypes object failed validation.  See messages for possible reasons.")
	}
	
	message("Done.")
	return(calls)
	
}

## process files from BeadStudio into a dataframe (of samples) and data.table (of calls/intensities)
.read.illumina.raw <- function(prefix, in.path = ".", colmap = NULL, check.dims = FALSE, ...) {
	
	if (length(paste(prefix, in.path)) > 1)
		stop("Only read one batch at a time.  To handle multiple batches, wrap with an lapply().")
	
	## Get the sample IDs from the Sample_Map.txt file.
	rawfile <- dir(path = in.path, pattern = "Sample_Map", full.names = TRUE)
	infile <- rawfile[ grep("zip$", rawfile) ]
	as.zip <- TRUE
	## If not found, then quit.
	if (length(infile) == 0) {
		infile <- rawfile[ grep("txt$", rawfile) ]
		if (length(infile) == 0) {
			stop(paste("No file with 'Sample_Map' in the filename was found in directory",
					   in.path))
		}
		else {
			as.zip <- FALSE
		}
	}
	samplefile <- infile
	if (as.zip)
		samplefile <- unz(infile, "Sample_Map.txt")
	message(paste("Reading sample manifest from <", infile, "> ..."))
	samples.df <- read.delim(samplefile, stringsAsFactors = FALSE)
	## handle case of duplicated IDs
	renamer <- make.unique(as.character(samples.df$Name))
	rownames(samples.df) <- renamer
	nsamp <- length(renamer)
	
	## Find a file with "FinalReport" in the filename.
	rawfile <- dir(path = in.path, pattern = "FinalReport", full.names = TRUE)
	infile <- rawfile[ grep("zip$", rawfile) ]
	as.zip <- TRUE
	## If not found, then quit.
	if (length(infile) == 0) {
		infile <- rawfile[ grep("txt$", rawfile) ]
		if (length(infile) == 0) {
			stop(paste("No file with 'FinalReport' in the filename was found in directory",
					   in.path, ".  Please make sure that the FinalReport file is",
					   "in the specified directory, and that there is exactly one."))
		}
		else {
			as.zip <- FALSE
		}
	}
	## If there is more than one FinalReport file, then quit.
	if (length(infile) > 1) {
		stop(paste("There is more than one file with FinalReport in the filename.",
				   "Please place only one data set in each directory. If you have",
				   "unzipped the file, keep the original in a different directory."))
	}
	
	## try to unzip on fly
	piper <- infile
	if (as.zip) {
		## check that system supports unzip
		piper <- paste0("unzip -ap '", infile, "'")
		rez <- system(paste(piper, "| head -n 10"), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
		if (rez) {
			## nope, system can't unzip
			stop(paste("Looks like this system doesn't support on-the-fly decompression: please unzip the",
				 		"FinalReport file manually and try again."))
		}
	}
	
	## first check header to get # of SNPs
	genofile <- infile
	if (as.zip)
		genofile <- unz(infile, gsub("\\.zip$",".txt", basename(infile)))
	header <- read.delim(genofile, header = FALSE, sep = "\t", nrows = 8, skip = 1,
						 colClasses = "character", stringsAsFactors = FALSE)
	if (!all(ncol(header) >= 2, "NUM SNPS" %in% toupper(header[,1])))
		stop("Header of FinalReport file should be key-value pairs, one of which is 'Total SNPs {x}'")
	header <- setNames( header[,2], toupper(header[,1]) )
	nsnps <- as.integer(header["NUM SNPS"])
	if (!(is.numeric(nsnps) && nsnps > 0))
		stop("Can't understand number of markers to read; header says '", nsnps, "'")
	
	## slurp file into a data.table, skipping the 9 header lines
	message(paste("Reading genotypes and intensities for", nsnps, "markers x",
				  nsamp, "samples from <", infile, "> ..."))
	data <- data.table::fread(piper, skip = 9, showProgress = interactive(), stringsAsFactors = FALSE, sep = "\t")
	
	## Construct the column-naming map
	cols.needed <- c("marker","iid","x","y","call1","call2")
	if (is.null(colmap)) {
		colmap <-  c("SNP Name" = "marker",
				   "Sample ID" = "iid",
				   "X" = "x","Y" = "y",
				   "Allele1 - Forward" = "call1",
				   "Allele2 - Forward" = "call2")
	}
	
	## Check that all required columns are specified in the map
	columns <- match(cols.needed, colmap)
	if (any(is.na(columns))) {
		stop(paste0("The column-name map must specify columns mapping to each of the following:\n",
					"'iid' (individual ID), 'marker' (SNP/marker name), 'x' (x-intensity), 'y' (y-intensity),\n",
					"''call1' (allele 1, in same strand as specified in marker map), and 'call2' (allele 2).\n\n",
					"You are missing the following: ", paste(colmap[is.na(columns)], sep = ", ")))
	}
	
	## Check that all mapped columns are present in the input
	columns <- match(names(colmap), names(data))  
	if (any(is.na(columns))) {
		stop(paste("All of the required column names were not found in the",
				   "FinalReport file. The missing column(s) are:", paste(
				   	names(colmap)[is.na(columns)], collapse = ", ")))
	}
	
	## assign new column names using column map
	data.table::setnames(data, names(colmap), colmap)
	
	## check that data is of expected size
	if (check.dims) {
		if (nrow(data) != (nsnps*nsamp))
			stop(paste0("Row count of genotype data (", nrow(data), ") doesn't match what was expected (", nsnps, " x ", nsamp, " = ", nsnps*nsamp, ").",
						" FinalReport file might be corrupt."))
	}
	else {
		# live dangerously
	}
	
	## rename samples by index
	nsamples <- length(samples)
	newids <- rep(renamer, 1, each = nsnps)
	#print(tail(cbind(data, newid = newids)))
	data.table::set(data, i = NULL, "iid", newids)
	
	## convert 2-column allele calls to single column; mark hets, missing, etc.
	data.table::set(data, i = NULL, "call", paste0(data$call1, data$call2))
	data.table::set(data, i = NULL, "is.het", (data$call1 != data$call2))
	data.table::set(data, i = NULL, "is.na", (data$call1 == "-" | data$call2 == "-"))
	data.table::set(data, i = which(data$is.het), "call", "H")
	data.table::set(data, i = which(data$is.na), "call", "N")
	data.table::set(data, i = NULL, "call", substr(data$call, 1, 1))
	
	## pre-key by marker (SNP name) for next step
	data.table::setkey(data, marker)
	
	return( list(samples = samples.df, intens = data) )
	
}

## convert data.table of calls/intensities to a (sites x samples) matrix
.raw.to.matrix <- function(data, snps, keep.map = FALSE, make.names = FALSE, verbose = FALSE,
						   sample.id.col = "iid", value.col = "call", ...) {
	
	if (!inherits(data, "data.table"))
		stop("Input should be an object of class 'data.table'.")
	
	## strip column names which might conflict between input and marker map
	if ("cM" %in% colnames(data))
		data.table::set(data, i = NULL, "cM", NULL)
	if ("chr" %in% colnames(data))
		data.table::set(data, i = NULL, "chr", NULL)
	if ("pos" %in% colnames(data))
		data.table::set(data, i = NULL, "pos", NULL)
	if (make.names)
		data.table::set(data, i = NULL, "marker", make.names(data$marker))
	
	## reshape to big matrix
	fm <- paste("marker ~", sample.id.col)
	gty.mat <- data.table::dcast.data.table(data, as.formula(fm), value.var = value.col)
	data.table::setkey(gty.mat, marker)
	
	before <- unique(gty.mat$marker)
	nsnps.before <- length(unique(gty.mat$marker))
	.map <- data.table::data.table(snps[ ,c("chr","marker","cM","pos") ])
	data.table::setkey(.map, marker)

	if (verbose)
		message(paste("Attaching map: started with", nsnps.before, "markers"))
	gty.mat <- data.table:::merge.data.table(gty.mat, .map)
	nsnps.after <- length(unique(gty.mat$marker))
	if (verbose)
		message(paste("Done attaching map: ended with", nsnps.after, "markers"))
	
	if (nsnps.before != nsnps.after && verbose) {
		if (nsnps.after - nsnps.before < 100) {
			message(paste("Dropped the following markers:"))
			message( paste(setdiff(before, gty.mat$marker) , collapse = ",") )
		}
		else {
			message(paste("\tdropped", (nsnps.after-nsnps.before), "markers."))
		}
	}
	
	## sort by position
	data.table::setorder(gty.mat, chr, pos, cM, marker)
	cn <- names(gty.mat)
	cols <- c("chr","marker","cM","pos")
	oth <- setdiff(cn, cols)
	data.table::setcolorder(gty.mat, c(cols, oth))
	
	## demote back to dataframe
	gty.mat <- as.data.frame(gty.mat)
	newmap <- gty.mat[ ,1:4, drop = FALSE ]
	rownames(newmap) <- as.character(newmap$marker)
	newmap <- data.frame(newmap, snps[ rownames(newmap),!(colnames(snps) %in% c("chr","marker","cM","pos")) ])
	gty.mat <- as.matrix(gty.mat[ ,-(1:4), drop = FALSE ])
	
	rownames(gty.mat) <- as.character(newmap$marker)
	colnames(gty.mat) <- gsub(" ","", colnames(gty.mat))
	
	if (keep.map)
		attr(gty.mat, "map") <- newmap
	
	return(gty.mat)
	
}

#' Export genotyping result in format suitable for DOQTL
#'
#' @param gty a \code{genotypes} object with intensity data attached
#' @param where name of output file, including path (else it goes in working directory)
#' @param recode if \code{TRUE}, genotype calls will be recoded 0/1/2 with respect to reference alleles
#' 	before the genotypes matrix is saved
#' @param ... ignored
#' 
#' @return Returns \code{TRUE} on completion.  The Rdata file at \code{where} contains the following
#' 	objects:
#' 	\itemize{
#' 	\code{G} -- genotype calls matrix
#' 	\code{x} -- matrix of x-intensities
#' 	\code{y} -- matrix of y-intensities
#' 	\code{sex} -- named vector of sample sexes (\code{NA} if missing)
#' 	\code{snps} -- the marker map attached to the input object
#' 	}
#' 	All matrices are sites x samples, following the convention of this pacakge, and have both row and column names.
#'
#' @references
#' DOQTL home: \url{http://cgd.jax.org/apps/doqtl/DOQTL.shtml}
#' 
#' Gatti DM et al. (2014) Quantitative trait locus mapping methods for Diversity Outbred mice. G3 4(9): 1623-1633. doi:10.1534/g3.114.013748.
#' 
#' Svenson KL et al. (2012) High-resolution genetic mapping using the mouse Diversity Outbred population. Genetics 190(2): 437-447. doi:10.1534/genetics.111.132597.
#'
#' @export
export.doqtl <- function(gty, where = "doqtl.Rdata", recode = FALSE, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.intensity(gty)))
		stop("To export stuff for DOQTL, need (1) genotype matrix; (2) intensity matrices; (3) marker map.")
	
	message("Preparing objects...")
	where <- file.path(where)
	x <- attr(gty, "intensity")$x
	y <- attr(gty, "intensity")$y
	if (!recode)
		G <- .copy.matrix.noattr(gty)
	else
		G <- .copy.matrix.noattr(recode.genotypes(gty, "01"))
	sex <- setNames( rep(NA, ncol(G)), colnames(G) )
	if (.has.valid.ped(gty) && "sex" %in% colnames(attr(gty, "ped")))
		sex[ rownames(attr(gty, "ped")) ] <- as.character(attr(gty, "ped")$sex)
	sex[ sex == "1" ] <- "M"
	sex[ sex == "2" ] <- "F"
	sex[ sex == "0" ] <- NA
	#print(sex)
	snps <- attr(gty, "map")
	message(paste("Saving DOQTL input objects in <", where, ">..."))
	save(x, y, G, sex, snps, file = where)
	
	message("Done.")
	invisible(TRUE)
	
}

#' Convert a \code{genotypes} object to an \code{R/qtl} object
#' 
#' @param gty a \code{genotypes} object
#' @param type cross type for \code{R/qtl} (only \code{"bc"} [backcross] and \code{"f2"} [F2 intercross] currently supported)
#' @param chroms vector of chromosome names to include in output (in order)
#' @param ... ignored
#' 
#' @return an object of class \code{cross}, with the specified cross type
#' 
#' @details Karl Broman's \code{R/qtl} is a widely-used package for mapping quantiative traits
#' 	in experimental crosses of laboratory organisms and crop plants.  It expects genotypes to
#' 	be coded with respect to parental lines: eg. AA, AB, BB for an F2 cross between (true-breeding)
#' 	lines A and B.  Be sure to recode genotypes in that mannyer way before passing them to this function.
#' 	
#' 	Marker positions in \code{R/qtl} are expressed in centimorgans, not basepairs, so only markers with
#' 	non-zero, non-missing genetic positions will be included in the output of this function.
#' 
#' @references
#' \code{R/qtl}: \url{http://www.rqtl.org}
#' 
#' Broman KW, Wu H, Sen S, Churchill GA. (2003) R/qtl: QTL mapping in experimental crosses.
#' 	Bioinformatics 19:889-890. doi:10.1093/bioinformatics/btg112.
#' 
#' Broman KW, Sen S. (2009) A Guide to QTL Mapping with R/qtl. Springer, New York.
#' 
#' @seealso \code{\link[qtl]{read.cross}}, \code{\link{as.genotypes.cross}} (for inverse operation)
#' 
#' @export as.rqtl
as.rqtl <- function(gty, type = c("f2","bc"), chroms = paste0("chr", c(1:19,"X")), ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.map(gty) && .has.valid.ped(gty)))
		stop("Please supply an object of class 'genotypes' with valid marker map and sample metadata.")
	
	if (!(is.numeric(gty) && attr(gty, "alleles") == "parent"))
		warning(paste("For export to R/qtl, genotypes should be encoded numerically, and by reference to",
					  "the parental strains of a cross.  See ?recode.to.parent."))
	
	## dump intensity data
	gty <- drop.intensity(gty)
	
	## drop unfamiliar chromosomes and positionless markers
	gty <- subset(gty, chr %in% chroms & !is.na(cM) & cM > 0)
	map <- attr(gty, "map")
	map$chr <- factor(map$chr, levels = chroms)
	map <- droplevels(map)
	
	message(paste("Exporting genotypes at", nrow(gty), "markers on",
				  length(unique(map$chr)), "chromosomes."))
	
	## make sure sex is right
	sex <- .fix.sex(attr(gty, "ped")$sex)
	sex <- factor(c("male","female",NA)[ sex+1 ],
				  levels = c("female","male"))
	
	.make.rqtl <- function(cc) {
		g <- .copy.matrix.noattr(subset(gty, chr == cc))
		g <- matrix(t(g)+1,
					nrow = ncol(g), ncol = nrow(g),
					dimnames = list(colnames(g), rownames(g)))
		m <- subset(map, chr == cc)
		this.chr <- list(data = g,
			 			map = setNames( as.vector(m$cM), as.character(m$marker) ))
		if (grepl("X", cc))
			class(this.chr) <- "X"
		else
			class(this.chr) <- "A"
		return(this.chr)
	}
	
	## loop on chromosomes
	message("Converting genotypes...")
	geno <- lapply(levels(map$chr), .make.rqtl)
	names(geno) <- gsub("^chr","", levels(map$chr))
	rez <- list(geno = geno)
	
	## convert PLINK to R/qtl style of phenotype definition
	pheno <- attr(gty, "ped")
	pheno$pheno[ pheno$pheno == -9 ] <- NA
	newsex <- pheno$sex
	newsex[ newsex == 0 ] <- NA
	newsex[ newsex == 2 ] <- 0
	pheno$sex <- newsex
	
	rez$pheno <- pheno[ ,c("pheno","sex",setdiff(colnames(pheno), c("pheno","sex"))) ]
	class(rez) <- c("f2","cross")
	attr(rez, "alleles") <- c("A","B")
	
	message("Done.")
	return(rez)
	
	
}

#' Convert an \code{R/qtl} object to a \code{genotypes} object
#' 
#' @param x a \code{qtl::cross} object
#' @param ... ignored
#' 
#' @return an object of class \code{genotypes}
#' 
#' @details Karl Broman's \code{R/qtl} is a widely-used package for mapping quantiative traits
#' 	in experimental crosses of laboratory organisms and crop plants.  It expects genotypes to
#' 	be coded with respect to parental lines: eg. AA, AB, BB for an F2 cross between (true-breeding)
#' 	lines A and B.  Be sure to recode genotypes in that mannyer way before passing them to this function.
#' 	
#' 	Marker positions in \code{R/qtl} are expressed in centimorgans, not basepairs.  On conversion,
#' 	physical positions are faked by assuming recombination rate 1 cM per 1 Mbp and rounding to
#' 	the next-lowest integer.
#' 
#' 	Only crosses of type \code{"f2"} (F2 intercross) or \code{"bc"} are supported, and partially-
#' 	informative genotypes will probably be mangled.
#' 
#' @references
#' \code{R/qtl}: \url{http://www.rqtl.org}
#' 
#' Broman KW, Wu H, Sen S, Churchill GA. (2003) R/qtl: QTL mapping in experimental crosses.
#' 	Bioinformatics 19:889-890. doi:10.1093/bioinformatics/btg112.
#' 
#' Broman KW, Sen S. (2009) A Guide to QTL Mapping with R/qtl. Springer, New York.
#' 
#' @seealso \code{\link[qtl]{read.cross}}, \code{\link{as.rqtl.genotypes}} (for inverse operation)
#' 
#' @export as.genotypes
as.genotypes <- function(x, ...) {
	
	if (!inherits(x, "cross") && any(inherits(x, "f2"), inherits(x, "bc")))
		stop("Please supply an object of class 'cross' (from KW Broman's R/qtl package).")
	
	geno <- do.call( rbind, lapply(x$geno, function(chr) t(chr$data)) )
	
	geno <- geno-1
	map <- do.call( rbind, lapply(names(x$geno), function(chr) data.frame(chr = chr, marker = names(x$geno[[chr]]$map),
																		  cM = x$geno[[chr]]$map, pos = floor(x$geno[[chr]]$map*1e6))) )
	map$A1 <- attr(x, "alleles")[1]
	map$A2 <- attr(x, "alleles")[2]
	if (!any(grepl("^chr", map$chr)))
		map$chr <- paste0("chr", map$chr)
	
	fam <- make.fam(rownames(x$pheno))
	if (!is.null(x$pheno$sex))
		fam$sex <- .fix.sex(x$pheno$sex)
	for (col in setdiff(colnames(x), c("sex","pgm","id","ID")))
		fam[ ,col ] <- x$pheno[ ,col ]
	fam$pheno <- x$pheno[,1]
	colnames(geno) <- as.character(rownames(fam))
	
	rez <- genotypes(geno, map = map, ped = fam, alleles = "01")
	return(rez)
	
}

#' Export genotypes in Stanford HGDP format
#' 
#' @param gty a \code{genotypes} object
#' @param prefix filename prefix for output; result will be two files, \code{{prefix}.geno} and
#' 	with genotypes and \code{{prefix}.map} with marker map.
#' @param ... ignored
#' 
#' @details Write genotypes to disk in the Stanford HGDP format, which can be read by (among
#' 	others) the PGDSpider format-conversion suite.
#' 
#' @references
#' Lischer HEL and Excoffier L (2012) PGDSpider: An automated data conversion tool for connecting
#' 	population genetics and genomics programs. Bioinformatics 28: 298-299.
#' 
#' @export write.hgdp
write.hgdp <- function(gty, prefix, ...) {
	
	if (!(inherits(gty, "genotypes") && .has.valid.map(gty)))
		stop("Please supply an object of class 'genotypes' which includes a valid marker map.")
	
	## convert genotypes to numeric
	gty <- recode.genotypes(gty, "01")
	
	## now convert to HGDP style (AA, AB, BB, -)
	message("Converting genotypes to HGDP encoding (AA/AB/BB/-)...")
	out <- matrix("-", nrow = nrow(gty), ncol = ncol(gty),
					 dimnames = dimnames(gty))
	for (i in seq_len(ncol(gty))) {
		majr <- with(attr(gty, "map"), paste0(A1, A1))
		hetr <- with(attr(gty, "map"), paste0(A1, A2))
		minr <- with(attr(gty, "map"), paste0(A2, A2))
		is.miss <- is.na(gty[ ,i ])
		is.maj <- gty[ ,i ] == 0 & !is.miss
		is.het <- gty[ ,i ] == 1 & !is.miss
		is.min <- gty[ ,i ] == 2 & !is.miss
		out[ is.maj,i ] <- majr[is.maj]
		out[ is.het,i ] <- hetr[is.het]
		out[ is.min,i ] <- minr[is.min]
		out[ is.miss,i ] <- "-"
	}
	colnames(out)[1] <- paste0("\t", colnames(out)[1])
	write.table(out, paste0(prefix, ".geno"), col.names = TRUE, row.names = TRUE,
				quote = FALSE, sep = "\t")
	
	## prepare marker map
	message("Writing marker map...")
	write.table(attr(gty, "map")[ ,c("marker","chr","pos") ], paste0(prefix, ".map"),
				col.names = FALSE, row.names = FALSE,
				quote = FALSE, sep = "\t")
	
	message("Done.")
	invisible(TRUE)
	
}

#' Convert genotypes to a dataframe
#' 
#' @param gty a \code{genotypes} object
#' @param ... ignored
#' 
#' @return a \code{data.frame} with marker information in the leftmost columns, followed by a
#' 	matrix of genotypes with samples in columns
#' 	
#' @details In general the dataframe will be a less-efficient way to store genotypes, but is a
#' 	useful intermediate for writing genotypes to disk in a human-readable format.
#' 
#' @export
as.data.frame.genotypes <- function(gty, ...) {

	if (!(inherits(gty, "genotypes") && .has.valid.map(gty)))
		stop("Please supply an object of class 'genotypes' which includes a valid marker map.")
		
	map <- attr(gty, "map")[ ,c("chr","pos","cM","marker","A1","A2") ]
	geno <- .copy.matrix.noattr(gty)
	
	df <- as.data.frame(geno)
	colnames(df) <- colnames(geno)
	df <- cbind(map, df)
	return(df)

}
#as.data.frame <- function(x, ...) UseMethod("as.data.frame")
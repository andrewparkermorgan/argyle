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
#' 	The behavior of this function with respect to missing data in the genotypes versus the contents
#' 	of \code{snps} is asymmetric.  Markers in \code{snps} which are absent in the input files will
#' 	be present in the output, but with missing calls and intensities.  Markers in the input files
#' 	which are missing from \code{snps} will simply be dropped.  If that occurs, check that the marker
#' 	names in \code{snps} match exactly those in the input file.
#'
#' @references
#' 	Inspiration from Dan Gatti's DOQTL package: <https://github.com/dmgatti/DOQTL/blob/master/R/extract.raw.data.R>
#'
#' @export
read.beadstudio <- function(prefix, snps, in.path = ".", keep.intensity = TRUE, ...) {

	## stop here if marker map is not well-formed
	if (!.is.valid.map(snps)) {
		if (!all(rownames(snps) == snps$marker))
			stop(paste("Marker manifest is not well-formed.  It should follow the format of a PLINK",
					   "*.bim file: a dataframe with columns <chr,marker,cM,pos> with rownames",
					   "same as 'marker' column.  If genetic positions are unknown, set them to zero."))
	}
	
	## read files from Illumina BeadStudio
	data <- .read.illumina.raw(prefix, in.path)
	rownames(data$samples) <- gsub(" ","", rownames(data$samples))
	
	## convert to matrices using data.table's optimized code
	message("Constructing genotype matrix...")
	calls <- .raw.to.matrix(data$intens, snps, keep.map = TRUE, value.col = "call")
	if (keep.intensity) {
		message("Constructing intensity matrices...")
		x <- .raw.to.matrix(data$intens, snps, value.col = "x")
		y <- .raw.to.matrix(data$intens, snps, value.col = "y")
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
	
	class(calls) <- c("genotypes", class(calls))
	attr(calls, "ped") <- fam
	if (keep.intensity) {
		attr(calls, "intensity") <- list(x = x, y = y)
		attr(calls, "normalized") <- FALSE
	}
	attr(calls, "alleles") <- "native"
	attr(calls, "filter.samples") <- rep(FALSE, ncol(calls))
	attr(calls, "filter.sites") <- rep(FALSE, nrow(calls))
	
	message("Done.")
	return(calls)
	
}

.read.illumina.raw <- function(prefix, in.path = ".", ...) {
	
	if (length(paste(prefix, in.path)) > 1)
		stop("Only read one batch at a time.  To handle multiple batches, wrap with an lapply().")
	
	## Get the sample IDs from the Sample_Map.txt file.
	samplefile <- dir(path = in.path, pattern = "Sample_Map.zip", full.names = TRUE)
	## If not found, then quit.
	if(length(samplefile) == 0) {
		stop(paste("No file called 'Sample_Map.txt' was found in directory",
				   in.path[i], ".  Please make sure that the Sample_Map file is unzipped and",
				   "in the specified directory."))
	}
	message(paste("Reading sample manifest from <", samplefile, "> ..."))
	samples.df <- read.delim(unz(samplefile, "Sample_Map.txt"), stringsAsFactors = FALSE)
	rownames(samples.df) <- as.character(samples.df$Name)
	
	## Find a file with "FinalReport" in the filename.
	rawfile <- dir(path = in.path, pattern = "FinalReport", full.names = TRUE)
	infile <- rawfile[grep("zip$", rawfile)]
	as.zip <- TRUE
	## If not found, then quit.
	if (length(infile) == 0) {
		infile <- rawfile[grep("txt$", rawfile)]
		if (length(infile) == 0) {
			stop(paste("No file with 'FinalReport' in the filename was found in directory",
					   in.path[i], ".  Please make sure that the FinalReport file is",
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
	
	## try to read file without having to unzip; fallback to plain sample name
	piper <- infile
	if (as.zip) {
		## check that system supports unzip
		piper <- paste("unzip -ap", infile)
		rez <- system(paste(piper, "| head -n 10"), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
		if (rez) {
			## nope, system can't unzip
			stop(paste("Looks like this system doesn't support on-the-fly decompression: please unzip the",
				 		"FinalReport file manually and try again."))
		}
	}
	
	## slurp file into a data.table, skipping the 9 header lines
	message(paste("Reading genotypes and intensities from <", infile, "> ..."))
	data <- fread(piper, skip = 9)
	
	## Verify that we have all of the column names that we expect.
	column.names <- c("SNP Name", "Sample ID", "X", "Y", "Allele1 - Forward",
					  "Allele2 - Forward")
	columns <- match(column.names, names(data))  
	if (any(is.na(columns))) {
		stop(paste("All of the expected column names were not found in the",
				   "FinalReport file. The missing column(s) are:", paste(
				   	colnames[is.na(columns)], collapse = ",")))
	}
	nsamples <- length(samples)
	setnames(data, c("marker","iid","call1","call2","x","y","gc","theta","x.raw","y.raw","R"))
	data[ , call := paste0(call1, call2) ]
	data[ , is.het := (call1 != call2) ]
	data[ , is.na := (call1 == "-" | call2 == "-") ]
	setkey(data, is.het)
	data[ is.het == TRUE, call := "H" ]
	setkey(data, is.na)
	data[ is.na == TRUE, call := "N" ]
	data[ , call := substr(call, 1, 1) ]
	setkey(data, marker)
	
	return( list(samples = samples.df, intens = data) )
	
}

## convert data.table of calls/intensities to a (sites x samples) matrix
.raw.to.matrix <- function(data, snps, keep.map = FALSE, make.names = FALSE, verbose = FALSE,
						   sample.id.col = "iid", value.col = "call", ...) {
	
	if (!inherits(data, "data.table"))
		stop("Input should be an object of class 'data.table'.")
	
	## strip column names which might conflict between input and marker map
	if ("cM" %in% colnames(data))
		data[ ,cM := NULL ]
	if ("chr" %in% colnames(data))
		data[ ,chr := NULL ]
	if ("pos" %in% colnames(data))
		data[ ,pos := NULL ]
	if (make.names)
		data$marker <- make.names(data$marker)
	
	## reshape to big matrix
	fm <- paste("marker ~", sample.id.col)
	gty.mat <- dcast.data.table(data, as.formula(fm), value.var = value.col)
	setkey(gty.mat, "marker")
	
	before <- unique(gty.mat$marker)
	nsnps.before <- length(unique(gty.mat$marker))
	.map <- data.table(snps[ ,c("chr","marker","cM","pos") ])
	setkey(.map, "marker")
	if (verbose)
		message(paste("Attaching map: started with", nsnps.before, "markers"))
	gty.mat <- gty.mat[.map]
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
	setorder(gty.mat, chr, pos, cM, marker)
	cn <- names(gty.mat)
	cols <- c("chr","marker","cM","pos")
	oth <- setdiff(cn, cols)
	setcolorder(gty.mat, c(cols, oth))
	
	## demote back to dataframe
	gty.mat <- as.data.frame(gty.mat)
	newmap <- gty.mat[ ,1:4 ]
	rownames(newmap) <- as.character(newmap$marker)
	gty.mat <- as.matrix(gty.mat[ ,-(1:4) ])
	
	rownames(gty.mat) <- as.character(newmap$marker)
	colnames(gty.mat) <- gsub(" ","", colnames(gty.mat))
	
	if (keep.map)
		attr(gty.mat, "map") <- newmap
	
	return(gty.mat)
	
}
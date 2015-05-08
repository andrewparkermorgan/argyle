## io.R

## read calls, intensities, samples from Illumina BeadStudio output
## taken from Dan Gatti's DOQTL package: <https://github.com/dmgatti/DOQTL/blob/master/R/extract.raw.data.R>
read.beadstudio <- function(prefix, snps, in.path = ".", keep.intensity = TRUE, ...) {

	if (any(is.na(as.numeric(rownames(snps)))))
		stop(paste("It looks like the <snps> object has bogus rownames. Rownames should match the 'marker' column,"
				   "and the marker names used on the array."))
	
	if (length(paste(prefix, in.path)) > 1)
		stop("Only read one batch at a time.  To handle multiple batches, wrap with an lapply().")
		
	## Get the sample IDs from the Sample_Map.txt file.
	samplefile <- dir(path = in.path, pattern = "Sample_Map.txt", full.names = TRUE)
	message(paste("Reading sample manifest from <", samplefile, "> ..."))
	## If not found, then quit.
	if(length(samplefile) == 0) {
		stop(paste("No file called 'Sample_Map.txt' was found in directory",
				   in.path[i], ".  Please make sure that the Sample_Map file is unzipped and",
				   "in the specified directory."))
	}
	samples.df <- read.delim(samplefile, stringsAsFactors = FALSE)
	samples <- samples.df$Name
	sexes <- samples.df$Gender
	## Find a file with "FinalReport" in the filename.
	rawfile <- dir(path = in.path, pattern = "FinalReport", full.names = TRUE)
	rawfile <- rawfile[grep("txt", rawfile)]
	## If not found, then quit.
	if (length(rawfile) == 0) {
		stop(paste("No file with 'FinalReport' in the filename was found in directory",
				   in.path[i], ".  Please make sure that the FinalReport file is unzipped and",
				   "in the specified directory."))
	}
	## If there is more than one FinalReport file, then quit.
	if (length(rawfile) > 1) {
		stop(paste("There is more than one file with FinalReport in the filename.",
				   "Please place only one data set in each directory."))
	}
	
	## Read in the first sample.  The current format requires us to skip 9 lines
	## because there is no comment delimiter at the top of the file.
	message(paste("Reading genotypes and intensities from <", rawfile, "> ..."))
	rawfile <- file(rawfile, open = "r")
	data <- readLines(con = rawfile, n = 10)
	hdr <- strsplit(data, split = "\t")
	cn <- hdr[[10]]
	
	## Verify that we have all of the column names that we expect.
	column.names = c("SNP Name", "Sample ID", "X", "Y", "Allele1 - Forward",
					 "Allele2 - Forward")
	columns <- match(column.names, cn)  
	if (any(is.na(columns))) {
		stop(paste("All of the expected column names were not found in the",
				   "FinalReport file. The missing column(s) are:", paste(
				   	colnames[is.na(columns)], collapse = ",")))
	}
	nsamples <- length(samples)
	cr <- rep(0, nsamples)
	samples.in.data <- rep("", nsamples)
	
	## Pre-allocate genotype and intensity matrices.
	## This eats memory but runs faster than repeated *bind()-ing.
	geno <- matrix(NA, nrow = nsnps, ncol = nsamples)
	if (keep.intensity) {
		x <- matrix(0, nrow = nsnps, ncol = nsamples)
		y <- matrix(0, nrow = nsnps, ncol = nsamples)
	}
	
	markers <- rownames(snps)
	
	## Read in the data per sample.
	for(j in seq_along(samples)) {
		
		message("\tsample ", j, "...")
		data <- read.table(rawfile, nrows = nsnps, as.is = TRUE)
		colnames(data) <- c("marker","id","x","y","a1","a2")
		
		## convert intensities
		this.x <- as.double(data$x)
		this.y <- as.double(data$y)
		
		## convert genotype calls
		data$a1 <- toupper(data$a1)
		data$a2 <- toupper(data$a2)
		this.call <- paste0(data$a1, data$a2)
		this.call[ data$a1 == data$a2 ] <- data$a1
		this.call[ data$a1 != data$a2 ] <- "H"
		this.call[ this.call == "--" ] <- NA
		
		## reorder markers by order of snps object
		mk <- match(data$marker, rownames(snps))
		mk <- mk[ !is.na(mk) ]
		
		if (length(mk) != length(data$markers))
			warning("Markers in the input and in the array manifest don't match.")
		
		## update matrices
		geno[ ,j ] <- this.call[mk]
		if (keep.intensity) {
			x[ ,j ] <- this.x[mk]
			y[ ,j ] <- this.y[mk]
		}
		
		samples.in.data[j] <- data$id
		
	}
	close(rawfile)
	
	message("Adding sample and marker metadata...")
	colnames(geno) <- samples.in.data
	rownames(geno) <- markers[mk]
	if (keep.intensity) {
		colnames(x) <- colnames(y) <- colnames(geno) 
		rownames(x) <- rownames(y) <- rownames(geno)
	}
	fam <- make.fam(samples, sex = sexes)
	
	## sort genotypes by karyotype order
	o <- with(snps, order(chr, pos, cM))
	snps <- snps[ o, ]
	geno <- geno[ o, ]
	if (keep.intensity) {
		x <- x[ o, ]
		y <- y[ o, ]
	}
	
	class(geno) <- c("genotypes", class(geno))
	attr(geno, "map") <- snps
	attr(geno, "ped") <- fam
	if (keep.intensity) {
		attr(geno, "intensity") <- list(x = x, y = y)
		attr(geno, "normalized") <- FALSE
	}
	attr(geno, "alleles") <- "native"
	attr(geno, "filter.samples") <- rep(FALSE, ncol(geno))
	attr(geno, "filter.sites") <- rep(FALSE, nrow(geno))
	
	message("Done.")
	return(geno)
	
}
## plink.R
## wrappers for interacting with plink

#' Read a PLINK binary fileset into a \code{genotypes} object
#' 
#' @param prefix path to a PLINK fileset, excluding \code{*.bed} suffix
#' 
#' @return a \code{genotypes} object, with alleles in the "01" encoding (see \code{\link{recode.genotypes}})
#' 
#' @details Reads a PLINK binary fileset using pure \code{R}.  The fileset should be generated in 
#' 	the "variant-major" format (the default in all recent versions of PLINK.) Missing genotypes are
#' 	represented as \code{NA}s.
#' 	
#' @seealso \code{\link{write.plink}}
#' 
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @export
read.plink <- function(prefix, ...) {
	
	## construct filenames; check that all are present and accounted for
	prefix <- gsub("\\.bed$", "", prefix)
	bed <- paste0(prefix, ".bed")
	bim <- paste0(prefix, ".bim")
	fam <- paste0(prefix, ".fam")
	
	if (!all(sapply(c(bed, bim, fam), file.exists)))
		stop("One or more of the bed/bim/fam files for this plink fileset doesn't exist.")
	
	message(paste0("Reading family info from: <", fam,">"))
	message(paste0("Reading marker info from: <", bim,">"))
	message(paste0("Reading binary genotypes from: <", bed,">"))
	
	## read map and family file to calculate number of markers and samples
	bim <- read.map(bim)
	fam <- read.fam(fam)
	
	nr <- nrow(bim)
	ni <- nrow(fam)
	## block size in bytes: (number of individuals)/4, to nearest byte
	bsz <- ceiling(ni/4)
	
	## open bed file and check its magic number
	bed <- file(bed, open = "rb")
	magic <- readBin(bed, "raw", 3)
	if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01"))
		stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
	
	## now actually read genotypes block by block
	intr <- interactive()
	if (intr)
		pb <- txtProgressBar(min = 0, max = nr, style = 3)
	gty <- matrix(NA, nrow = ni, ncol = nr)
	for (i in 1:nr) {
		geno.raw <- as.logical(rawToBits(readBin(bed, "raw", bsz)))
		j <- seq(1,length(geno.raw),2)
		## express genotypes as minor allele dosage (0,1,2)
		geno <- geno.raw[ j ] + geno.raw[ j+1 ]
		## recall that 0/1 is het, but 1/0 is missing
		geno[ geno.raw[ j ] == 1 & geno.raw[ j+1 ] == 0 ] <- NA
		gty[ ,i ] <- geno[1:ni]
		if (intr)
			setTxtProgressBar(pb, i)
	}
	close(bed)
	
	## add rownames (markers) and colnames (samples)
	gty <- t(gty)
	rownames(gty) <- as.character(bim$marker)
	colnames(gty) <- as.character(fam$iid)
	
	## make it a 'genotypes' object
	gty <- genotypes(gty, map = bim, ped = fam, alleles = "01")	
	return(gty)
	
}

#' Write a \code{genotypes} object as a PLINK binary fileset
#' 
#' @param gty an \code{genotypes} object
#' @param prefix path to the PLINK fileset to generate, excluding \code{*.bed} suffix
#' @param map a valid marker map; overrides existing map
#' @param fam sample metadata to override any existing metadata
#' 
#' @return \code{TRUE} on completion
#' 
#' @details Writes a PLINK binary fileset using pure \code{R}.  The fileset is generated in the
#' 	"variant-major" format (the default in all recent versions of PLINK), with appropriate magic
#' 	number.  Complete sample metadata and a complete marker map, including reference alleles (columns "A1"
#' 	and "A2"), are both required.
#' 	
#' 	This function only supports writing from character (allele enconding \code{"native"}) genotypes,
#' 	in order to avoid ambiguity around the correspondence between alleles and numeric genotypes.  But
#' 	numeric genotypes can be converted back to character via \code{recode.genotypes(,"native")}.
#' 	
#' @seealso \code{\link{read.plink}}
#' 
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @export
write.plink <- function(gty, prefix, map = NULL, fam = NULL, ...) {
	
	if (!inherits(gty, "genotypes"))
		warning("Input has not been blessed as genotype data; proceeding with skepticism.")
	
	if (!is.matrix(gty))
		gty <- as.matrix(gty)
	
	if (!is.character(gty))
		stop("Only character genotypes allowed for now, to avoid problems with allele coding.")
	
	if (is.null(map))
		if (!is.null(attr(gty, "map")))
			map <- attr(gty, "map")
	else
		stop("Must supply a marker map (columns: chr, marker, cM, pos, A1, A2).")
	
	if (is.null(fam))
		if (!is.null(attr(gty, "ped")))
			fam <- attr(gty, "ped")
	else
		stop("Must supply a family file (columns: fid, iid, mom, dad, sex, pheno).")
	
	fam.cols <- c("fid","iid","mom","dad","sex","pheno")
	if (!all(fam.cols %in% colnames(fam)))
		stop("Must supply a family file (columns: fid, iid, mom, dad, sex, pheno).")
	fam <- fam[ ,c(fam.cols, setdiff(colnames(fam), fam.cols)) ]
	
	map[,1] <- convert.names(map[,1], to = "plink")
	map[is.na(map[,3]),3] <- 0
	map[is.na(map[,4]),4] <- 0
	
	o <- with(map, order(chr, pos, cM))
	map <- map[o,]
	gty <- gty[o,]
	
	if (nrow(gty) != nrow(map))
		stop(paste0("Dimensions of genotype matrix (", nrow(gty) ,") and marker map (", nrow(map) ,") don't match."))
	if (ncol(gty) != nrow(fam))
		stop(paste0("Dimensions of genotype matrix (", ncol(gty) ,") and family file (", nrow(fam) ,") don't match."))
	
	message("Creating plink fileset <", prefix ,".*>...")
	message(paste0("\twriting family file (", ncol(gty)," individuals) ..."))
	write.plink.file(fam, paste0(prefix, ".fam"))
	message(paste0("\twriting marker map (", nrow(gty)," markers) ..."))
	write.plink.file(map, paste0(prefix, ".bim"))
	
	message(paste0("\tpreparing to convert genotypes to binary format..."))
	nr <- nrow(gty)
	ni <- ncol(gty)
	## block size in bytes: (number of individuals)/4, to nearest byte
	bsz <- ceiling(ni/4)
	
	## open bed file and write magic number
	bed <- file(paste0(prefix, ".bed"), open = "wb")
	magic <- packBits( c( rev(as.logical(c(0,1,1,0,1,1,0,0))),
						  rev(as.logical(c(0,0,0,1,1,0,1,1))),
						  rev(as.logical(c(0,0,0,0,0,0,0,1)))) )
	writeBin(magic, bed)
	
	intr <- interactive()
	if (intr)
		pb <- txtProgressBar(min = 0, max = nrow(gty), style = 3)
	for (i in seq_len(nrow(gty))) {
		row <- gty[ i, ]
		row[ is.na(row) ] <- "N"
		a1 <- row == map[i,5]
		a2 <- row == map[i,6]
		het <- row == "H"
		miss <- row == "N" | is.na(row)
		row2 <- logical(2*length(row))
		if (any(a2)) {
			j <- which(a2)-1
			row2[ 2*j + 1 ] <- TRUE
			row2[ 2*j + 2 ] <- TRUE
		}
		if (any(miss))
			row2[ 2*(which(miss)-1) + c(1) ] <- TRUE
		if (any(het))
			row2[ 2*(which(het)-1) + c(2) ] <- TRUE
		towrite <- logical(8*bsz)
		towrite[ 1:length(row2) ] <- row2
		writeBin(packBits(towrite, "raw"), bed)
		if (intr)
			setTxtProgressBar(pb, i)
	}
	close(bed)
	
	return(TRUE)
	
}

### utilities and internals

## convert plink-style to UCSC-style chromosome identifiers
unplink.chroms <- function(x, ...) {
	
	mm <- paste0("chr", c(1:19,"X","Y","M"))
	chrs <- setNames( mm, c(1:19,23,24,26) )
	chrs <- c(chrs, "X" = "chrX", "Y" = "chrY", "MT" = "chrM")
	factor(chrs[ as.character(x) ], levels = mm)
	
}

## create a fam-formatted dataframe for plink
make.fam <- function(ids, fid = NULL, sex = 0, pheno = -9, ...) {
	
	if (is.character(ids)) {
		resex <- rep(0, length(ids))
		resex[ grepl("^[mM]", sex) ] <- 1
		resex[ grepl("^[fF]", sex) ] <- 2
		resex[ is.na(resex) ] <- 0
	}
	else {
		resex <- factor(sex, levels = c(1,2))
		resex <- as.numeric(resex)
		resex[ is.na(resex) ] <- 0
	}
	
	pheno[ is.na(pheno) ] <- -9
	ids <- gsub(" ","", as.character(ids))
	if (is.null(fid))
		fid <- ids
	if (length(ids) != length(unique(ids)))
		warning("IDs are not unique; plink might complain.")
	
	fam <- data.frame(fid = fid, iid = ids,
					  mom = 0, dad = 0, sex = resex, pheno = pheno)
	rownames(fam) <- as.character(fam$iid)
	return(fam)
	
}

## shortcut to write tab-separated output files without rownames, colnames, quotes
write.plink.file <- function(x, ..., fix.missing = NA) {
	
	## convert NAs to specific missing-value code for plink
	if (!is.na(fix.missing)) {
		for (i in colnames(x))
			x[ is.na(x[,i]),i ] <- fix.missing
	}
	
	write.table(x, ..., col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

#' Create a pointer to a PLINK fileset
#' 
#' @param gty a length-1 character vector containing the path to an existing PLINK fileset, OR a
#' 	\code{genotypes} object which will be converted to one
#' @param map if \code{gty} a \code{genotypes} object, a marker map to override its current one
#' @param map if \code{gty} a \code{genotypes} object, sample metadata to override any existing metadata
#' @param where if \code{gty}, directory in which to write the new PLINK fileset
#' @param prefix basename of the fileset to be created
#' @param flags (ignored)
#' 
#' @return a length-1 character vector with path to a PLINK binary fileset, of class \code{"plink"}
#' 
#' @details Command-line PLINK takes its input from filesets on disk. This function creates a pointer to such
#' 	a fileset after checking that all its members (\code{*.bed}, \code{*.bim}, \code{*.fam}) are present.
#' 	Wrapper functions for various PLINK analyses accept this pointer as input.
#' 	
#' 	When the fileset alraedy exists, an additional attribute (\code{"where"}) is added to the result
#' 	specifying where to place the output of any future PLINK calls.  By default such files are routed
#' 	to the temporary directory corresponding to the current \code{R} session, to avoid cluttering the
#' 	directory containing the source fileset.
#' 	
#'	When the fileset does not exist and \code{gty} is a \code{genotypes} object, a new fileset is created
#'	in the location specified by \code{where} via a call to \code{write.plink(gty)}.  Then the path to
#'	that fileset is returned. 
#'
#' @export
plinkify <- function(gty, map = NULL, ped = NULL, where = tempdir(), prefix = "stuff", flags = "", ...) {
	
	## allow shortcut to bless paths of existing plink files
	if (inherits(gty, "character") && is.atomic(gty)) {
		gty <- gsub("\\.bed$", "", gty)
		paths <- sapply(c("bim","bed","fam"), function(x) file.exists(paste0(gty, ".", x)))
		if (all(paths)) {
			prefix <- as.character(gty)
			attr(prefix, "working") <- where
			class(prefix) <- c("plink", class(prefix))
			return(prefix)
		}
		else {
			stop("Can't find a plink fileset by this name in specified location.")
		}
	}
	
	## only move forward if input is in right form
	if (!inherits(gty, "genotypes"))
		stop("Plinkification of non-'genotypes' objects is not yet implemented.")
	
	nind <- ncol(gty)
	if (is.null(ped))
		ped <- make.fam(1:nind)
	
	prefix <- file.path(where, prefix)
	system(paste0("rm ", prefix, "*"))
	write.plink(gty, paste0(prefix), fam = ped)
	
	class(prefix) <- c("plink", class(prefix))
	return(prefix)
	
}

#' @export
summary.plink <- function(prefix, ...) {
	
	cat("-- Pointer to a PLINK fileset -- ","\n")
	cat("Source:", normalizePath(paste0(prefix, ".bed")), "\n")
	cat("Ouput dir:", normalizePath(attr(prefix, "working")), "\n")
	
}

#' @export
print.plink <- function(prefix, ...) {
	summary.plink(prefix, ...)
}

## read a plink fam/tfam file
read.fam <- function(f, pheno.missing = c(0,-9), ...) {
	
	ff <- f
	if (!file.exists(ff)) {
		ff <- paste0(f, ".fam")
		if (!file.exists(ff)) {
			ff <- paste0(f, ".tfam")
			if (!file.exists(ff))
				stop("Oops -- can't find a plink family file by that name.")
		}
	}
	fam <- read.table(ff, header = FALSE)
	colnames(fam) <- c("fid","iid","mom","dad","sex","pheno")
	rownames(fam) <- as.character(fam$iid)
	fam$iid <- as.character(fam$iid)
	
	## convert plink's missing-data characters to proper NAs
	fam$pheno[ fam$pheno %in% c(pheno.missing) ] <- NA
	
	return(fam)
	
}

## read plink bim/map file
read.map <- function(f, ...) {
	
	ff <- f
	if (!file.exists(ff)) {
		ff <- paste0(f, ".bim")
		if (!file.exists(ff)) {
			ff <- paste0(f, ".map")
			if (!file.exists(ff))
				stop("Oops -- can't find a plink bim/map file by that name.")
		}
	}
	fam <- read.table(ff, header = FALSE)
	colnames(fam) <- c("chr","marker","cM","pos","A1","A2")
	fam$chr <- unplink.chroms(fam$chr)
	fam$cM[ fam$cM == 0 ] <- NA
	rownames(fam) <- as.character(fam$marker)
	
	return(fam)
	
}


### wrappers for specific plink commands

## run a plink command; check that expected output is present
.plink.command <- function(prefix, command, expected = list(), where = NULL,
						   suffix = NULL, overwrite = TRUE, capture = FALSE, 
						   remove = NULL, keep = NULL, flags = "", ...) {
	
	## require that <prefix> has been blessed as a plink fileset
	if (!inherits(prefix, "plink"))
		stop("Not convinced the input is a pointer to a plink fileset.")
	if (is.null(where))
		if (!is.null(attr(prefix, "working")))
			where <- attr(prefix, "working")
	else
		where <- dirname(prefix)
	outprefix <- file.path(where, basename(prefix))
	if (!is.null(suffix))
		outprefix <- paste0(outprefix, ".", gsub("^\\.", "", suffix))
	
	if (!file.exists(paste0(prefix, ".bed")) & file.exists(paste0(prefix, ".tped"))) {
		cmd <- paste0("plink --tfile ", prefix, " --make-bed --out ", outprefix)
		system(cmd, intern = FALSE)
	}
	
	expected <- sapply(expected, function(f) paste0(outprefix, ".", f))
	for (f in expected) {
		if (file.exists(f)) {
			warning(paste0("File <", f, "> already exists; this command will overwrite it."))
			if (!overwrite)
				return(FALSE)
		}
	}
	
	if (!is.null(remove)) {
		rm.ff <- tempfile()
		write.table(remove, rm.ff, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
		flags <- paste("--remove", rm.ff, flags)
	}
	
	if (!is.null(keep)) {
		keep.ff <- tempfile()
		write.table(keep, keep.ff, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
		flags <- paste("--keep", keep.ff, flags)
	}
	
	cmd <- paste0("plink --bfile ", prefix, " ", command, " --out ", outprefix, " --allow-extra-chr ", flags)
	rez <- system(cmd, intern = capture)
	success <- sapply(expected, file.exists)
	if (capture) {
		return(rez)
	}
	else {
		if (all(success) | !length(success)) {
			return(outprefix)
		}
		else {
			warning( paste("Call to plink failed; command was '", cmd, "'") )
			warning( paste("Couldn't find output file(s):", paste(expected[!success], collapse = ",")) )
			return(FALSE)
		}
	}
	
}

#' Prune markers by pairwise LD with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param prune.by command-line flags to control LD-pruning
#' 
#' @return a pointer to a new, pruned PLINK fileset
#' 
#' @details See the relevant PLINK documentation for details of the underlying calculations.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @export
prune.plink <- function(prefix, prune.by = "--make-founders --indep 50 5 2", ...) {
	
	success <- .plink.command(prefix, prune.by, list("prune.in"))
	keep <- paste0(success, ".prune.in")
	cmd <- paste("--extract", keep, "--make-bed")
	success <- .plink.command(prefix, cmd, list("bed","bim","fam"), suffix = "pruned", ...)
	if (is.character(success)) {
		class(success) <- c("plink", class(success))
		return(success)
	}
	else {
		stop( paste("LD pruning failed.") )
	}
	
}

#' Perform a genome-wide association scan with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param model statistical model for association between genotype and phenotype (\code{"assoc"} for
#' case-control contingency table; \code{"linear"} for linear model on quantitative trait;
#' \code{"logistic"} for case-control logistic regression)
#' @param test family of hypothesis test(s) to perform
#' @param perms logical: establish significance thresholds using adaptive permutation, or skip it?
#' @param geno.missing drop markers with call rate lower than this threshold
#' @param ind.missing drop samples with call rate lower than this threshold
#' @param maf drop markers with minor-allele frequency lower than this threhsold
#' @param hwe drop markers p-value less than this threshold for test of Hardy-Weinbery equilibrium
#' @param flags additional command-line flags passed directly to PLINK call
#' 
#' @return a dataframe with association results, having the following columns:
#' \itemize{
#' \item chr
#' \item pos
#' \item marker
#' \item A1 -- reference allele (interpretation depends on dataset)
#' \item p.value -- p-value for hypothesis test
#' \item OR -- odds ratio (in case-control context)
#' \item n -- count of non-missing genotypes at this marker
#' \item test -- which family of hypothesis test this is
#' \item p.value.perm -- empirical p-value from permutations, if applicable
#' }
#' 
#' @details See the relevant PLINK documentation for details of the underlying calculations. Note
#' 	that PLINK doesn't perform kinship correction.  There now exists other, smarter software for
#' 	GWAS including \code{EMMAX}, \code{FastLMM} and \code{LIMIX}.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @export
assoc.plink <- function(prefix, model = c("assoc","linear","logistic"),
						test = c("additive","genotypic","hethom","dominant","recessive"),
						perms = TRUE, geno.missing = 0, ind.missing = 0, maf = 0, hwe = 1,
						flags = "--keep-allele-order --allow-no-sex --nonfounders", ...) {
	
	model <- match.arg(model)
	test <- match.arg(test)
	if (test == "additive") {
		test <- ""
	}
	if (model == "assoc") {
		suffix <- ".assoc"
		cmd <- paste("--assoc", flags)
	}
	else if (model == "linear") {
		suffix <- ".assoc.linear"
		cmd <- paste("--linear", test, flags)
	}
	else if (model == "logistic") {
		suffix <- ".assoc.logistic"
		cmd <- paste("--logistic", test, flags)
	}
	
	if (perms)
		cmd <- paste(cmd, "--perm")
	
	where <- dirname(prefix)
	if (geno.missing > 0)
		cmd <- paste(cmd, "--geno", geno.missing)
	if (ind.missing > 0)
		cmd <- paste(cmd, "--mind", ind.missing)
	if (maf > 0)
		cmd <- paste(cmd, "--maf", maf)
	if (hwe < 1)
		cmd <- paste(cmd, "--hwe", hwe, "midp include-nonctrl")
	
	expect <- list(gsub("^\\.", "", suffix))
	if (perms)
		expect <- c(expect, paste0(expect[[1]], ".perm"))
	success <- .plink.command(prefix, cmd, expect, ...)
	if (is.character(success)) {
		.df <- read.table(paste0(success, suffix), header = TRUE)
		if (test == "assoc") {
			df <- with(.df, data.frame(chr = CHR, pos = BP, marker = SNP, A1 = A1,
									   p.value = P, OR = OR))
		}
		else {
			df <- with(.df, data.frame(chr = CHR, pos = BP, marker = SNP, A1 = A1,
									   p.value = P, OR = OR, n = NMISS, test = TEST))
		}
		df$chr <- unplink.chroms(df$chr)
		if (perms) {
			pdf <- read.table(paste0(success, paste0(suffix, ".perm")), header = TRUE)
			rownames(pdf) <- as.character(pdf$SNP)
			df$p.value.perm <- NA
			df$p.value.perm <- pdf[ as.character(df$marker),"EMP1" ]
		}
		return(df)
	}
	else {
		stop("Association-mapping failed.")
	}
	
}

## compute genome-wide IBD estimate (\hat{pi}) with plink, optionally with initial LD pruning step
ibd.plink <- function(prefix, flags = "--nonfounders", prune = FALSE, ...) {
	
	require(Matrix)
	
	prefix.new <- prefix
	if (prune) {
		prefix.new <- prune.plink(prefix, ...)
	}
	
	cmd <- paste("--genome --min 0.0", flags)
	success <- .plink.command(prefix.new, cmd, list("genome"), ...)
	if (is.character(success)) {
		ibd.file <- paste0(success, ".genome")
		ibd <- read.table(ibd.file, header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
		iid <- union( ibd$IID1, ibd$IID2 )
		#nr <- nrow(ibd)
		#nind <- round( (1 + sqrt(1 + 4*2*nr))/2 )
		nind <- length(iid)
		K <- matrix(0, nrow = nind, ncol = nind)
		colnames(K) <- iid
		rownames(K) <- iid
		for (i in 1:nrow(ibd)) {
			K[ ibd$IID1[i], ibd$IID2[i] ] <- ibd$PI_HAT[i]
		}
		return(Matrix(K, sparse = TRUE))
	}
	else {
		stop( paste0("IBD detection failed.") )
	}
	
}

#' Compute realized kinship matrix (aka GRM) with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param flags command-line flags passed directly to underlying PLINK call
#' 
#' @return a square matrix of genetic distances (or similarities), with row and column names
#' 	equal to the sample IDs
#' 
#' @details See the relevant PLINK documentation for details of the underlying calculations.
#' 	The default is to compute the 'Mahnattan distance' between samples, ie. 1 minus the mean
#' 	number alleles shared identical-by-state across all markers.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @export
kinship.plink <- function(prefix, method = "--distance square 1-ibs", flags = "", prune = FALSE, ...) {
	
	prefix.new <- prefix
	if (prune) {
		prefix.new <- prune.plink(prefix, ...)
	}
	
	suffix <- "mdist"
	if (grepl("--make-rel", method))
		suffix <- "rel"
	
	cmd <- paste(method, flags)
	success <- .plink.command(prefix.new, cmd, list(suffix, paste0(suffix, ".id")), ...)
	if (is.character(success)) {
		kin.file <- paste0(success, ".", suffix)
		kin.ids <- paste0(success, ".", suffix, ".id")
		K <- as.matrix(read.table(kin.file))
		colnames(K) <- rownames(K) <- read.table(kin.ids, stringsAsFactors = FALSE)[,2]
		return(K)
	}
	else {
		stop("Kinship estimation failed.")
	}
	
}

#' Compute Weir & Cockerham's F_st estimator using PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param by a vector of population labels; if not specified, family ID is used
#' @param per.locus logical: if \code{TRUE}, return one value per marker; else return aggregate
#' @param chr restrict analysis to this chromosome (if \code{NULL}, all chromosomes)
#' @param flags command-line flags passed directly to underlying PLINK call
#' 
#' @return a square matrix of pairwise F_st values; diagonal has mean heterozygosity within
#' 	each population as defined by the labels in \code{by}
#' 
#' @details The F_st statistic as proposed by Sewall Wright is a measure of allele-frequency
#' 	differentiation between two populations. See the relevant PLINK documentation for details
#' 	of the underlying calculations, and Weir & Cockerham (1984) for derivation of the estimator.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 	
#' Wright S (1931) Evolution in Mendelian populations. Genetics 16: 97-159.
#' 
#' Weir BS & Cockerham CC (1984) Estimating F-statistics for the analysis of population
#' 	structure. Evolution 38(6): 1358-1370.
#' 
#' @export
weir.fst.plink <- function(prefix, by = NULL, per.locus = FALSE, chr = NULL, flags = "", ...) {
	
	if (!inherits(prefix, "plink"))
		stop("Not convinced the input is a pointer to a plink fileset.")
	
	if (!is.null(chr))
		prefix <- filter.plink(prefix, chr = chr)
	print(prefix)
	
	fam <- read.fam(prefix)
	by.flags <- "--family"
	if (!is.null(by)) {
		ff <- tempfile()
		if (length(by) != nrow(fam))
			stop("Group labels and family file have different lengths.")
		write.table(data.frame(fid = fam$fid, iid = fam$iid, group = by), ff,
					col.names = FALSE, row.names = FALSE, quote = FALSE)
		by.flags <- paste0("--within ", ff)
	}
	else {
		by <- fam$fid
	}
	
	## compute off-diagonals: Fst
	cl <- levels(factor(by))
	pairs <- combn(cl, 2)
	fst <- matrix(0, nrow = length(cl), ncol = length(cl))
	loci <- NULL
	rownames(fst) <- cl
	colnames(fst) <- cl
	for (i in seq_len(ncol(pairs))) {
		by.flags.pair <- paste(by.flags, "--keep-cluster-names", paste(as.vector(pairs[,i]), collapse = " "))
		cmd <- paste("--fst", by.flags.pair, flags)
		print(cmd)
		rez <- .plink.command(prefix, cmd, list("fst"), capture = TRUE, ...)
		if (is.character(rez)) {
			fst[ pairs[1,i], pairs[2,i] ] <- as.numeric(unlist(strsplit(rez[ length(rez) ], ":"))[2])
			fst[ pairs[2,i], pairs[1,i] ] <- fst[ pairs[1,i], pairs[2,i] ]
			if (per.locus) {
				rez <- read.table(paste0(file.path(attr(prefix, "working"), basename(prefix)), ".fst"), header = TRUE)
				colnames(rez) <- c("chr","marker","pos","n","fst")
				rez$pop1 <- factor(pairs[1,i], levels = cl)
				rez$pop2 <- factor(pairs[2,i], levels = cl)
				rez$chr <- unplink.chroms(rez$chr)
				loci <- rbind(loci, rez)
			}
			
		}
		else
			stop("Calculation of F_st failed.")
	}
	
	## compute on-diagonals: mean F
	cmd <- paste("--het")
	success <- .plink.command(prefix, cmd, list("het"), ...)
	if (is.character(success)) {
		rez <- read.table(paste0(success, ".het"), header = TRUE)
		colnames(rez) <- c("fid","iid","hom.o","hom.e","nmar","fhat")
	}
	else
		stop("Estimation of inbreeding coefficients failed.")
	
	fhat <- tapply(rez$fhat, by, mean)
	for (i in cl) {
		fst[ i,i ] <- fhat[i]
	}
	
	if (per.locus)
		return( list(matrix = fst, loci = loci) )
	else
		return(fst)
}

## get count of Hs and Ns by sample with plink
## NOT WORKING
qc.plink <- function(prefix, flags = "--nonfounders", ...) {
	
	# missingness
	missfile <- paste0(prefix, ".imiss")
	cmd <- paste("--missing --out", prefix, flags)
	success <- .plink.command(prefix, cmd, list(missfile))
	if (!success)
		stop("Attempt to compute missing-genotype rate failed.")
	
	# heterozygosity
	hetfile <- paste0(prefix, ".het")
	cmd <- paste("--het --out", prefix, flags)
	success <- .plink.command(prefix, cmd, list(hetfile), ...)
	if (!success)
		stop("Attempt to compute heterozygosity failed.")
	
	nmiss <- read.table(missfile, header = TRUE)
	nhom <- read.table(hetfile, header = TRUE)
	nmiss <- cbind(nmiss, nhom[ ,-(1:2) ])
	nmiss$H <- with(nmiss, (N_GENO - O.HOM. - N_MISS + 1))
	nmiss$N <- nmiss$N_MISS
	nmiss$valid <- nmiss$N_GENO
	rez <- nmiss[ ,c("FID","IID","H","N","valid") ]
	
	return(rez)
	
}

#' Compute pairwise LD between markers with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param index.snp name of the marker used to anchor the calculation
#' @param markers compute all pairwise LD values between markers in this list
#' @param chr limit analysis to this chromosome
#' @param from with \code{chr}, limit analysis to window starting at this position
#' @param to with \code{chr}, limit analysis to window ending at this position
#' @param window (ignored)
#' @param window.r2 keep calculating until pairwise LD falls below this threshold
#' @param flags command-line flags passed directly to underlying PLINK call
#' 
#' @return a \code{data.table} of pairwise LD values, one marker pair per row
#' 
#' @details Linkage disequlibrium (LD) is estimated here via the squared genetic correlation (r^2)
#' 	between loci.  Since the size of the output scales with the square of the number of markers and
#' 	can easily grow to gigabytes for a dataset with >10k markers, it is highly recommended
#' 	to limit the calculation using the parameters above.
#' 	
#' 	See the relevant PLINK documentation for details of the underlying calculations.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 	
#' Lewontin RC (1960) The interaction of selection and linkage. I. General considerations; heterotic
#' 	models. Genetics 49(1): 49-67.
#' 
#' @export
ld.plink <- function(prefix, index.snp = NULL, markers = NULL, chr = NULL, from = NULL, to = NULL,
					 window = NULL, window.r2 = 0, flags = "", ...) {
	
	if (!inherits(prefix, "plink"))
		stop("Not convinced the input is a pointer to a plink fileset.")
	where <- dirname(prefix)
	
	cmd <- paste0("--r2 inter-chr")
	if (!is.null(index.snp))
		cmd <- paste0(cmd, " --ld-snp ", index.snp)
	if (!is.null(markers)) {
		mkfile <- tempfile()
		write.table(data.frame(markers = markers), mkfile, quote = FALSE, col.names = FALSE, row.names = FALSE)
		cmd <- paste0(cmd, " --ld-snp-list ", mkfile)
	}
	if (!is.null(chr))
		cmd <- paste0(cmd, " --chr ", gsub("^chr","", as.character(chr)))
	if (!is.null(from))
		cmd <- paste0(cmd, " --from-bp ", formatC(from, format = "d"))
	if (!is.null(to))
		cmd <- paste0(cmd, " --to-bp ", formatC(to, format = "d"))
	if (!is.null(window.r2))
		cmd <- paste0(cmd, " --ld-window-r2 ", window.r2)
	cmd <- paste(cmd, flags)
	
	success <- .plink.command(prefix, cmd, list("ld"), ...)
	if (is.character(success)) {
		ld.file <- paste0(success, ".ld")
		ld <- data.table::data.table(read.table(ld.file, header = TRUE, strip.white = TRUE))
		return(ld)
	}
	else {
		stop( paste0("LD estimation failed.") )
	}
	
}

#' Filter markers and samples from a dataset with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param chr keep only markers on this chromosome
#' @param from with \code{chr}, keep only markers with position greater than this
#' @param to with \code{chr}, keep only markers with position less than this
#' @param maf drop markers with minor-allele frequency lower than this threhsold
#' @param hwe drop markers p-value less than this threshold for test of Hardy-Weinbery equilibrium
#' @param geno.missing drop markers with call rate lower than this threshold
#' @param ind.missing drop samples with call rate lower than this threshold
#' @param attrib with \code{attrib.file}, keep only samples with these attribute(s) (labels)
#' @param attrib.file with \code{attrib}, a file containing family IDs, sample IDs and attribute(s)
#' @param remove character vector of samples to exclude
#' @param keep character vector of samples to keep
#' @param remove.fam character vector of family IDs to exclude
#' @param keep.fam character vector of family IDs to keep
#' @param flags additional command-line flags passed directly to PLINK call
#' 
#' @return a pointer (of class \code{plink}) to a new PLINK binary fileset after application
#' 	of the above filters
#' 	
#' @details See the relevant PLINK documentation for details of the filters and precedence rules
#' 	governing how they are applied.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @export
filter.plink <- function(prefix, chr = NULL, from = NULL, to = NULL,
						 maf = 0, hwe = 0, geno.missing = 0,
						 ind.missing = 0, attrib = NULL, attrib.file = NULL,
						 remove = NULL, keep = NULL, remove.fam = NULL, keep.fam = NULL,
						 flags = "", ...) {
	
	if (!inherits(prefix, "plink"))
		stop("Not convinced the input is a pointer to a plink fileset.")
	where <- dirname(prefix)
	
	cmd <- "--make-bed"
	if (!is.null(attrib)) {
		if (file.exists(attrib.file))
			cmd <- paste(cmd, "--attrib-indiv", attrib.file, paste(attrib, collapse = ","))
		else
			stop("Filtering of individuals by attributes was requested but attribute file couldn't be found.")
	}
	if (!is.null(chr))
		cmd <- paste0(cmd, " --chr ", gsub("^chr","", as.character(chr)))
	if (!is.null(from))
		cmd <- paste0(cmd, " --from-bp ", from)
	if (!is.null(to))
		cmd <- paste0(cmd, " --to-bp ", to)
	if (ind.missing > 0)
		cmd <- paste(cmd, "--mind", ind.missing)
	if (maf > 0)
		cmd <- paste(cmd, "--maf", maf)
	if (hwe < 1)
		cmd <- paste(cmd, "--hwe", hwe, "midp include-nonctrl")
	cmd <- paste(cmd, flags)
	
	if (!is.null(remove.fam)) {
		fam <- read.fam(prefix)
		ff <- tempfile()
		write.plink.file(subset(fam, fid %in% remove.fam)[,1:2], ff)
		cmd <- paste(cmd, "--remove", ff)
	}
	
	if (!is.null(remove)) {
		fam <- read.fam(prefix)
		ff <- tempfile()
		write.plink.file(subset(fam, iid %in% remove.fam)[,1:2], ff)
		cmd <- paste(cmd, "--remove", ff)
	}
	
	if (!is.null(keep.fam)) {
		fam <- read.fam(prefix)
		ff <- tempfile()
		write.plink.file(subset(fam, fid %in% remove.fam)[,1:2], ff)
		cmd <- paste(cmd, "--keep", ff)
	}
	
	if (!is.null(keep)) {
		fam <- read.fam(prefix)
		ff <- tempfile()
		write.plink.file(subset(fam, iid %in% remove.fam)[,1:2], ff)
		cmd <- paste(cmd, "--keep", ff)
	}
	
	success <- .plink.command(prefix, cmd, list("bed","bim","fam"), suffix = "filtered", ...)
	if (is.character(success)) {
		class(success) <- c("plink", class(success))
		if (!is.null(attr(prefix, "working")))
			attr(success, "working") <- attr(prefix, "working")
	}
	else {
		stop( paste0("Filtering of plink dataset failed.") )
	}
	
	return(success)
	
}

#' Perform classical multidimensional scaling (MDS) with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param flags command-line flags passed directly to underlying PLINK call
#' @param K project samples onto this many dimensions
#' 
#' @return a dataframe with individual IDs, family IDs, and then projections in columns
#' 	"MDS1"..."MDS{k}".
#' 
#' @details See the relevant PLINK documentation for details of the underlying calculations.
#' 	The default is to perform MDS on autosomal genotypes only.  Scaled eigenvalues (ie. percent
#' 	variance explained by each dimension) are provided in \code{attr(,"eigvals")}.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 
#' @seealso \code{\link{pca.plink}}
#' 
#' @export
mds.plink <- function(prefix, flags = "--autosome", K = 3, ...) {
	
	if (!inherits(prefix, "plink"))
		stop("Not convinced the input is a pointer to a plink fileset.")
	where <- dirname(prefix)
	
	if (K < 1 || !is.numeric(K))
		stop("Can't do MDS with no dimensions; specify a number K > 0.")
	
	cmd <- paste("--cluster --mds-plot", K, "eigvals", flags)
	success <- .plink.command(prefix, cmd, list("mds","mds.eigvals"), ...)
	if (is.character(success)) {
		mds.file <- paste0(success, ".mds")
		mds <- read.table(mds.file, header = TRUE, strip.white = TRUE)
		colnames(mds) <- c("fid","iid","cluster",paste0("MDS", seq_len(K)))
		eigs <- read.table(paste0(success, ".mds.eigvals"), header = FALSE)[,1]
		attr(mds, "eigvals") <- eigs/sum(eigs)
		return(mds)
	}
	else {
		stop( paste0("MDS procedure failed.") )
	}
	
}

#' Perform PCA on genotypes with PLINK
#' 
#' @param prefix a pointer to a PLINK fileset (of class \code{plink})
#' @param flags command-line flags passed directly to underlying PLINK call
#' @param K return the projection of samples onto the top K PCs
#' @param by project individuals (\code{"indiv"}) or markers (\code{"var"}) onto PCs?
#' 
#' @return When \code{by = "indiv"} (the default), a dataframe with individual IDs, family IDs,
#' and then projections in columns "PC1"..."PC{k}".  When \code{by == "var"}, a dataframe with
#' columns chromosome and marker name followed by PCs.  The object inherits from \code{pca.result}
#' so it can be autoploted with \code{plot()}.
#' 
#' @details See the relevant PLINK documentation for details of the underlying calculations.
#' 	The default is to perform PCA on autosomal genotypes only.  Scaled eigenvalues (ie. percent
#' 	variance explained by each dimension) are provided in \code{attr(,"eigvals")}.
#' 	
#' 	A recently-proposed test for natural selection uses the squared loadings of each marker
#' 	against the top PCs as an analog of F_st.  To obtain those loadings use \code{by = "var"}.
#'
#' @references
#' PLINK v1.9: \url{https://www.cog-genomics.org/plink2}
#' 
#' Purcell S et al. (2007) PLINK: a toolset for whole-genome association and population-based 
#' 	linkage analysis. Am J Hum Genet 81(3): 559-575. doi:10.1086/519795.
#' 	
#' Duforet-Frebourg N et al. (2015) Detecting genomic signatures of natural selection with
#' 	principal componentanalysis: application to the 1000 Genomes data
#' 	\url{http://arxiv.org/abs/1504.04543}
#' 
#' @seealso \code{\link{mds.plink}}, \code{\link{plot.pca.result}}
#' 
#' @export
pca.plink <- function(prefix, flags = "--autosome", K = 20, by = c("indiv","var"), ...) {
	
	if (!inherits(prefix, "plink"))
		stop("Not convinced the input is a pointer to a plink fileset.")
	where <- dirname(prefix)
	
	if (K < 1 || !is.numeric(K))
		stop("You probably don't want a PCA result with no dimensions; specify a number K > 0.")
	
	by <- match.arg(by)
	cmd <- paste(flags, "--pca", K)
	if (by == "indiv") {
		cols <- c("fid","iid", paste0("PC", seq_len(K)))
		expect <- list("eigenvec","eigenval")
	}
	else if (by == "var") {
		cmd <- paste(cmd, "var-wts")
		cols <- c("chr","marker", paste0("PC", seq_len(K)))
		expect <- list("eigenvec.var","eigenval")
	}
	else {
		stop()
	}
	
	success <- .plink.command(prefix, cmd, expect, ...)
	if (is.character(success)) {
		pcs.file <- paste0(success, ".", expect[[1]])
		pcs <- read.table(pcs.file, header = FALSE, strip.white = TRUE)
		colnames(pcs) <- cols[1:ncol(pcs)]
		eigs <- read.table(paste0(success, ".", expect[[2]]), header = FALSE)[,1]
		attr(pcs, "eigvals") <- eigs/sum(eigs)
		attr(pcs, "explained") <- attr(pcs, "eigvals")
		class(pcs) <- c("pca.result", class(pcs))
		return(pcs)
	}
	else {
		stop( paste0("PCA procedure failed.") )
	}
	
}
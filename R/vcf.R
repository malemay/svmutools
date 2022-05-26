# These scripts extract the sequences of the variants found
# by svmu and prepare metadata for outputting in VCF format.
# Other functions take care of formatting the VCF header and
# body for output.

#' Add sequences to svmu output
#'
#' This function takes a data.frame of svmu SVs as input and returns
#' the same data.frame with the columns POS, REF and ALT appended to
#' it, thus preparing it for formatting to VCF
#'
#' @param svmu_data A data.frame of variants discovered by svmu that
#'   should have been processed by the functions \code{\link{process_insdup}}
#'   and \code{\link{filter_svmu}}. At this time, only deletions, insertions
#'   and insertions/duplications (\code{svmu_data$SV_TYPE %in% c("DEL", "INS", "INSDUP")})
#'   are supported. This means that the data.frame used for this parameter
#'   should have been filtered accordingly.
#' @param ref_fasta A character. The path to the file in fasta format
#'   containing the sequence of the reference used for the alignment.
#'   It will be automatically indexed by \code{\link[Rsamtools]{indexFa}}
#'   if its .fai index does not already exist.
#' @param query_fasta A character. Same as ref_fasta but for the query sequence
#'   used in the alignment.
#' @param mc.cores A numeric. The number of cores passed to \code{\link[parallel]{mclapply}}
#'   for processing in parallel. The default is 1 (no parallelization).
#'
#' @return A data.frame with the same number of records (rows) as the input
#'   svmu_data, but with added POS, REF and ALT columns. These columns correspond
#'   to the columns with the same name in a VCF file and will eventually be used
#'   as such in the final VCF file.
#'
#' @export
#' @examples
#' NULL
add_sequences <- function(svmu_data, ref_fasta, query_fasta, mc.cores = 1) {

	# Checking the validity of the variant types in svmu_data
	stopifnot(all(svmu_data$SV_TYPE %in% c("DEL", "INS", "INSDUP")))

	# Checking the validity of the fasta inputs and indexing if necessary
	stopifnot(file.exists(ref_fasta, query_fasta))

	if(!file.exists(paste0(ref_fasta, ".fai"))) {
		message("Indexing ", ref_fasta)
		Rsamtools::indexFa(ref_fasta)
	}

	if(!file.exists(paste0(query_fasta, ".fai"))) {
		message("Indexing ", query_fasta)
		Rsamtools::indexFa(query_fasta)
	}

	# Processing the rows in parallel on mc.cores cores
	output <- parallel::mclapply(split(svmu_data, 1:nrow(svmu_data)),
				     FUN = extract_sequences,
				     ref_fasta = ref_fasta,
				     query_fasta = query_fasta,
				     mc.cores = mc.cores)

	# Rebinding the rows together before returning
	do.call("rbind", output)
}

#' Extract the REF and ALT sequences of an SV
#'
#' This function takes a single SV (one row of a data.frame describing
#' svmu SVs) and returns the same SV with POS, REF and ALT columns
#' added to it. This function is meant to be called in parallel by
#' \code{\link{add_sequences}}.
#'
#' @param svmu_data A single row of a data.frame containing variants
#'  identified by svmu. This function processes a single row because it
#'  is meant to be called by \code{\link[base]{lapply}} or in parallel
#'  by \code{\link[parallel]{mclapply}} on a completed data.frame. Typically
#'  this function will not be called directly but instead will be called
#'  by the wrapper \code{\link{add_sequences}}.
#' @param ref_fasta A character. The path to the file containing the reference
#'   sequence in fasta format. The existence of an index is not checked
#'   by this function because this is typically done by the wrapper
#'   \code{\link{add_sequences}}. 
#' @param query_fasta A character. Same as ref_fasta but for the query
#'   sequence used in the alignment.
#'
#' @export
#' @examples
#' NULL
extract_sequences <- function(svmu_data, ref_fasta, query_fasta) {
	# Checking the validity of the svmu_data input
	stopifnot(nrow(svmu_data) == 1)

	# We assume that REF_START is never smaller than REF_END so let's check it
	stopifnot(all(svmu_data$REF_END >= svmu_data$REF_START))

	# Reading the .fai index of the reference and query
	ref_fai <- Rsamtools::scanFaIndex(ref_fasta)
	query_fai <- Rsamtools::scanFaIndex(query_fasta)
	# Extracting the maximum position on the reference and query
	ref_max <- GenomicRanges::end(ref_fai[GenomicRanges::seqnames(ref_fai) == svmu_data$REF_CHROM])
	query_max <- GenomicRanges::end(query_fai[GenomicRanges::seqnames(query_fai) == svmu_data$Q_CHROM])

	# We must first ensure that REF_END is smaller than the maximum length of the sequence
	if(svmu_data$REF_END > ref_max) {
		warning("End position of variant ", svmu_data$ID, " greater than end position of reference sequence. Setting END position to ", ref_max)
		svmu_data$REF_END <- ref_max
	}

	# We check the same thing for the query sequence (both END and START)
	if(svmu_data$Q_END > query_max) {
		warning("End position of variant ", svmu_data$ID, " greater than end position of query sequence. Setting END position to ", query_max)
		svmu_data$Q_END <- query_max
	}

	if(svmu_data$Q_START > query_max) {
		warning("Start position of variant ", svmu_data$ID, " greater than end position of query sequence. Setting START position to ", query_max)
		svmu_data$Q_START <- query_max
	}

	if(svmu_data$SV_TYPE == "DEL") {
		# --- The case if the SV is a deletion
		# For deletions, the reference range REF_START -- REF_END strictly includes
		#   the deleted nucleotides. Therefore POS needs to be at REF_START - 1
		#   if we are to include the last non-deleted nucleotide in the REF/ALT sequences
		svmu_data$POS <- svmu_data$REF_START - 1

		# This means that REF must be the sequence from POS to the last deleted nucleotide (REF_END)

		refseq <- Rsamtools::scanFa(ref_fasta,
					    param = GenomicRanges::GRanges(seqnames = svmu_data$REF_CHROM,
									   ranges = IRanges::IRanges(start = svmu_data$POS,
												     end = svmu_data$REF_END)))
		svmu_data$REF <- as.character(refseq)

		# The first nucleotide of ALT must be the same as the first of REF
		svmu_data$ALT <- substr(svmu_data$REF, 1, 1)

		# The query sequence coordinates, on the other hand, includes the aligned nucleotides
		# on either side of the deletion. The positions must therefore be excluded

		# The query start and end positions will depend on the INVERTED and Q_OVERLAP status
		# These define the beginning and end of the inserted sequence and will be used
		#   to generate the GRanges object used for extracting the sequence
		if(!svmu_data$INVERTED && !svmu_data$Q_OVERLAP) {
			qstart <- svmu_data$Q_START + 1
			qend   <- svmu_data$Q_END - 1
		} else if(!svmu_data$INVERTED && svmu_data$Q_OVERLAP) {
			qstart <- svmu_data$Q_END
			qend   <- svmu_data$Q_START - 1
		} else if(svmu_data$INVERTED && !svmu_data$Q_OVERLAP) {
			qstart <- svmu_data$Q_END + 1
			qend   <- svmu_data$Q_START - 1
		} else if(svmu_data$INVERTED && svmu_data$Q_OVERLAP) {
			qstart <- svmu_data$Q_START + 1
			qend   <- svmu_data$Q_END
		} else {
			stop("Unexpected INVERTED and Q_OVERLAP condition.")
		}
		qseq <- Rsamtools::scanFa(query_fasta,
					  param = GenomicRanges::GRanges(seqnames = svmu_data$Q_CHROM,
									 ranges = IRanges::IRanges(start = qstart,
												   end = qend)))

		if(svmu_data$INVERTED) qseq <- Biostrings::reverseComplement(qseq)

		svmu_data$ALT <- paste0(svmu_data$ALT, as.character(qseq))
	
	
	} else if (svmu_data$SV_TYPE %in% c("INS", "INSDUP")) {
		# --- The case if the SV is an insertion
		
		# In this case REF_START indicates the first nucleotide before the insertion
		# This means we can set this to POS and use this as the first nucleotide in the sequence
		svmu_data$POS <- svmu_data$REF_START

		# The end position of the deleted part is not included in the deletion and therefore we must
		# query until REF_END - 1
		# One exception to this is the special case where REF_START == REF_END ;
		#   in this case REF_END itself is the end of the query range
		refend <- if(svmu_data$REF_END > svmu_data$REF_START) svmu_data$REF_END - 1 else svmu_data$REF_END
		refseq <- Rsamtools::scanFa(ref_fasta,
					    param = GenomicRanges::GRanges(seqnames = svmu_data$REF_CHROM,
									   ranges = IRanges::IRanges(start = svmu_data$POS,
												     end = refend)))
		svmu_data$REF <- as.character(refseq)

		# The first nucleotide of ALT must be the same as the first of REF
		svmu_data$ALT <- substr(svmu_data$REF, 1, 1)
	
		# The query sequence coordinates, on the other hand, strictly include the inserted sequence
		# Coordinates on the query sequence for insertions should be strictly increasing
		# Let us verify that by making sure that !Q_END < Q_START
		if(svmu_data$Q_END < svmu_data$Q_END) {
			stop("Inserted sequence has end position smaller than start position")
		}

		# The query start and end coordinates correspond directly to Q_START and Q_END
		qseq <- Rsamtools::scanFa(query_fasta,
					  param = GenomicRanges::GRanges(seqnames = svmu_data$Q_CHROM,
									 ranges = IRanges::IRanges(start = svmu_data$Q_START,
												   end = svmu_data$Q_END)))

		# We invert the sequence if it is located in an inversion
		if(svmu_data$INVERTED) qseq <- Biostrings::reverseComplement(qseq)

		svmu_data$ALT <- paste0(svmu_data$ALT, as.character(qseq))

	} else {
		stop("SV type ", svmu_data$SV_TYPE, " not supported by svmutools")
	}

	return(svmu_data)
}

#' Format svmu data for VCF output
#'
#' This function takes a data.frame returned by \code{\link{add_sequences}}
#' as input and returns a data.frame compliant with VCF 4.3 format. This function
#' generates only the data records part of the VCF; the meta-informationheader 
#' is generated by the function \code{\link{vcf_header}}.
#'
#' @param svmu_data A data.frame with explicit sequences and position
#'   as generated by the function \code{\link{add_sequences}}, with
#'   POS, REf, and ALT columns.
#' @param sample_name A character. The name of the sample, for naming the
#'   sample in the header line.
#'
#' @return A data.frame that constitutes the body of the VCF to output,
#'   with the header line as column names and one record per row of the
#'   data.frame.
#'
#' @export
#' @examples
#' NULL
format_vcf <- function(svmu_data, sample_name) {
	# Creating the skeleton of the output VCF data
	output <- data.frame(CHROM = svmu_data$REF_CHROM,
			     POS = svmu_data$POS,
			     ID = paste0(sample_name, "_", svmu_data$ID),
			     REF = svmu_data$REF,
			     ALT = svmu_data$ALT,
			     QUAL = ".",
			     FILTER = "PASS",
			     INFO = "",
			     FORMAT = "GT",
			     stringsAsFactors = FALSE)
	
	# ----- Now formatting the INFO field
	# Starting with SVTYPE
	output$INFO <- paste0("SVTYPE=", ifelse(svmu_data$SV_TYPE == "INSDUP", "INS", svmu_data$SV_TYPE), ";")

	# Adding END (corresponds to the end position of the SV relative to the reference)
	output$INFO <- paste0(output$INFO, "END=", output$POS + nchar(output$REF) - 1, ";")

	# Adding SVLEN (difference between ALT and REF alleles)
	output$INFO <- paste0(output$INFO, "SVLEN=", nchar(output$ALT) - nchar(output$REF), ";")

	# Adding the coordinates of the alignment on the query sequence
	output$INFO <- paste0(output$INFO, "Q_CHROM=", svmu_data$Q_CHROM, ";")
	output$INFO <- paste0(output$INFO, "Q_START=", svmu_data$Q_START, ";")
	output$INFO <- paste0(output$INFO, "Q_END=", svmu_data$Q_END, ";")

	# Adding the coverage on the reference and query (COV_REF and COV_Q columns of svmu output)
	output$INFO <- paste0(output$INFO, "COVREF=", svmu_data$COV_REF, ";")
	output$INFO <- paste0(output$INFO, "COVQ=", svmu_data$COV_Q)

	# Adding optional flags indicating whether :
	# - the variant is an INSDUP
	# - the query alignments on either side of the SV overlapped
	# - the SV is located in an inverted region
	output$INFO <- paste0(output$INFO, ifelse(svmu_data$SV_TYPE == "INSDUP", ";INSDUP", ""))
	output$INFO <- paste0(output$INFO, ifelse(svmu_data$Q_OVERLAP, ";QOVERLAP", ""))
	output$INFO <- paste0(output$INFO, ifelse(svmu_data$INVERTED, ";INVERTED", ""))

	# Adding the genotype of the sample; simply the sample name with homozygous ALT genotypes at all loci
	output[[sample_name]] <- "1/1"

	return(output)
}

#' Format VCF header for svmu data
#'
#' This function generates meta-information lines for output to VCF
#' format including all the relevant information for a VCF file
#' generated from svmu using svmutools. A reference fasta file is
#' needed to generate the contig meta-information lines.
#'
#' @param ref_fasta A character. The path to the reference fasta
#'   file used for the original alignment with MUMmer. The .fai
#'   index will be generated on the fly if it does not already exists.
#'
#' @export
#' @examples
#' NULL
vcf_header <- function(ref_fasta) {
	# Creating the first few general information line
	fileformat <- "fileformat=VCFv4.3"
	filedate <- paste0("fileDate=", format(Sys.Date(), "%Y%m%d"))
	filesource <- "source=svmu+svmutools"
	reference <- paste0("reference=", basename(ref_fasta))

	# Generating the contig lines from the reference fasta
	stopifnot(file.exists(ref_fasta))

	if(!file.exists(paste0(ref_fasta, ".fai"))) {
		message("Indexing ", ref_fasta)
		Rsamtools::indexFa(ref_fasta)
	}

	# Reading the fai index
	fai_index <- Rsamtools::scanFaIndex(ref_fasta)
	
	# Creating the contig lines
	contigs <- paste0("contig=<ID=", as.character(GenomicRanges::seqnames(fai_index)),
			  ",length=", as.character(GenomicRanges::width(fai_index)),
			  ">")

	# Creating the FILTER lines
	filter_lines <- "FILTER=<ID=PASS,Description=\"All filters passed\">"

	# Creating the INFO lines
	info <- c("INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
		  "INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
		  "INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
		  "INFO=<ID=Q_CHROM,Number=1,Type=String,Description=\"Chromosome where the SV is located in the query\">",
		  "INFO=<ID=Q_START,Number=1,Type=Integer,Description=\"Start position of the SV on the query sequence\">",
		  "INFO=<ID=Q_END,Number=1,Type=Integer,Description=\"End position of the SV on the query sequence\">",
		  "INFO=<ID=COVREF,Number=1,Type=Float,Description=\"Coverage on the reference sequence as reported by svmu\">",
		  "INFO=<ID=COVQ,Number=1,Type=Float,Description=\"Coverage on the query sequence as reported by svmu\">",
		  "INFO=<ID=INSDUP,Number=0,Type=Flag,Description=\"Indicates that the insertion comprises a duplication\">",
		  "INFO=<ID=QOVERLAP,Number=0,Type=Flag,Description=\"Indicates that the query alignments overlapped at SV location\">",
		  "INFO=<ID=INVERTED,Number=0,Type=Flag,Description=\"Indicates that SV is located in an inverted alignment\">")

	# Creating the FORMAT lines
	format_lines <- "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"

	# Assmbling the lines for output
	output <- c(fileformat, filedate, filesource, reference, contigs, filter_lines, info, format_lines)

	# Adding the ## before each line and returning
	paste0("##", output)

}

#' Convert and filter SVs output by svmu to VCF format
#'
#' This function takes a file of SVs output by svmu (typically
#' called "sv.txt") and converts it to VCF format. This function
#' performs processing of combined insertion/duplications and
#' inversion flagging by default. No filters are applied by
#' default (i.e. \code{\link{filter_svmu}}) is not applied)
#' but is nevertheless strongly recommended.
#'
#' @param svmu_file A character. The path to the file containing the output
#'   of svmu (typically "sv.txt")
#' @param sample_name A character. The name of the sample. Will be used as the
#'   genotype column in the output VCF file and as a prefix for the VCF ID column.
#' @param ref_fasta A character. The path to the reference fasta file. Will be
#'   indexed by \code{\link[Rsamtools]{faiIndex}} if the corresponding .fai index
#'   does not already exist.
#' @param query_fasta A character. Same as \code{ref_fasta} but for the query sequence.
#' @param output_file A character. The VCF file to which output is directed.
#' @param mc.cores An integer. The number of cores to use for processing INSDUP
#'   records and extracting REF/ALT sequences in parallel. By default, mc.cores = 1
#'   which means that processing is not done in parallel.
#' @param logging Logical. Whether the user should be informed at various processing
#'   steps (\code{TRUE} by default).
#' @param remove_duplicates Logical. Whether duplicate records in the input file
#'   should be removed (\code{TRUE} by default). See \code{\link{read_svmu}} for more
#'   details.
#' @param insdup Logical. Whether insertions that co-occur with duplications should
#'   be processed and flagged as INSDUP (\code{TRUE} by default and recommended).
#' @param min_insdup An integer of length one. The minimum size of the insertion
#'   part of an INSDUP record for it to be kept (50 by default). This filter is
#'   only applied if \code{insdup} is \code{TRUE}.
#' @param process_inversions Logical. Whether inversions should be used to
#'   flag SVs as located in inverted regions (\code{TRUE} by default and recommended).
#' @param min_inv_distance An integer. The minimum distance (in bp) between two inverted regions
#'   for them to be merged together as a single one (see \code{\link{inverted_regions}})
#'   for more details.
#' @param inv_maxgap An integer. The maximum distance separating the location of an SV from an
#'   inverted region for this SV to be considered as located in an inverted region.
#'   See \code{\link{inverted_regions}} for more details.
#' @param apply_filters A logical. Whether to apply the filters implemented by
#'   \code{\link{filter_svmu}}. These comprise the parameters min_size, sv_types,
#'   maxcov and max_size.
#' @param min_size An integer. The minimum size of an SV for it to be kept for further
#'   processing. See \code{\link{filter_svmu}} for more details.
#' @param sv_types A character vector. The SV types to be kept for further processing.
#'   See \code{\link{filter_svmu}} for more details.
#' @param maxcov A numeric. The maximum coverage on the reference and query for the SV
#'   to be kept for further processing. See \code{\link{filter_svmu}} for more details.
#' @param max_size An integer. The maximum size for an SV to be kept for further processing.
#'   See \code{\link{filter_svmu}} for more details.
#'
#' @return
#'
#' @export
#' @examples
#' NULL
svmu_to_vcf <- function(svmu_file, sample_name, ref_fasta, query_fasta, output_file,
			mc.cores = 1, logging = TRUE,
			remove_duplicates = TRUE, insdup = TRUE, min_insdup = 50,
			process_inversions = TRUE, min_inv_distance = 5000, inv_maxgap = 10,
			apply_filters = FALSE, min_size = NULL, sv_types = NULL,
			maxcov = NULL, max_size = NULL) {

	# Some logging information
	if(logging) message("Processing sample ", sample_name)

	# Reading in the input
	svs <- read_svmu(svmu_file, remove_duplicates = remove_duplicates)
	if(logging) message(nrow(svs), " SVs in raw dataset")

	# Merging combined insertions/duplications
	if(insdup) {
		svs <- process_insdup(svs, min_size = min_insdup, mc.cores = mc.cores)
		if(logging) message(nrow(svs), " SVs after INSDUP processing")
	}

	# Reading the inverted regions and using them to flag the SVs
	if(process_inversions) {
		inversions <- inverted_regions(svmu_file, min_distance = min_inv_distance)
		if(logging) message(length(inversions), " inverted regions identified in the dataset with a total of ",
				    sum(GenomicRanges::width(inversions)) / 10^6, " Mb")
		svs <- flag_inversions(svs, inversions = inversions, maxgap = inv_maxgap)
		if(logging) message(sum(svs$INVERTED), " SVs are located in inverted regions")
	}

	if(apply_filters) {
		# Filtering out some of the SVs
		svs <- filter_svmu(svmu_data = svs,
				   min_size = min_size,
				   sv_types = sv_types,
				   maxcov = maxcov,
				   max_size = max_size)

		if(logging) message(nrow(svs), " SVs in filtered dataset")
	} else {
		if(!is.null(min_size)) warning("min_size is specified but apply_filters is FALSE")
		if(!is.null(sv_types)) warning("sv_types is specified but apply_filters is FALSE")
		if(!is.null(maxcov)) warning("maxcov is specified but apply_filters is FALSE")
		if(!is.null(max_size)) warning("max_size is specified but apply_filters is FALSE")
	}

	# Adding the explicit sequences to the data.frame (with POS, REF and ALT columns)
	if(logging) message("Extracting explicit SV sequences on ", mc.cores, " core", ifelse(mc.cores > 1, "s", ""))

	svs <- add_sequences(svmu_data = svs,
			     ref_fasta = ref_fasta, 
			     query_fasta = query_fasta, 
			     mc.cores = mc.cores)

	# Converting the records to VCF format
	if(logging) message("Converting ", nrow(svs), " records to VCF format.")

	vcf_records <- format_vcf(svs, sample_name)
	header <- vcf_header(ref_fasta)

	# Now we need to order the rows of the VCF files according to the order of the sequences in the fasta file
	fai_seqnames <- as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(ref_fasta)))
	vcf_records$CHROM <- factor(vcf_records$CHROM, levels = fai_seqnames)
	vcf_records <- vcf_records[order(vcf_records$CHROM, vcf_records$POS), ]
	vcf_records$CHROM <- as.character(vcf_records$CHROM)

	# Writing the VCF file to disk
	if(logging) message("Writing ", nrow(vcf_records), " VCF records to ", output_file)

	output_connection <- file(output_file, open = "w")
	on.exit(close(output_connection), add = TRUE)

	# Writing the meta-information lines
	writeLines(header, con = output_connection)
	# Writing the '#' in from the VCF records header line
	cat("#", file = output_connection)
	# Writing the VCF records themselves, including the header line
	write.table(vcf_records, file = output_connection, append = TRUE, quote = FALSE,
		    sep = "\t", row.names = FALSE, col.names = TRUE)

	# Return the data.frame of VCF records, invisibly
	return(invisible(vcf_records))
}


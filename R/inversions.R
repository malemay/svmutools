# This script contains function to identify inverted regions
# in the alignment between two assemblies and flag SVs located
# in inverted regions so that their sequence can be correctly
# determined when converting to VCF

#' Get GRanges of inversions from svmu output
#'
#' This function takes a file of SVs output by svmu (typically "sv.txt")
#' and returns a GRanges object of inversions found in the data. Inversion
#' are reduced such that the returned object broadly describes regions that
#' are inverted between the two aligned genomes
#'
#' @param filename a character, the name of the file of SVs output
#'   by svmu. Typically "sv.txt"
#' @param a single numeric value. The minimum distance above which
#'   two inversion ranges are not merged together.
#'
#' @return A GRanges object describing where inversions occur relative
#'   to the sequence used as a reference in the alignment
#' 
#' @export
#' @examples
#' NULL
inverted_regions <- function(filename, min_distance) {
	# Reading the data in and keeping only inversions
	svmu_data <- read_svmu(filename)
	inversions <- svmu_data[svmu_data$SV_TYPE == "INV", ]

	# Converting it to a GenomicRanges object
	inversions <- GenomicRanges::makeGRangesFromDataFrame(inversions,
							      ignore.strand = TRUE,
							      seqnames.field = "REF_CHROM",
							      start.field = "REF_START",
							      end.field = "REF_END")

	# Reducing the ranges by merging them if they are closer than min_distance
	inversions <- GenomicRanges::reduce(inversions, min.gapwidth = min_distance)
	
	inversions
}

#' Flag SVs located in inverted regions
#'
#' This function takes a GRanges object of inverted regions between
#' two alignments and a data.frame of SVs identified by svmu and adds
#' a column indicating whether this SV is located within an inverted region.
#' Columns are also added indicating whether the alignments on either
#' side of the SV overlap or not.
#'
#' @param svmu_data a data.frame returned by \code{\link{read_svmu}} or another
#'   function returning a similar data.frame, such as \code{\link{process_insdup}}
#' @param inversions a GRanges object of inversions found in the dataset, as returned
#'   by \code{\link{inverted_regions}}.
#' @param maxgap a numeric of length one. The maximum width of the gap between an SV
#'   and an inversion for the SV to be considered as overlapping the SV region.
#'   the value is 0 be default (immediately adjacent ranges are considered overlapping)
#'   but this value should probably be a small positive integer in most cases to account
#'   for the fact that e.g. some deletions may not overlap any of the ranges defined
#'   as deletions.
#'
#' @return a data.frame similar to the input svmu_data but with an additional column
#'   INVERTED which is a logical vector indicating whether the SV is located within
#'   an inversion. A column called Q_OVERLAP is also added which indicates whether the
#'   query alignments on either side of the SV overlap (TRUE) or NOT (FALSE).
#'
#' @export
#' @examples
#' NULL
flag_inversions <- function(svmu_data, inversions, maxgap = 0) {
	# First we create a GRanges version of the svmu data so we can apply
	# GRanges overlap operations on it
	svmu_granges <- GenomicRanges::makeGRangesFromDataFrame(svmu_data,
								ignore.strand = TRUE,
								seqnames.field = "REF_CHROM",
								start.field = "REF_START",
								end.field = "REF_END")

	# Checking the overlap
	svmu_data$INVERTED <- IRanges::overlapsAny(svmu_granges, inversions, maxgap = maxgap)

	# Checking whether the query alignments on either side of the SV overlap
	svmu_data$Q_OVERLAP  <- FALSE
	svmu_data[!svmu_data$INVERTED & (svmu_data$Q_END <= svmu_data$Q_START), "Q_OVERLAP"] <- TRUE
	svmu_data[svmu_data$INVERTED & (svmu_data$Q_END >= svmu_data$Q_START), "Q_OVERLAP"]  <- TRUE

	svmu_data
}


#' Generate a mummerplot of the alignments at a variant's position
#'
#' This function takes a single row from a data.frame of svmu
#' SVs as returned by \code{\link{read_svmu}} and plots the
#' alignment between the reference and query genomes at that
#' position using the mummerplot command of the MUMmer suite
#' of tools. The plot is output to a PNG file and its name
#' contains information about the variant.
#'
#' @param svmu_data a single row from a data.frame returned by
#'   \code{\link{read_svmu}} or \code{\link{process_insdup}}.
#'   The function accepts a single row because it is meant to
#'   be called on a full data.frame using \code{\link{lapply}}.
#' @param delta_file a character. The path to the delta file
#'   describing the alignment as produced by MUMmer.
#' @param expansion The expansion factor for the plotting window around
#'   the variant. By default (expansion = 1), the plotting window is 
#'   expanded by the length of the variant in both the x and y directions 
#'   along the reference and query sequences from the start and end positions
#'   of the variant. This parameter can be adjusted to enlarge or shrink
#'   the plotting window.
#' @param prefix character. A prefix that will be added to the output
#'   file name in addition to the metadata on the variant. By default
#'   the prefix is the empty string.
#' @param mp_path a character. The path to the mummerplot executable.
#'
#' @return \code{NULL}, invisibly. This function is called for its
#'   side effect of writing a plot to disk.
#'
#' @export
#' @examples
#' NULL
mummerplot <- function(svmu_data, delta_file, expansion = 1, prefix = "", mp_path) {

	# Checking that the input is just one row
	stopifnot(nrow(svmu_data) == 1)

	# We compute the start and end coordinates on the reference sequence
	if(svmu_data$REF_END >= svmu_data$REF_START) {
		xstart <- svmu_data$REF_START - svmu_data$LEN * expansion
		xend   <- svmu_data$REF_END + svmu_data$LEN * expansion
		xdir   <- "for"
	} else {
		xstart <- svmu_data$REF_END - svmu_data$LEN * expansion
		xend   <- svmu_data$REF_START + svmu_data$LEN * expansion
		xdir <- "rev"
	}

	# Then we compute the start and end coordinates on the query sequence
	if(svmu_data$Q_END >= svmu_data$Q_START) {
		ystart <- svmu_data$Q_START - svmu_data$LEN * expansion
		yend   <- svmu_data$Q_END + svmu_data$LEN * expansion
		ydir <- "for"
	} else {
		ystart <- svmu_data$Q_END - svmu_data$LEN * expansion
		yend   <- svmu_data$Q_START + svmu_data$LEN * expansion
		ydir   <- "rev"
	}

	# Making sure that the start positions are positive
	xstart <- max(1, xstart)
	ystart <- max(1, ystart)

	# Generating the name of the output png file
	output_name <- suppressWarnings(
					paste0(prefix,
					       xdir,
					       svmu_data$REF_CHROM, "_",
					       prettyNum(svmu_data$REF_START, big.mark = "."), "-",
					       prettyNum(svmu_data$REF_END, big.mark = "."), "_",
					       ydir,
					       svmu_data$Q_CHROM, "_",
					       prettyNum(svmu_data$Q_START, big.mark = "."), "-",
					       prettyNum(svmu_data$Q_END, big.mark = "."), "_",
					       svmu_data$SV_TYPE, "_", svmu_data$LEN, "_",
					       "rcov", as.character(svmu_data$COV_REF), "_",
					       "qcov", as.character(svmu_data$COV_Q),
					       if("INVERTED" %in% names(svmu_data) && svmu_data$INVERTED) "_INV" else "")
					)

	# Assembly the command that till be passed to the system function
	command <- paste0(mp_path,
			  " -terminal png",
			  " -r ", svmu_data$REF_CHROM,
			  " -q ", svmu_data$Q_CHROM,
			  " -x [", xstart, ":", xend, "]",
			  " -y [", ystart, ":", yend, "]",
			  " --prefix ", output_name, 
			  " -title ", output_name,
			  " ", delta_file)

	system(command)

	invisible(NULL)
}




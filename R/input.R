# This file contains input functions for files produced
# by svmu or MUMmer

#' Reading the SVs output by svmu
#'
#' This function simply takes the sv.txt file output by svmu and
#' returns a data.frame of that data
#'
#' @param filename a character, the name of the file (should be "sv.txt")
#' @param remove_duplicates logical. Whether duplicate entries in the
#'   SV dataset should be removed. svmu appears to return at least some
#'   CNVs in a duplicated fashion by outputting a record both for the
#'   first part of the query sequence that overlaps it and the second part
#'   of the query sequence. This option allows to keep only records that
#'   have a unique combination of REF_CHROM, REF_START, REF_END, SV_TYPE,
#'   ID and LEN. It therefore considers almost all fields but those related
#'   to the query sequence as its coordinates may be different.
#'
#' @return a data.frame describing the SVs found by svmu
#'
#' @export
#' @examples
#' NULL
read_svmu <- function(filename, remove_duplicates = TRUE) {
	svmu_data <- read.table(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

	# Removing the duplicates if remove_duplicates
	if(remove_duplicates) {
		svmu_data <- svmu_data[!duplicated(svmu_data[, c("REF_CHROM", "REF_START", "REF_END", "SV_TYPE", "ID", "LEN")]), ]
	}

	svmu_data
}


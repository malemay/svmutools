# This script contains functions to filter out SVs discovered
# by svmu based on various criteria

#' Filter out SVs discovered by svmu
#'
#' This function filters out SVs discovered by svmu according to
#' various criteria in order to keep those that are trustworth
#' and can be easily converted into sequence-explicit format.
#' What SVs to consider is left to the user, but includes only
#' DEL, INS and INSDUP by default.
#'
#' @param svmu_data a data.frame of SVs identified by svmu, as returned by
#'   functions such as \code{\link{read_svmu}} and \code{\link{process_insdup}}
#' @param min_size a numeric of length one. The minimum length of an
#'   SV to be kept (based on the LEN column).
#' @param max_size a numeric of length one. The maximum length of an
#'   SV to be kept (based on the LEN column).
#' @param sv_types a character vector. The SV types to keep. By default
#'   this includes deletions, insertions and /insertions duplications:
#'   \code{c("DEL", "INS", "INSDUP")}
#' @param maxcov a numeric of length one. Some events show increased
#'   coverage both on the reference sequence (COV_REF) and the query
#'   sequence (COV_Q) and could be difficult to convert to explicit
#'   sequence. This parameter provides a maximum value that must not
#'   be exceeded by both COV_REF and COV_Q simultaneously for the SV
#'   to be kept.
#'
#' @return a data.frame similar to the input svmu_data, but with
#'   some SVs filtered out according to the parameters.
#'
#' @export
#' @examples
#' NULL
filter_svmu <- function(svmu_data, min_size = 50, max_size = 500000,
			sv_types = c("DEL", "INS", "INSDUP"),
			maxcov = 2) {
	# Filtering for SV size
	svmu_data <- svmu_data[svmu_data$LEN >= min_size & svmu_data$LEN <= max_size, ]
	# Also removing any SV that has very large REF or Q ranges
	svmu_data <- svmu_data[!(svmu_data$REF_END - svmu_data$REF_START > max_size), ]
	svmu_data <- svmu_data[!(abs(svmu_data$Q_END - svmu_data$Q_START) > max_size), ]

	# Filtering for SV type
	svmu_data <- svmu_data[svmu_data$SV_TYPE %in% sv_types, ]

	# Filtering for coverage
	svmu_data <- svmu_data[!(svmu_data$COV_REF > maxcov & svmu_data$COV_Q > maxcov), ]

	svmu_data
}


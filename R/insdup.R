# Many insertions output by svmu are flanked by duplications
# These are represented in the sv.txt file as separate insertions
# and CNVs. However, in order to represent both breakpoints properly
# for genotyping, we need to consider them together. The functions in
# this file provide functionality for representing those situations
# together as an "INSDUP" type which eventually will be considered
# as just an INS in the final VCF

#' process_insdup
#'
#' @param svmu_data A data.frame returned by \code{\link{read_svmu}}
#'   containing a set of SVs output by svmu.
#' @param min_size a numeric. Insertions smaller than this size are flagged
#'   with INSDUP_SIZE as their SV_TYPE value.
#' @param mc.cores The number of cores to use for processing the 
#'   insertions/duplications in parallel. By default the INSDUP records
#'   are not processed in parallel (mc.cores = 1).
#'
#' @return A data.frame in which insertions that co-occur with duplications
#'   (CNVs) have their coordinates updated so that they also include the
#'   duplicated range. These insertions will be returned with the INSDUP
#'   SV_TYPE value. Records for which no unique CNV was associated with the
#'   insertion are returned with the INSDUP_FAIL SV_TYPE value. Records
#'   for which the Q_END position of the insertion was not greater or equal
#'   to the Q_START are also flagged with this SV_TYPE value.
#' @export
#' @examples
#' NULL
process_insdup <- function(svmu_data, min_size, mc.cores = 1) {

	# Splitting svmu_data into INSDUP components
	svmu_split <- split_insdup(svmu_data)

	# The "0" element represents all the records that are not INSDUP
	main_records <- svmu_split[["0"]]

	# The other elements reprsent the INSDUP records
	insdup_records <- svmu_split[names(svmu_split) != "0"]

	# For these we apply the merge_insdup function on each
	insdup_records <- parallel::mclapply(insdup_records, merge_insdup, min_size = 50, mc.cores = mc.cores)
	stopifnot(all(sapply(insdup_records, nrow) == 1))
	# And put them back into a data.frame
	insdup_records <- do.call("rbind", insdup_records)

	# We merge those results back with the "regular" records and reorder according to position
	svmu_data <- rbind(main_records, insdup_records)
	svmu_data <- svmu_data[order(svmu_data$REF_CHROM, svmu_data$REF_START), ]
	svmu_data
}

#' split_insdup
#'
#' This function accepts a data.frame of SVs output by svmu, as returned
#' by the function \code{\link{read_svmu}}. Its purpose is to identify
#' the insertion records that are associated with CNVs (insertion duplications)
#' and output a list of data.frames bundling together the insertion and its
#' associated CNV(s). The list element with the name "0" contains records that
#' are not insertion/duplication events. These should not be processed further.
#' This function is typically called at the beginning of the \code{\link{process_insdup}}
#' function.
#'
#' @param svmu_data a data.frame of SVs output by svmu, as returned
#'   by the function \code{\link{read_svmu}}.
#'
#' @return A list of data.frames bundling together the insertion and its
#'   associated CNV(s). The list element with the name "0" contains records that
#'   are not insertion/duplication events.
#' @export
#' @examples
#' NULL
split_insdup <- function(svmu_data) {
	# Identifying the insertions that correspond to the INSDUP case
	# These typically have an insertion at a single location and
	# there should be a matching CNV with its end position at that
	# location as well

	# We create an identifier that will link INS/CNV that belong
	# to a given INSDUP
	svmu_data$insdup <- 0
	index <- 1

	# Then we loop over the rows of the file to identify those pairs
	for(i in 1:nrow(svmu_data)) {
		i_row <- svmu_data[i, ]

		if(!(i_row$SV_TYPE == "INS" && i_row$REF_START == i_row$REF_END)) next

		i_cnv <- which(svmu_data$REF_CHROM == i_row$REF_CHROM & 
			       svmu_data$SV_TYPE == "CNV" &
			       (svmu_data$REF_START == i_row$REF_START | svmu_data$REF_END == i_row$REF_START))

		if(length(i_cnv)) {
			svmu_data[i, "insdup"] <- index
			svmu_data[i_cnv, "insdup"] <- index
			index <- index + 1
		}
	}

	# Splitting each such events into a data.frame each
	svmu_data$insdup <-  as.character(svmu_data$insdup)
	split(svmu_data, svmu_data$insdup)
}

#' Merge an insertion and its associated CNV record(s)
#'
#' This function processes a single INSDUP record. It is meant to be called
#' from \code{\link{process_insdup}} using lapply  on the output of \code{\link{split_insdup}}
#' to process all such records simultaneously.
#'
#' @param svmu_data A data.frame containing an INS record and at least one
#'  CNV record associated with it
#'
#' @param min_size a numeric. Insertions smaller than this size are flagged
#'   with INSDUP_SIZE as their SV_TYPE value.
#'
#' @return A single INSDUP record that represents the merged INS/CNV
#'   record. Records that could be successfully merged are returned
#'   with INSDUP as their SVTYPE value. Records that could not be successfully
#'   merged are returned with either INSDUP_FAIL or INSDUP_SIZE depending
#'   on the reason why it was not processed.
#' @export
#' @examples
#' NULL
merge_insdup <- function(svmu_data, min_size) {
	# First checking the validity of the input
	stopifnot(nrow(svmu_data) > 1)

	# Extracting the INS record from the data.frame
	ins_record <- svmu_data[svmu_data$SV_TYPE == "INS", ]
	stopifnot(nrow(ins_record) == 1)

	# If there are more than 2 records then we need to check
	# that all CNV sub-records contain similar data
	if(nrow(svmu_data) > 2) {
		cnv_record <- svmu_data[svmu_data$SV_TYPE == "CNV", ]
		cnv_record <- cnv_record[!duplicated(cnv_record[, c("REF_CHROM", "REF_START", "REF_END", "LEN")]), ]

		# If there is more than a single CNV associated with the INS, we do not keep it for downstream analyses
		if(nrow(cnv_record) > 1) {
			ins_record$SV_TYPE <- "INSDUP_FAIL"
			return(ins_record)
		}
	} else {
		cnv_record <- svmu_data[svmu_data$SV_TYPE == "CNV", ]
	}

	# We also exclude INSDUP records for which ! Q_END >= Q_START (it looks like these are mostly INS of length 1)
	if(! ins_record$Q_END >= ins_record$Q_START) {
		ins_record$SV_TYPE <- "INSDUP_FAIL"
		return(ins_record)
	}

	# We exclude insertions whose INS part is smaller than min_size
	if(ins_record$LEN < min_size) {
		ins_record$SV_TYPE <- "INSDUP_SIZE"	
		return(ins_record)
	}

	# Otherwise the new Q_END can be computed by appending the CNV range to the INS range
	# We do not add +1 because it looks like svmu ranges for insertions include only the inserted part
	ins_record$Q_END <- ins_record$Q_END + cnv_record$LEN
	ins_record$LEN   <- ins_record$LEN + cnv_record$LEN
	ins_record$SV_TYPE <- "INSDUP"

	return(ins_record)
}


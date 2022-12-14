% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insdup.R
\name{process_insdup}
\alias{process_insdup}
\title{process_insdup}
\usage{
process_insdup(svmu_data, min_size, mc.cores = 1)
}
\arguments{
\item{svmu_data}{A data.frame returned by \code{\link{read_svmu}}
containing a set of SVs output by svmu.}

\item{min_size}{a numeric. Insertions smaller than this size are flagged
with INSDUP_SIZE as their SV_TYPE value.}

\item{mc.cores}{The number of cores to use for processing the
insertions/duplications in parallel. By default the INSDUP records
are not processed in parallel (mc.cores = 1).}
}
\value{
A data.frame in which insertions that co-occur with duplications
(CNVs) have their coordinates updated so that they also include the
duplicated range. These insertions will be returned with the INSDUP
SV_TYPE value. Records for which no unique CNV was associated with the
insertion are returned with the INSDUP_FAIL SV_TYPE value. Records
for which the Q_END position of the insertion was not greater or equal
to the Q_START are also flagged with this SV_TYPE value.
}
\description{
process_insdup
}
\examples{
NULL
}

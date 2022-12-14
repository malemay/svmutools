% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter_svmu}
\alias{filter_svmu}
\title{Filter out SVs discovered by svmu}
\usage{
filter_svmu(
  svmu_data,
  min_size = 50,
  max_size = 5e+05,
  sv_types = c("DEL", "INS", "INSDUP"),
  maxcov = 2
)
}
\arguments{
\item{svmu_data}{a data.frame of SVs identified by svmu, as returned by
functions such as \code{\link{read_svmu}} and \code{\link{process_insdup}}}

\item{min_size}{a numeric of length one. The minimum length of an
SV to be kept (based on the LEN column).}

\item{max_size}{a numeric of length one. The maximum length of an
SV to be kept (based on the LEN column).}

\item{sv_types}{a character vector. The SV types to keep. By default
this includes deletions, insertions and /insertions duplications:
\code{c("DEL", "INS", "INSDUP")}}

\item{maxcov}{a numeric of length one. Some events show increased
coverage both on the reference sequence (COV_REF) and the query
sequence (COV_Q) and could be difficult to convert to explicit
sequence. This parameter provides a maximum value that must not
be exceeded by both COV_REF and COV_Q simultaneously for the SV
to be kept.}
}
\value{
a data.frame similar to the input svmu_data, but with
some SVs filtered out according to the parameters.
}
\description{
This function filters out SVs discovered by svmu according to
various criteria in order to keep those that are trustworth
and can be easily converted into sequence-explicit format.
What SVs to consider is left to the user, but includes only
DEL, INS and INSDUP by default.
}
\examples{
NULL
}

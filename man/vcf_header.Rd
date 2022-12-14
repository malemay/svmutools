% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf.R
\name{vcf_header}
\alias{vcf_header}
\title{Format VCF header for svmu data}
\usage{
vcf_header(ref_fasta)
}
\arguments{
\item{ref_fasta}{A character. The path to the reference fasta
file used for the original alignment with MUMmer. The .fai
index will be generated on the fly if it does not already exists.}
}
\description{
This function generates meta-information lines for output to VCF
format including all the relevant information for a VCF file
generated from svmu using svmutools. A reference fasta file is
needed to generate the contig meta-information lines.
}
\examples{
NULL
}

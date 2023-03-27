# Overview

The `svmutools` package includes functions to convert the output of
[svmu](https://github.com/mahulchak/svmu) (`sv.txt`) to VCF.

# Installation

## Dependencies

Some Bioconductor packages will need to be installed for this package to work:

	* Biostrings
	* GenomicRanges
	* IRanges
	* Rsamtools

Any other required packages should be pulled in automatically from CRAN.

## Installing the package

Package sources can be downloaded from GitHub by running `git clone https://github.com/malemay/svmutools`.
Then, running the following command in `R` should install the package from source:

	install.packages("svmutools", repos = NULL, type = "source")

# Documentation

There is no vignette availble at the moment for `svmutools`, however all functions
included in the package are individually documented. A complete list of the
available functions can be found here:

* `add_sequences`: Add sequences to svmu output
* `extract_sequences`: Extract the REF and ALT sequences of an SV
* `filter_svmu`: Filter out SVs discovered by svmu
* `flag_inversions`: Flag SVs located in inverted regions
* `format_vcf`: Format svmu data for VCF output
* `inverted_regions`: Get GRanges of inversions from svmu output
* `merge_insdup`: Merge an insertion and its associated CNV record(s)
* `mummerplot`: Generate a mummerplot of the alignments at a variant's position
* `process_insdup`: process_insdup
* `read_svmu`: Reading the SVs output by svmu
* `split_insdup`: split_insdup
* `svmu_to_vcf`: Convert and filter SVs output by svmu to VCF format
* `vcf_header`: Format VCF header for svmu data

# Citation

The reference to the `svmutools` package will be uploaded shortly.


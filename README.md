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

If you use this software, plase cite our publication:

Lemay, M.-A., de Ronne, M., Bélanger, R., & Belzile, F. (2023). k-mer-based GWAS enhances the discovery of causal variants and candidate genes in soybean. *The Plant Genome*, 16, e20374. [doi:10.1002/tpg2.20374](https://doi.org/10.1002/tpg2.20374)

You should also cite the original SVMU publication, as requested by the author of that software:

Chakraborty, M., Emerson, J.J., Macdonald, S.J. et al. Structural variants exhibit widespread allelic heterogeneity and shape variation in complex traits. *Nat Commun* 10, 4872 (2019). [doi:10.1038/s41467-019-12884-1](https://doi.org/10.1038/s41467-019-12884-1)

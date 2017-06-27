#!/usr/bin/env Rscript

# Command parser
doc <- "
Usage:
    extraction.R [options] <vcf_file> <sample> <output_file>

Options:
    -h, --help                      show this help message
    -p, --python                    use python 
                                        [default: false]
    -f <depth>, --filter <depth>    skip variants below <filter> depth 
                                        [default: 10]
"
opts <- docopt::docopt(doc)

# Load package
suppressPackageStartupMessages(library("CellAuthentication"))

# Extract variants
extract_variants(opts$vcf_file,
                 opts$sample,
                 opts$output_file,
                 as.numeric(opts$filter),
                 opts$python)

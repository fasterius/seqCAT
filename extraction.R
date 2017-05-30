#!/usr/bin/env Rscript

# Load package
library("CellAuthentication")

# Command parser
doc <- "
Usage:
    extraction.R [options] <vcf_file> <sample> <output_file>

Options:
    -h, --help                      show this help message
    -f <depth>, --filter <depth>    skip variants below <filter> depth
"
opts <- docopt::docopt(doc)

# Extract variants
message(paste0("reading '", basename(opts$vcf_file), "' ..."))
extract_variants(opts$vcf_file, opts$sample, opts$output_file, opts$filter)
message("finished.")

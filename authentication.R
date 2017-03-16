#!/usr/bin/env Rscript

# Load package
library("CellAuthentication")

# Command parser
doc = "
Usage: 
    authentication.R [options] <input_1> <input_2> <output>

Options:
    -h, --help                      show this help message
    -s <name>, --sample_1 <name>    name of input sample 1 [default: sample_1]
    -S <name>, --sample_2 <name>    name of input sample 2 [default: sample_2]
"
opts = docopt::docopt(doc)

# Read variant data
data_1 = read_variants(opts$input_1, opts$sample_1)
data_2 = read_variants(opts$input_2, opts$sample_2)

# Find overlaps between the variant sets
data = variant_overlaps(data_1, data_2)

# Filter variants
data = filter_variants(data)

# Write output to file
write.table(data, opts$output, sep='\t', na='', row.names=FALSE)
message('finished.')


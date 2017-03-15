#!/usr/bin/env Rscript

# Sources
source('filter_variants.R', chdir = TRUE)
source('load_packages.R', chdir = TRUE)
source('read_variants.R', chdir = TRUE)
source('variant_overlaps.R', chdir = TRUE)

# Load wrapper-specific packages
loadPackages('docopt')

# Command parser
doc = "
Usage: 
    authentication.R [options] <input_1> <input_2> <output>

Options:
    -h, --help                      show this help message
    -s <name>, --sample_1 <name>    name of input sample 1 [default: sample_1]
    -S <name>, --sample_2 <name>    name of input sample 2 [default: sample_2]
"
opts = docopt(doc)

# Read variant set 1
message(paste0('reading sample data ', basename(opts$input_1), ' ...'))
data_1 = readVariants(opts$input_1, opts$sample_1)

# Read variant set 2
message(paste0('reading sample data ', basename(opts$input_2), ' ...'))
data_2 = readVariants(opts$input_2, opts$sample_2)

# Find overlaps between the variant sets
message('finding variant overlaps ...'
data = variantOverlaps(data_1, data_2)

# Filter variants
message('filtering variants ...')
data = filterVariants(data)

# Write output to file
write.table(data, opts$output, sep='\t', na='', row.names=FALSE)
message('finished.')


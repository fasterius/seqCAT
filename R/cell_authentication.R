#!/usr/bin/env rscript

# Command parser and packages -------------------------------------------------

# Install missing packages (if applicable)
packages = c("docopt", "GenomicRanges", "plyr")
if ( length(setdiff(packages, rownames(installed.packages()))) > 0 ) {
  cat('installing missing packages ...\n')
  tryCatch (silent=TRUE,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos='http://cran.us.r-project.org'),
            warning=function(bc) {
              source("http://bioconductor.org/biocLite.R")
              biocLite(setdiff(packages, rownames(installed.packages())))
            },
            error=function(bc) {
              source("http://bioconductor.org/biocLite.R")
              biocLite(setdiff(packages, rownames(installed.packages())))
            })
}

# Command parser
suppressPackageStartupMessages(library('docopt'))
doc = "
Usage: 
    cell_authentication.R [options] <input> <output> 

Options:
    -h, --help                              show this help message
    -P, --print                             print stats of existing auth file
    -a <assembly>, --assembly <assembly>    GRCh37 or GRCh38 [default: GRCh38]
    -c <cell>, --comparison <cell>          COSMIC cell line to compare with
                                                    [default: all]
    -p <file>, --paired <file>              input a second file to compare with
    -s <name>, --sample <name>              name of input sample 
                                                    [default: sample.1]
    -S <file>, --sample.2 <name>            name of input sample 2
                                                    [default: sample.2]
"
opts = docopt(doc)

# show statistics for finished auth file --------------------------------------

if ( opts$print ) {
  
  # Read data
  data = read.table(opts$input, sep='\t', header=TRUE, stringsAsFactors=FALSE,
                    comment='', quote='\"')
  
  # Check if file is COSMIC or pairwise comparison
  if ( any(grepl('\\bmatch.cosmic\\b', names(data))) ) {  # COSMIC comparison

    # Calculate statistics
    data.cosmic = data[data$match.cosmic != '', ]
    n.cosmic = dim(data.cosmic)[1]
    n.cosmic.match = dim(data.cosmic[data.cosmic$match.cosmic == 'yes', ])[1]
    conc.cosmic = round(n.cosmic.match / n.cosmic * 100, 1)
    data.yu = data[data$match.yu != '', ]
    n.yu = dim(data.yu)[1]
    n.yu.match = dim(data.yu[data.yu$match.yu == 'yes', ])[1]
    
    # Sample and cell line names
    sample = unique(data[data$sample != '', 'sample'])
    cell.line = unique(data[data$sample.name != '', 'sample.name'])
    
    # Print statistics
    cat(paste('Sample', 'Cell line', 'Total SNPs', 'COSMIC Calls', 'Matches', 
                'Concordance', 'Yu calls', 'Yu matches', sep='\t'), '\n')
    cat(paste(sample, cell.line, dim(data)[1], n.cosmic, n.cosmic.match, 
              conc.cosmic, n.yu, n.yu.match, sep='\t'), '\n')
    
  } else {  # Pairwise comparison

    # Calculate statistics
    data.all = suppressWarnings(data[!is.na(data$match), ])
    n.all = dim(data.all)[1]
    n.all.match = dim(data.all[data.all$match=='yes', ])[1]
    conc.all = round(n.all.match / n.all * 100, 1)
    
    # Change concordance to '-' if concordance = NA
    if ( is.na(conc.all) ) {
      conc.all = 0
    }
    
    # Sample names
    sample.1 = unique(data$sample.input.1)
    sample.2 = unique(data$sample.input.2)
    
    # Print statistics
    cat(paste('Sample 1', 'Sample 2', 'Calls', 'Matches', 'Concordance',
              sep='\t'), '\n')
    cat(paste(sample.1, sample.2, n.all, n.all.match, conc.all, sep='\t'), '\n')
    
  }
  
  # Exit
  quit()
}

# Metadata function -----------------------------------------------------------

# Function for adding subject metadata columns to query metadata
addMetadata = function(query, subject, column.suffix) {
  
  # Find overlapping ranges
  hits = findOverlaps(query, subject)
  
  for ( column in names(mcols(subject)) ) {
    
    # Create empty metadata column to be filled
    mcols(query)[paste(column, column.suffix, sep='')] = NA
    
    # Convert DNAStringSet and DNAStringSetList columns to character vectors
    if (class(mcols(subject)[[column]])[1] == 'DNAStringSet') {
      mcols(subject)[column] = as.character(mcols(subject)[[column]])
    } else if (class(mcols(subject)[[column]])[1] == 'DNAStringSetList') {
      mcols(subject)[column] = 
        unstrsplit(CharacterList(mcols(subject)[[column]]))
    }
    
    # Add subject metadata to query
    mcols(query)[queryHits(hits), paste(column, column.suffix, sep='')] = 
      mcols(subject)[subjectHits(hits), column]
  }
  return(query)
}

# -----------------------------------------------------------------------------

# Read input data
cat(paste('reading sample data "', basename(opts$input), '" ...\n', sep=''))
data = read.table(opts$input, sep='\t', quote='\"', comment='', 
                  stringsAsFactors=FALSE)
colnames = c('chr', 'pos', 'rsID', 'REF', 'ALT', 'gene', 'ENSTID', 'ENSGID', 
             'impact', 'effect', 'feature.type', 'transcript.biotype', 
             'allele_1', 'allele_2','DP', 'filter', 'AD', 'nucleotide',
             'amino.acid', 'warnings')
names(data) = colnames
data$ALT = gsub('\\[|\\]', '', data$ALT)
data = data[!duplicated(data[, c('chr', 'pos', 'ENSGID')]), ]

# Add sample name
data$sample = opts$sample

# Remove indels
data = data[(nchar(data$allele_1) == 1 & nchar(data$allele_2) == 1), ]

# Convert to GRanges object
suppressPackageStartupMessages(library("GenomicRanges"))
data.gr = 
  makeGRangesFromDataFrame(data, keep.extra.columns=TRUE, 
                           ignore.strand=TRUE, seqinfo=NULL, 
                           seqnames.field='chr', start.field="pos",
                           end.field="pos", starts.in.df.are.0based=FALSE)

# Rename and remove seqlevels
seqlevels(data.gr) = gsub("chr", "", seqlevels(data.gr))
data.gr = suppressWarnings(keepSeqlevels(data.gr, 
                                         c(as.character(1:22), 'X', 'Y')))

# -----------------------------------------------------------------------------

# Optional secondary input
if ( !is.null(opts$paired) ) {
  
  cat(paste('reading sample data "', basename(opts$paired), '" ...\n', sep=''))
  
  # Read data
  data2 = read.table(opts$paired, sep='\t', quote='\"', comment='', 
                     stringsAsFactors=FALSE)
  names(data2) = colnames
  data2$ALT = gsub('\\[|\\]', '', data2$ALT)
  data2 = data2[!duplicated(data2[, c('chr', 'pos', 'ENSGID')]), ]
  
  # Remove unneeded columns
  data2 = data2[c('chr', 'pos', 'REF', 'ALT', 'allele_1', 'allele_2', 'DP',
                  'filter', 'AD', 'warnings')]
  
  # Add sample name
  data2$sample = opts$sample.2
  
  # Remove indels
  data2 = data2[(nchar(data2$allele_1) == 1 & nchar(data2$allele_2) == 1), ]
  
  data2.gr = 
    makeGRangesFromDataFrame(data2, keep.extra.columns=TRUE, 
                             ignore.strand=TRUE, seqinfo=NULL, 
                             seqnames.field='chr', start.field='pos',
                             end.field="pos", starts.in.df.are.0based=FALSE)
  
  seqlevels(data2.gr) = gsub('chr', '', seqlevels(data2.gr))
  data2.gr = suppressWarnings(keepSeqlevels(data2.gr, 
                                            c(as.character(1:22), 'X', 'Y')))
  
  # Names of files
  name.1 = strsplit(opts$sample,'/')[[1]][1]
  name.2 = strsplit(opts$sample.2,'/')[[1]][1]
  cat(paste(name.1, ' vs. ', name.2, ' ...\n', sep=''))
  
  # Merge metadata with first sample data (keeping all records)
  union.gr = union(data.gr, data2.gr)
  union.gr = addMetadata(union.gr, data.gr, '.input.1')
  union.gr = addMetadata(union.gr, data2.gr, '.input.2')
  data = as.data.frame(union.gr)
  
  # Filter variants on depth and GATK filters
  data = data[data$DP.input.1 >= 10 & data$DP.input.2 >= 10, ]
  data = data[data$filter.input.1 == 'None' & data$filter.input.2 == 'None', ]
  
  # # Remove differing REFs
  # data = data[data$REF.input.1 == data$REF.input.2, ]

  # Find variants that have genotypes in both samples
  alleles = c('allele_1.input.1', 'allele_2.input.1', 
              'allele_1.input.2', 'allele_2.input.2')
  idx.notna = row.names(subset(data, rowSums(is.na(data[, alleles])) == 0))
  
  # Check if there are overlapping variants
  if ( length(idx.notna) != 0 ) {  # At least one overlapping variant
    
    # Remove non-overlapping variants
    data = data[idx.notna, ]
    
    # Number of overlapping variants
    variants.overlap = dim(data)[1]
    
    # Check for genotype matches
    data$match = 'no'
    data.alleles = data[alleles]
    data.alleles$in.1 = paste(data.alleles[, 1], data.alleles[, 2], sep=':')
    data.alleles$in.2 = paste(data.alleles[, 3], data.alleles[, 4], sep=':')
    data.alleles$in.1.rev = paste(data.alleles[, 2], data.alleles[, 1], sep=':')
    
    idx.ok.1 = apply(data.alleles, 1, function(x) x['in.1'] %in% x['in.2'])
    idx.ok.2 = apply(data.alleles, 1, function(x) x['in.1.rev'] %in% x['in.2'])
    data[idx.ok.1, 'match'] = 'yes'
    data[idx.ok.2, 'match'] = 'yes'
    
    # Number of matching variants
    variants.match = dim(data[data$match=='yes', ])[1]
    
  } else {  # No overlapping variants
    
    # Manually set statistics to zero
    variants.overlap = 0
    variants.match = 0
    
    # Re-introduce sample names
    data[1, 'sample.input.1'] = opts$sample
    data[1, 'sample.input.2'] = opts$sample.2
    
  }
  
  # Concordance
  variants.concordance = round(variants.match/variants.overlap * 100, 1)
  
  # Change concordance to '-' if concordance = NA
  if ( is.na(variants.concordance) ) {
    variants.concordance = 0
  }

  # Print statistics
  cat(sprintf('%10s\t%8s\t%10s\n', 'Overlap', 'Matches', 'Concordance'))
  cat(sprintf('%10s\t%8s\t%10s\n', variants.overlap, variants.match, 
              variants.concordance))
  
  # Output
  write.table(data, opts$output, sep='\t', row.names=FALSE, quote=TRUE)
  quit()
}

# -----------------------------------------------------------------------------

# Check assembly version
if ( tolower(opts$assembly) == 'grch37' ) {
  assembly = '~/local/scripts/script_data/CosmicCLP_MutantExport.GRCh37.txt'
  cat('[ using GRCh37 assembly ]\n')
} else if ( tolower(opts$assembly) == 'grch38') {
  assembly = '~/local/scripts/script_data/CosmicCLP_MutantExport.GRCh38.txt'
  cat('[ using GRCh38 assembly ]\n')
} else {
  cat('ERROR: no assembly provided!\n')
  cat('Please specify either GRCh37 or GRCh38 assembly.\n')
  quit()
}

# Read COSMIC mutational data
cat('loading COSMIC mutational data ...\n')
cosmic.all = read.table(assembly, header=TRUE, sep='\t', quote='\"', 
                        comment='', stringsAsFactors=FALSE)

cosmic.all = cosmic.all[c('Gene.name', 'Accession.Number', 'Sample.name', 
                  'Mutation.ID', 'Mutation.CDS', 'Mutation.AA', 
                  'Mutation.Description', 'Mutation.zygosity', 
                  'Mutation.genome.position', 'strand', 
                  'Mutation.somatic.status', 'Mutation.verification.status')]
names(cosmic.all) = tolower(names(cosmic.all))
cosmic.all$sample.name = tolower(gsub('[-_. ]', '', cosmic.all$sample.name))

# Remove sites without listed position
cosmic.all = cosmic.all[cosmic.all$mutation.genome.position != '', ]

# Compare all or a single cell lines
if ( opts$comparison == 'all' ) {
  comparisons = unique(cosmic.all$sample.name)
} else {
  comparisons = opts$comparison
}

# Perform comparison for all selected cell lines
for ( line in comparisons ) {
  
  # Grab COSMIC data for selected cell line
  cell.line = tolower(gsub('[\\[\\]-_. ]', '', line))
  if ( cell.line %in% cosmic.all$sample.name ) {
    cosmic = cosmic.all[grep(cell.line, cosmic.all$sample.name), ]
  } else {
    cat('cell line', cell.line, 'not present in COSMIC data; aborting.\n')
    alarm()
    quit()
  }

  # Remove duplicate positions
  cosmic$gene.name = gsub('_ENST\\d+', '', cosmic$gene.name)
  cosmic = cosmic[!duplicated(cosmic[, c('gene.name', 
                                         'mutation.genome.position')]), ]
  # Only substitutions
  cosmic = cosmic[grep('Substitution', cosmic$mutation.description), ]
  n.unique.snps = dim(cosmic)[1]

  # Separate chromosome and ranges
  positions = data.frame(do.call(rbind, strsplit(as.vector(
    cosmic$mutation.genome.position), split = ":|-")))
  names(positions) = c('chr', 'start', 'end')
  cosmic = cbind(cosmic, positions)
  cosmic$start = as.numeric(as.character(cosmic$start))
  cosmic$end = as.numeric(as.character(cosmic$end))
  
  # Get reference and alternative nucleotides
  mut = as.character(cosmic$mutation.cds)
  cosmic$REF.cosmic = substr(cosmic$mutation.cds, nchar(mut)-2, nchar(mut)-2)
  cosmic$ALT.cosmic = substr(cosmic$mutation.cds, nchar(mut), nchar(mut))
  
  # Complement reverse strand mutations for comparison with all-forward VCF file
  idx = row.names(subset(cosmic, strand=='-'))
  cosmic[idx, 'REF.cosmic'] = chartr('ATCG','TAGC', cosmic[idx, 'REF.cosmic'])
  cosmic[idx, 'ALT.cosmic'] = chartr('ATCG','TAGC', cosmic[idx, 'ALT.cosmic'])
  cosmic$strand = NULL
  
  # List COSMIC alleles
  cosmic$cosmic_1 = cosmic$REF.cosmic
  cosmic[cosmic$mutation.zygosity == 'hom', 'cosmic_1'] = 
    cosmic[cosmic$mutation.zygosity == 'hom', 'ALT.cosmic']
  cosmic$cosmic_2 = cosmic$ALT.cosmic
  
  # Convert to GRanges object
  cosmic.gr = 
    makeGRangesFromDataFrame(cosmic, keep.extra.columns=TRUE, 
                             ignore.strand=TRUE, seqinfo=NULL, 
                             seqnames.field='chr', start.field="start",
                             end.field="end", starts.in.df.are.0based=FALSE)
  
  # Rename seqlevels; (23, 24) --> (X, Y)
  seqlevels(cosmic.gr, force=TRUE) = c(as.character(1:22), 'X', 'Y')
  
  # Merge COSMIC metadata to VCF GRanges object
  data = as.data.frame(addMetadata(data.gr, cosmic.gr, ''))
  
  # Filter
  data = data[data$filter=='None' & data$DP > 9, ]

  # Check for matching genotypes
  alleles = c('cosmic_1', 'cosmic_2', 'allele_1', 'allele_2')
  idx.notna = row.names(subset(data, rowSums(is.na(data[, alleles])) == 0))
  data[idx.notna, 'match.cosmic'] = 'no'
  
  data[idx.notna, 'for.1'] = paste(data[idx.notna, 'cosmic_1'], 
                                   data[idx.notna, 'cosmic_2'], sep=':')
  data[idx.notna, 'rev.1'] = paste(data[idx.notna, 'cosmic_2'], 
                                   data[idx.notna, 'cosmic_1'], sep=':')
  data['for.2'] = paste(data[, 'allele_1'], data[, 'allele_2'], sep=':')
  
  idx.ok.1 = apply(data, 1, function(x) x['for.1'] %in% x['for.2'])
  idx.ok.2 = apply(data, 1, function(x) x['rev.1'] %in% x['for.2'])
  data[idx.ok.1, 'match.cosmic'] = 'yes'
  data[idx.ok.2, 'match.cosmic'] = 'yes'
  data$for.1 = NULL
  data$rev.1 = NULL
  
  # ---------------------------------------------------------------------------

  # Concordance with Yu et al. SNPs
  yu = read.table('~/local/scripts/script_data/yu.snps.txt', header=TRUE, 
                  sep='\t')
  yu$cName = tolower(gsub('[-_. ]', '', yu$cName))
  
  # Add if cell line is available in Yu et al. dataset
  if ( cell.line %in% yu$cName ) {
    yu = t(subset(yu, cName==cell.line))
    yu = data.frame(alleles=yu[4:nrow(yu), ])
    yu$rsID = row.names(yu)
    yu = data.frame(yu[yu$alleles != '', ])
    yu$alleles = as.character(yu$alleles)
    
    # RefSNP orientation data
    strand = read.table('~/local/scripts/script_data/RefSNP.strand.txt', 
                        header=TRUE)
    yu = merge(yu, strand, by='rsID', all.x=TRUE)
    
    # Find reverse strand genotypes and make forward
    idx = row.names(subset(yu, strand=='-'))
    yu[idx, 'alleles'] = chartr('ATCG','TAGC', yu[idx, 'alleles'])
    
    # Add reverse genotype for comparison
    yu$alleles.rev = reverse(yu$alleles)
    yu$strand = NULL
    names(yu) = c('rsID', 'yu.genotype', 'rev.yu')
    
    # Merge SNP profiles with sample data
    suppressPackageStartupMessages(library(plyr))
    data = join(data, yu, by='rsID', type='left')
    
    alleles = c('for.2', 'yu.genotype', 'rev.yu')
    idx.notna = row.names(subset(data, rowSums(is.na(data[, alleles])) == 0))
    data[idx.notna, 'match.yu'] = 'no'
    
    idx.ok.1 = apply(data, 1, function(x) x['yu.genotype'] %in% x['for.2'])
    idx.ok.2 = apply(data, 1, function(x) x['rev.yu'] %in% x['for.2'])
    data[idx.ok.1, 'match.yu'] = 'yes'
    data[idx.ok.2, 'match.yu'] = 'yes'
    data$for.2 = NULL
    data$rev.yu = NULL
    
  } else {  # If cell line is not available in Yu et al. dataset
    
    data$for.2 = NULL
    data['match.yu'] = NA
    
  }

  # ---------------------------------------------------------------------------

  # Calculate and print statistics
  data.cosmic = data[!is.na(data$match.cosmic), ]
  n.cosmic = dim(data.cosmic)[1]
  n.cosmic.match = dim(data.cosmic[data.cosmic$match.cosmic=='yes', ])[1]
  conc.cosmic = round(n.cosmic.match / n.cosmic * 100, 1)
  
  data.yu = data[!is.na(data$match.yu), ]
  n.yu = dim(data.yu)[1]
  n.yu.match = dim(data.yu[data.yu$match.yu=='yes', ])[1]
  
  cat(sprintf('%11s\t%s\t%10s\t%11s\t%8s\t%10s\n', 'Unique SNPs', 'Calls', 
              'Matches', 'Concordance', 'Yu calls', 'Yu matches'))
  cat(sprintf('%11s\t%s\t%10s\t%11s\t%8s\t%10s\n', n.unique.snps, n.cosmic, 
              n.cosmic.match, conc.cosmic, n.yu, n.yu.match))
  
  # Write all variants to .txt-file for characterisation plots
  write.table(data, opts$output, sep='\t', na='', row.names=FALSE)
}


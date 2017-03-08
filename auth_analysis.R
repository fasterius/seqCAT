#!/usr/bin/Rscript

# Command parser and packages -------------------------------------------------

# Install missing packages (if applicable)
packages = c("argparse", "dplyr", "tidyr", "ggplot2")
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
suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser(epilog=gsub("[\r\n]", "", 'Analyses and visualises 
                                    cell authentication data.'))
parser$add_argument('input', type='character', help='input data file path')
parser$add_argument('output', type='character', help='output image file path')
parser$add_argument('-G', '--gene-level', action='store_true', dest='genes',
                    help='collapse variant calls to gene level')
parser$add_argument('-t', '--table', action='store_true', dest='table',
                    help='also output table with data')
parser$add_argument('-c', '--colour', type='character', default='blue',
                    help='colour palette ("blue", "grey")', dest='colour', 
                    metavar='')
parser$add_argument('-s', '--subset', type='character', default='none',
                    help='path to gene list with which to subset data',
                    metavar='', dest='subset')
args = parser$parse_args()

# Data ------------------------------------------------------------------------

# Read input data
cat(paste('reading sample data "', args$input, '" ...\n', sep=''))
data = read.table(args$input, sep='\t', quote=NULL, comment='', header=TRUE,
                  stringsAsFactors=FALSE)

# Set factor levels for matches and impacts
data$impact = data[[grep('impact', names(data), value=TRUE)[1]]]
data$impact = factor(data$impact, 
                     levels=c('HIGH', 'MODERATE', 'LOW', 'MODIFIER'))
data$match = factor(data$match, levels=c('yes', 'no'))
data$gene = data[[grep('gene', names(data), value=TRUE)[1]]]
data$ENSGID = data[[grep('ENSGID', names(data), value=TRUE)[1]]]

# Collapse calls to gene-level (if applicable)
if ( args$genes ) {
  data = data[!duplicated(data[c('match', 'impact', 'ENSGID')]), ]
}

# Subset from gene list (if applicable)
if ( args$subset != 'none' ) {
  subset = read.table(args$subset, header=TRUE, sep='\t', quote=NULL,
                      stringsAsFactors=FALSE)
  data = merge(data, subset, by='ENSGID')
}

# Calculate and print statistics
n.match = length(data[data$match == 'yes', 'match'])
n.mis = length(data[data$match == 'no', 'match'])
total = n.match + n.mis
match = round(n.match / (n.match + n.mis) * 100, 1)
cat('total\tmatches\tmismatches\tconcordance\n')
cat(paste(total, '\t', n.match, '\t', n.mis, '\t\t', match, ' %\n', sep=''))

# Impact distribution image ---------------------------------------------------
  
# Load packages
suppressPackageStartupMessages(suppressWarnings(library("dplyr")))
suppressPackageStartupMessages(suppressWarnings(library("tidyr")))
suppressPackageStartupMessages(suppressWarnings(library("ggplot2")))

# Impact groups
impact = data %>% group_by(match, impact) %>% 
  summarise(count=n()) %>% mutate(perc=round(count/sum(count) * 100, 1)) %>%
  complete(impact)
impact[is.na(impact)] = 0

gg.impact = ggplot(impact, aes(x=impact, y=perc, fill=match)) + 
  geom_bar(stat='identity', position='dodge', color='black', size=0.3) + 
  labs(x='\nImpact', y='Proportion of calls (%)\n') + 
  theme(axis.text=element_text(size=14), panel.grid.major.x=element_blank(),
        axis.title=element_text(size=17.5)) +
  geom_text(data=impact, aes(label=count), position=position_dodge(width=0.9),
            vjust=-0.5) + 
  ggtitle('SNP impact characterisation') +
  theme_bw()

# Add colours
if ( args$colour == 'blue' ) {
  gg.impact = gg.impact + 
    scale_fill_manual(values=c("#1D55A1", "#689BE3"), 
                      labels=c('Matches', 'Mismatches'), name='')
} else if ( args$colour == 'grey' ) {
  gg.impact = gg.impact + 
    scale_fill_grey(labels=c('Matches', 'Mismatches'), name='')
} else {
  cat('unknown colour choice; using greyscale.\n')
  gg.impact = gg.impact + 
    scale_fill_grey(labels=c('Matches', 'Mismatches'), name='')
}

# Save image to file
ggsave(args$output, gg.impact, dpi=300, width=7, height=7)

if ( args$table ) {
  name = gsub('.png', '.txt', args$output)
  write.table(data[data$match=='no', ], name, sep='\t', quote=TRUE, 
              row.names=FALSE)
}

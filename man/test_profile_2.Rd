% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{test_profile_2}
\alias{test_profile_2}
\title{SNV profile 2}
\format{A GRanges object with 382 elements and 17 metadata columns:
\describe{
    \item{rsID}{mutation ID, if available}
    \item{gene}{associated gene}
    \item{ENSGID}{ensembl gene ID}
    \item{ENSTID}{ensembl transcript ID}
    \item{REF}{reference allele}
    \item{ALT}{alternative allele}
    \item{impact}{putative variant impact}
    \item{effect}{variant effect}
    \item{feature}{transcript feature}
    \item{biotype}{transcript biotype}
    \item{DP}{total variant depth}
    \item{AD1}{allelic depth, allele 1}
    \item{AD2}{allelic depth, allele 2}
    \item{A1}{allele 1}
    \item{A2}{allele 2}
    \item{warnings}{warnings from variant calling}
    \item{sample}{sample name}
}}
\usage{
data(test_profile_2)
}
\description{
SNV profile in GRanges format from "sample2", originating from the
test_profile_2.txt in the inst/extdata directory, for use in unit tests.
}
\keyword{datasets}

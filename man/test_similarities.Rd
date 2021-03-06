% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{test_similarities}
\alias{test_similarities}
\title{Collated similarities object}
\format{A dataframe with 3 rows and 6 columns:
\describe{
     \item{sample_1}{name of sample 1}
     \item{sample_2}{name of sample 2}
     \item{overlaps}{the number of overlaps for the comparison}
     \item{matches}{the number of matches for the comparison}
     \item{concordance}{the concordance of the profiles}
     \item{similarity_score}{the similarity score of the profiles}
}}
\usage{
data(test_similarities)
}
\description{
Collated similarities of multiple sample comparisons from "sample1" and
"sample" from the example.vcf file, for use in unit tests.
}
\keyword{datasets}

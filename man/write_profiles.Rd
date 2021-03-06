% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_profile.R
\name{write_profiles}
\alias{write_profiles}
\title{Write SNV profiles}
\usage{
write_profiles(profile_list, format = "TXT", directory = "./")
}
\arguments{
\item{profile_list}{The SNV profiles to be written (list).}

\item{format}{The desired file format (character).}

\item{directory}{The directory to write to (path).}
}
\value{
None; writes to disk only.
}
\description{
Write several SNV profiles to file for later re-use.
}
\details{
This is a wrapper function for writing multiple SNV profiles
 present in a directory (and its sub-directories in recursive mode).
}
\examples{
# Load test profiles
data(test_profile_1)
data(test_profile_2)
profiles <- list(test_profile_1, test_profile_2)

# Write test profile to file
write_profiles(profiles, format = "TXT", directory = "./")
}

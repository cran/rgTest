% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permu_max.R
\name{permu_max}
\alias{permu_max}
\title{get the test statistic and p-value based on permutation using robust max-type edge-count test}
\usage{
permu_max(R1_list, R2_list, R1_test, R2_test, n1, n2, n_per)
}
\arguments{
\item{R1_list}{list of permuted weighted within-sample edge-counts of sample 1}

\item{R2_list}{list of permuted weighted within-sample edge-counts of sample 2}

\item{R1_test}{weighted within-sample edge-counts of sample 1}

\item{R2_test}{weighted within-sample edge-counts of sample 2}

\item{n1}{number of observations in sample 1}

\item{n2}{number of observations in sample 2}

\item{n_per}{number of permutations}
}
\value{
The p-value based on permutation distribution using robust max-type graph-based test.
}
\description{
get the test statistic and p-value based on permutation using robust max-type edge-count test
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asy_wei.R
\name{asy_wei}
\alias{asy_wei}
\title{get the approximate test statistic and p-value based on asymptotic theory using robust weighted edge-count test}
\usage{
asy_wei(asy_res, R1_test, R2_test, n1, n2)
}
\arguments{
\item{asy_res}{analytic expressions of expectations, variances and covariances}

\item{R1_test}{weighted within-sample edge-counts of sample 1}

\item{R2_test}{weighted within-sample edge-counts of sample 2}

\item{n1}{number of observations in sample 1}

\item{n2}{number of observations in sample 2}
}
\value{
A list containing the following components:
\item{test_statistic}{the asymptotic test statistic using robust weighted graph-based test.}
\item{p_value}{the asymptotic p-value using robust weighted graph-based test.}
}
\description{
get the approximate test statistic and p-value based on asymptotic theory using robust weighted edge-count test
}
\keyword{internal}

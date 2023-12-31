% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asy_gen.R
\name{asy_gen}
\alias{asy_gen}
\title{get the approximate test statistic and p-value based on asymptotic theory using robust generalized edge-count test}
\usage{
asy_gen(asy_res, R1_test, R2_test)
}
\arguments{
\item{asy_res}{analytic expressions of expectations, variances and covariances}

\item{R1_test}{weighted within-sample edge-counts of sample 1}

\item{R2_test}{weighted within-sample edge-counts of sample 2}
}
\value{
A list containing the following components:
\item{test_statistic}{the asymptotic test statistic using robust generalized graph-based test.}
\item{p_value}{the asymptotic p-value using robust generalized graph-based test.}
}
\description{
get the approximate test statistic and p-value based on asymptotic theory using robust generalized edge-count test
}
\keyword{internal}

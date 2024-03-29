% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpredBB.R
\name{dpredBB}
\alias{dpredBB}
\title{The Beta-Binomial Predictive Distribution}
\usage{
dpredBB(tpred, N, t, M, alpha = 1, beta = 1)
}
\arguments{
\item{tpred}{vector of integers: number(s) of successes for which prediction is desired}

\item{N}{number of observed independent binary variables}

\item{t}{number of observed successes out of N observations}

\item{M}{number of future observations (note tpred <= M)}

\item{alpha}{first shape factor for prior Beta distribution on theta}

\item{beta}{second shape factor for prior Beta distribution on theta}
}
\value{
Beta-Binomial predictive probability
}
\description{
dpredBB returns the predictive probability of x future successes out of M trials, given s observed successes out of a sample of size
N, and user input shape parameters for Beta prior on Pr(success)
}
\examples{
1
}

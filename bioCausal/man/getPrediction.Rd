% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFuncs.R
\name{getPrediction}
\alias{getPrediction}
\title{get prediction from multi elastic net regression}
\usage{
getPrediction(data, y, nrep = 100, seed = 123, family, alpha = 0.9, ...)
}
\arguments{
\item{data}{data (feature * samples)}

\item{y}{label}

\item{nrep}{total number of elastic net regressions, default nrep=100}

\item{seed}{random seed, default seed=123}

\item{family}{"binomial" or "gaussian" for y}

\item{alpha}{elastic net mixing parameter, default alpha=0.9}
}
\description{
get prediction from multi elastic net regression
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureSelection.R
\name{featureSelection}
\alias{featureSelection}
\title{Feature selection using chosen focal variable}
\usage{
featureSelection(
  data,
  focalVariable,
  contr = NULL,
  limmaCutoff = NULL,
  glmnetCutoff = NULL,
  class = NULL,
  nfolds = 3,
  bootstrap = T,
  alpha = 0.9,
  ...
)
}
\arguments{
\item{data}{data matrix (feature * sample)}

\item{focalVariable}{the interested focal variable used for feature selection}

\item{contr}{provide contrast for limma if focal variable is binomial}

\item{limmaCutoff}{number of top selected features when using limma for top differential features}

\item{glmnetCutoff}{number of top selected features when using multi-elastic net regression for top prediction features}

\item{class}{class of focal variables, only accepts "binomial" or "gaussian"}

\item{nfolds}{nfolds for runGnetMulti}

\item{bootstrap}{bootstrap for runGnetMulti}

\item{alpha}{alpha for runGnetMulti}
}
\value{
A list for feature selection results.
}
\description{
Feature selection using chosen focal variable
}

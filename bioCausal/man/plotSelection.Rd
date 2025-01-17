% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureSelection.R
\name{plotSelection}
\alias{plotSelection}
\title{Feature selection using chosen focal variable}
\usage{
plotSelection(labels, res, sign = F, yvalue = "pval")
}
\arguments{
\item{labels}{labels for the features}

\item{res}{the feature selection result list returned from featureSelection()}

\item{sign}{whether to plot signed values}

\item{yvalue}{pval or adj.pval for y-axis}
}
\value{
A list for feature selection results.
}
\description{
Feature selection using chosen focal variable
}

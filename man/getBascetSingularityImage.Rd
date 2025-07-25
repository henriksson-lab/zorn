% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/container.R
\name{getBascetSingularityImage}
\alias{getBascetSingularityImage}
\title{Get a Bascet image (singularity or docker).
It will be cached in the provided directory to avoid downloading it all the time}
\usage{
getBascetSingularityImage(storeAt = getwd(), tempdir = NULL)
}
\value{
A Bascet instance
}
\description{
Get a Bascet image (singularity or docker).
It will be cached in the provided directory to avoid downloading it all the time
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/debarcoding.R
\name{PrepareSharding}
\alias{PrepareSharding}
\title{Prepare to shard reads by collecting statistics about
each barcode, and filtering out cells with few reads}
\usage{
PrepareSharding(
  bascetRoot,
  inputName,
  minQuantile = 0.9,
  bascetInstance = GetDefaultBascetInstance(),
  verbose = TRUE
)
}
\arguments{
\item{minQuantile}{Read count-based cutoff for inclusion in final shards}

\item{bascetInstance}{A Bascet instance}
}
\value{
Statistics about the debarcoded reads
}
\description{
Prepare to shard reads by collecting statistics about
each barcode, and filtering out cells with few reads
}

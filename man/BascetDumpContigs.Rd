% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bascet_file.R
\name{BascetDumpContigs}
\alias{BascetDumpContigs}
\title{Store all contigs in an output directory, as cell_id.fa}
\usage{
BascetDumpContigs(
  bascetRoot,
  inputName = "skesa",
  listCells,
  outputDir,
  bascet_instance
)
}
\arguments{
\item{inputName}{Name of input shard}

\item{listCells}{List of cells to extract}

\item{outputDir}{Directory to store FASTA in}

\item{bascet_instance}{A Bascet instance}
}
\description{
Store all contigs in an output directory, as cell_id.fa
}

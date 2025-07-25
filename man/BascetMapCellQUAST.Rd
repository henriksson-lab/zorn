% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ext_tools.R
\name{BascetMapCellQUAST}
\alias{BascetMapCellQUAST}
\title{Run QUAST on reads of all cells.
This is a thin wrapper around BascetMapCell}
\usage{
BascetMapCellQUAST(
  bascetRoot,
  inputName = "filtered",
  outputName = "quast",
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
}
\arguments{
\item{bascetRoot}{The root folder where all Bascets are stored}

\item{runner}{The job manager, specifying how the command will be run (e.g. locally, or via SLURM)}

\item{bascetInstance}{Configuration for how to run the Bascet Rust API}

\item{db}{description}
}
\value{
TODO
}
\description{
Run QUAST on reads of all cells.
This is a thin wrapper around BascetMapCell
}

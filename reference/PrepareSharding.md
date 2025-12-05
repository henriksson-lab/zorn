# \############################################### \#' Take debarcoded reads and split them into suitable numbers of shards. \#' \#' The reads from one cell is guaranteed to only be present in a single shard. \#' This makes parallel processing simple as each shard can be processed on \#' a separate computer. Using more shards means that more computers can be used. \#' \#' If you perform all the calculations on a single computer, having more \#' than one shard will not result in a speedup. This option is only relevant \#' when using a cluster of compute nodes. \#' \#' \#' \#' TODO if we have multiple input samples, is there a way to group them? \#' otherwise we will be reading more input files than needed. that said, \#' if we got an index, so if list of cells specified, it is possible to quickly figure out \#' out if a file is needed at all for a merge \#' \#' TODO Figuring out if a file is needed can be done at "planning" (Zorn) stage \#' \#' TODO seems faster to have a single merger that writes multiple output files if \#' cell list is not provided. if the overhead is accepted then read all input files and \#' discard cells on the fly \#' \#' @param inputName Name of input file: Debarcoded reads \#' @param outputName Name of the output file: Properly sharded debarcoded reads \#' \#' @return A job to be executed, or being executed, depending on runner settings \#' @export BascetShardifyOld \<- function( bascetRoot, inputName="debarcoded", includeCells=NULL, \############# TODO: get rid of this parameter; only support direct merging numOutputShards=1, outputName="filtered", overwrite=FALSE, runner=GetDefaultBascetRunner(), bascetInstance=GetDefaultBascetInstance() ) input_shards \<- detectShardsForFile(bascetRoot, inputName) if(length(input_shards)==0) stop("Found no input files") \#Include all cells if nothing else provided if(is.null(includeCells)) includeCells \<- BascetCellNames(bascetRoot, inputName, bascetInstance=bascetInstance) includeCells \<- unique(includeCells\$cell) \#when shardifying, we expect cells to appear more than once – could warn for other commands! print(paste("Including all the",length(includeCells), "cells")) \#Figure out which cell goes into which shard list_cell_for_shard \<- shellscriptSplitArrayIntoListRandomly(includeCells, numOutputShards) \#Figure out input and output file names inputFiles \<- file.path(bascetRoot, input_shards) outputFiles \<- makeOutputShardNames(bascetRoot, outputName, "tirp.gz", numOutputShards) if(bascetCheckOverwriteOutput(outputFiles, overwrite)) \#Produce the script and run the job RunJob( runner = runner, jobname = "Z_shardify", bascetInstance = bascetInstance, cmd = c( \#shellscript_set_tempdir(bascetInstance), shellscriptMakeBashArray("files_out", outputFiles),

     ### Abort early if needed if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"), shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard), paste( bascetInstance@prependCmd, bascetInstance@bin, "shardify", "-t $BASCET_TEMPDIR", "-i",shellscriptMakeCommalist(inputFiles), #Need to give all input files for each job "-o ${files_out[$TASK_ID]}", #Each job produces a single output "--cells=$CELLFILE" #Each job takes its own list of cells ), "rm $CELLFILE" ), arraysize = numOutputShards ) 

else new_no_job() Prepare to shard reads by collecting statistics about
each barcode, and filtering out cells with few reads

\############################################### \#' Take debarcoded
reads and split them into suitable numbers of shards. \#' \#' The reads
from one cell is guaranteed to only be present in a single shard. \#'
This makes parallel processing simple as each shard can be processed on
\#' a separate computer. Using more shards means that more computers can
be used. \#' \#' If you perform all the calculations on a single
computer, having more \#' than one shard will not result in a speedup.
This option is only relevant \#' when using a cluster of compute nodes.
\#' \#' \#' \#' TODO if we have multiple input samples, is there a way
to group them? \#' otherwise we will be reading more input files than
needed. that said, \#' if we got an index, so if list of cells
specified, it is possible to quickly figure out \#' out if a file is
needed at all for a merge \#' \#' TODO Figuring out if a file is needed
can be done at "planning" (Zorn) stage \#' \#' TODO seems faster to have
a single merger that writes multiple output files if \#' cell list is
not provided. if the overhead is accepted then read all input files and
\#' discard cells on the fly \#' \#' @param inputName Name of input
file: Debarcoded reads \#' @param outputName Name of the output file:
Properly sharded debarcoded reads \#' \#' @return A job to be executed,
or being executed, depending on runner settings \#' @export
BascetShardifyOld \<- function( bascetRoot, inputName="debarcoded",
includeCells=NULL, \############# TODO: get rid of this parameter; only
support direct merging numOutputShards=1, outputName="filtered",
overwrite=FALSE, runner=GetDefaultBascetRunner(),
bascetInstance=GetDefaultBascetInstance() )input_shards \<-
detectShardsForFile(bascetRoot, inputName) if(length(input_shards)==0)
stop("Found no input files")#Include all cells if nothing else provided
if(is.null(includeCells)) includeCells \<- BascetCellNames(bascetRoot,
inputName, bascetInstance=bascetInstance) includeCells \<-
unique(includeCells\$cell) \#when shardifying, we expect cells to appear
more than once – could warn for other commands! print(paste("Including
all the",length(includeCells), "cells"))#Figure out which cell goes into
which shard list_cell_for_shard \<-
shellscriptSplitArrayIntoListRandomly(includeCells,
numOutputShards)#Figure out input and output file names inputFiles \<-
file.path(bascetRoot, input_shards) outputFiles \<-
makeOutputShardNames(bascetRoot, outputName, "tirp.gz",
numOutputShards)if(bascetCheckOverwriteOutput(outputFiles, overwrite))
\#Produce the script and run the job RunJob( runner = runner, jobname =
"Z_shardify", bascetInstance = bascetInstance, cmd = c(
\#shellscript_set_tempdir(bascetInstance),
shellscriptMakeBashArray("files_out", outputFiles),

        ### Abort early if needed
        if(!overwrite) shellscriptCancelJobIfFileExists("${files_out[$TASK_ID]}"),    shellscriptMakeFilesExpander("CELLFILE", list_cell_for_shard),
        paste(
          bascetInstance@prependCmd,
          bascetInstance@bin,
          "shardify",
          "-t $BASCET_TEMPDIR",
          "-i",shellscriptMakeCommalist(inputFiles), #Need to give all input files for each job
          "-o ${files_out[$TASK_ID]}",                 #Each job produces a single output
          "--cells=$CELLFILE"                          #Each job takes its own list of cells
        ),
        "rm $CELLFILE"
      ),
      arraysize = numOutputShards
    )

else new_no_job() Prepare to shard reads by collecting statistics about
each barcode, and filtering out cells with few reads

## Usage

``` r
PrepareSharding(
  bascetRoot,
  inputName = "debarcoded",
  minQuantile = 0.5,
  bascetInstance = GetDefaultBascetInstance(),
  verbose = TRUE
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- minQuantile:

  Read count-based cutoff for inclusion in final shards

- bascetInstance:

  A Bascet instance

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

Statistics about the debarcoded reads

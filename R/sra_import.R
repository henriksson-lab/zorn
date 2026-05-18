readSraList <- function(path) {
  runs <- trimws(readLines(path, warn = FALSE))
  runs <- runs[nzchar(runs) & !stringr::str_starts(runs, stringr::fixed("#"))]
  unique(runs)
}

normalizeSraRunTable <- function(runinfo) {
  stopifnot(is.data.frame(runinfo))
  stopifnot("Run" %in% names(runinfo))

  runinfo <- runinfo[!is.na(runinfo$Run) & nzchar(runinfo$Run), , drop = FALSE]
  runinfo <- runinfo[!duplicated(runinfo$Run), , drop = FALSE]
  runinfo$cell <- runinfo$Run
  runinfo <- runinfo[order(runinfo$cell, runinfo$Run), , drop = FALSE]
  rownames(runinfo) <- NULL
  runinfo
}

fetchSraRunInfoRentrez <- function(query, pageSize, sleep) {
  search <- rentrez::entrez_search(
    db = "sra",
    term = query,
    use_history = TRUE,
    retmax = 0
  )

  count <- as.integer(search$count)
  if(is.na(count) || count == 0) {
    stop(paste("No SRA runs found for query:", query))
  }

  chunks <- list()
  for(retstart in seq(0, count - 1L, by = pageSize)) {
    txt <- rentrez::entrez_fetch(
      db = "sra",
      web_history = search$web_history,
      rettype = "runinfo",
      retmode = "text",
      retstart = retstart,
      retmax = pageSize
    )
    chunks[[length(chunks) + 1L]] <- utils::read.csv(
      text = txt,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    Sys.sleep(sleep)
  }

  do.call(rbind, chunks)
}

###############################################
#' Resolve SRA inputs to a canonical run table
#'
#' Accepts direct SRR accessions, an `.sralist`, an existing RunInfo CSV/data
#' frame, or an NCBI SRA query. The returned table always contains `Run` and
#' `cell`, with `cell` defaulting to `Run` to avoid collisions. Other RunInfo
#' columns are retained as a mapping table for later renaming.
#'
#' @param runs Character vector of SRA run accessions.
#' @param sralist Path to a file with one SRA run accession per line.
#' @param runinfo SRA RunInfo data frame, or path to a RunInfo CSV file.
#' @param query NCBI SRA query passed to `rentrez`.
#' @param bioproject BioProject accession. Used to build a query when `query`
#'   is not provided.
#' @param terms Additional NCBI SRA query terms joined with `AND`.
#' @param pageSize Number of records fetched per E-utilities page.
#' @param sleep Seconds to sleep between paged E-utilities fetches.
#'
#' @return A canonical SRA run table.
#' @export
BascetResolveSra <- function(
    runs = NULL,
    sralist = NULL,
    runinfo = NULL,
    query = NULL,
    bioproject = NULL,
    terms = NULL,
    pageSize = 500,
    sleep = 1.2
) {
  sources <- c(
    !is.null(runs),
    !is.null(sralist),
    !is.null(runinfo),
    !is.null(query) || !is.null(bioproject)
  )
  if(sum(sources) != 1) {
    stop("Provide exactly one SRA source: runs, sralist, runinfo, query, or bioproject")
  }
  stopifnot(is.positive.integer(pageSize))
  stopifnot(is.numeric(sleep), length(sleep) == 1, !is.na(sleep), sleep >= 0)

  if(!is.null(runs)) {
    runinfo <- data.frame(Run = unique(as.character(runs)), stringsAsFactors = FALSE)
  } else if(!is.null(sralist)) {
    stopifnot(is.character(sralist), length(sralist) == 1, file.exists(sralist))
    runinfo <- data.frame(Run = readSraList(sralist), stringsAsFactors = FALSE)
  } else if(!is.null(runinfo)) {
    if(is.character(runinfo) && length(runinfo) == 1) {
      runinfo <- utils::read.csv(runinfo, stringsAsFactors = FALSE, check.names = FALSE)
    }
  } else {
    if(is.null(query)) {
      stopifnot(is.character(bioproject), length(bioproject) == 1, nzchar(bioproject))
      query <- paste(c(bioproject, terms), collapse = " AND ")
    }
    stopifnot(is.character(query), length(query) == 1, nzchar(query))
    runinfo <- fetchSraRunInfoRentrez(query, pageSize = pageSize, sleep = sleep)
  }

  normalizeSraRunTable(runinfo)
}

###############################################
#' Write sharded SRA accession lists for later import
#'
#' `BascetPrepareSraFetchLists()` reads an SRA RunInfo table, orders runs by
#' cell name, and writes one `.sralist` file per zorn shard. These files are
#' intended as inputs to external download jobs that later feed `bascet import-sra`.
#'
#' @param runinfo SRA RunInfo data frame, or path to a RunInfo CSV file.
#' @param bascetRoot The root folder where all Bascets are stored.
#' @param outputName Name of output shard prefix.
#' @param numShards Number of `.sralist` files to create. If NULL, this is
#'   derived from `runsPerShard`.
#' @param runsPerShard Target number of SRA runs per `.sralist` when
#'   `numShards` is NULL.
#' @param runColumn Column containing SRA run accessions.
#' @param cellColumn Column used for alphabetic cell ordering.
#' @param overwrite Overwrite existing output files.
#'
#' @return A data frame describing run-to-shard assignment.
#' @export
BascetPrepareSraFetchLists <- function(
    runinfo,
    bascetRoot,
    outputName = "tofetch",
    numShards = NULL,
    runsPerShard = 1000,
    runColumn = "Run",
    cellColumn = NULL,
    overwrite = FALSE
) {
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.null(numShards) || is.positive.integer(numShards))
  stopifnot(is.positive.integer(runsPerShard))
  stopifnot(is.character(runColumn), length(runColumn) == 1)
  stopifnot(is.null(cellColumn) || (is.character(cellColumn) && length(cellColumn) == 1))
  stopifnot(is.logical(overwrite), length(overwrite) == 1)

  if(is.character(runinfo) && length(runinfo) == 1) {
    runinfo <- utils::read.csv(runinfo, stringsAsFactors = FALSE, check.names = FALSE)
  }
  stopifnot(is.data.frame(runinfo))
  stopifnot(runColumn %in% names(runinfo))

  if(is.null(cellColumn)) {
    if("cell" %in% names(runinfo)) {
      cellColumn <- "cell"
    } else {
      cellColumn <- runColumn
    }
  }
  stopifnot(cellColumn %in% names(runinfo))

  runinfo <- runinfo[!is.na(runinfo[[runColumn]]) & nzchar(runinfo[[runColumn]]), , drop = FALSE]
  runinfo <- runinfo[order(runinfo[[cellColumn]], runinfo[[runColumn]]), , drop = FALSE]
  rownames(runinfo) <- NULL

  if(nrow(runinfo) == 0) {
    stop("RunInfo table contains no runs")
  }
  if(is.null(numShards)) {
    numShards <- ceiling(nrow(runinfo) / runsPerShard)
  }

  shard <- floor(((seq_len(nrow(runinfo)) - 1) * numShards) / nrow(runinfo)) + 1L
  outFiles <- makeOutputShardNames(bascetRoot, outputName, "sralist", numShards)
  manifestFile <- file.path(bascetRoot, paste0(outputName, ".runinfo.csv"))

  existing <- c(outFiles, manifestFile)
  if(any(file.exists(existing)) && !overwrite) {
    stop("Output files already exist; use overwrite=TRUE")
  }

  for(i in seq_len(numShards)) {
    writeLines(runinfo[[runColumn]][shard == i], con = outFiles[i])
  }

  runinfo$zorn_shard <- shard
  runinfo$zorn_sralist <- outFiles[shard]
  utils::write.csv(runinfo, file = manifestFile, row.names = FALSE)

  invisible(runinfo)
}

###############################################
#' Download sharded SRA run lists into Bascet TIRP files
#'
#' Runs `bascet import-sra` once per `.sralist` shard. The Bascet command calls
#' SRA Toolkit and writes indexed TIRP output directly.
#'
#' @param bascetRoot The root folder where all Bascets are stored.
#' @param inputName Name of input `.sralist` shard prefix.
#' @param outputName Name of output TIRP shard prefix.
#' @param runinfo Optional RunInfo CSV. Defaults to `<inputName>.runinfo.csv`
#'   if present in `bascetRoot`.
#' @param threads Threads passed to `fasterq-dump`. Defaults to runner ncpu.
#' @param runsAhead Maximum number of completed `fasterq-dump` outputs buffered
#'   ahead of TIRP writing.
#' @param sraWorkers Number of concurrent SRA Toolkit workers. If NULL, Bascet
#'   uses `threads`.
#' @param overwrite Whether to submit jobs when output files already exist.
#' @param keepTemp Keep per-run SRA/FASTQ scratch files.
#' @param prefetch Path to SRA Toolkit `prefetch` executable.
#' @param fasterqDump Path to SRA Toolkit `fasterq-dump` executable.
#' @param runner Bascet runner.
#' @param bascetInstance Bascet instance.
#'
#' @return A zorn job object.
#' @export
BascetDownloadSraRuns <- function(
    bascetRoot,
    inputName = "tofetch",
    outputName = "filtered",
    runinfo = NULL,
    threads = NULL,
    runsAhead = 10,
    sraWorkers = NULL,
    overwrite = FALSE,
    keepTemp = FALSE,
    prefetch = "prefetch",
    fasterqDump = "fasterq-dump",
    runner = GetDefaultBascetRunner(),
    bascetInstance = GetDefaultBascetInstance()
) {
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))
  stopifnot(is.null(threads) || is.positive.integer(threads))
  stopifnot(is.positive.integer(runsAhead))
  stopifnot(is.null(sraWorkers) || is.positive.integer(sraWorkers))
  stopifnot(is.logical(overwrite), length(overwrite) == 1)
  stopifnot(is.logical(keepTemp), length(keepTemp) == 1)
  stopifnot(is.character(prefetch), length(prefetch) == 1)
  stopifnot(is.character(fasterqDump), length(fasterqDump) == 1)
  if(Sys.which(prefetch) == "") {
    stop(paste0("Cannot find SRA Toolkit executable: ", prefetch))
  }
  if(Sys.which(fasterqDump) == "") {
    stop(paste0("Cannot find SRA Toolkit executable: ", fasterqDump))
  }

  sralists <- stringr::str_sort(
    list.files(
      bascetRoot,
      pattern = paste0("^", inputName, "\\.[0-9]+\\.sralist$")
    ),
    numeric = TRUE
  )
  if(length(sralists) == 0) {
    stop(paste0("No SRA list shards found for inputName=", inputName))
  }

  inputFiles <- file.path(bascetRoot, sralists)
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, "tirp.gz", length(inputFiles))

  if(is.null(runinfo)) {
    defaultRuninfo <- file.path(bascetRoot, paste0(inputName, ".runinfo.csv"))
    if(file.exists(defaultRuninfo)) {
      runinfo <- defaultRuninfo
    }
  }
  if(!is.null(runinfo)) {
    stopifnot(is.character(runinfo), length(runinfo) == 1, file.exists(runinfo))
  }

  if(is.null(threads)) {
    if("ncpu" %in% slotNames(runner)) {
      threads <- as.integer(runner@ncpu)
    } else {
      threads <- 1L
    }
  }
  if(is.na(threads) || threads < 1) {
    threads <- 1L
  }

  if(bascetCheckOverwriteOutput(outputFiles, overwrite)) {
    RunJob(
      runner = runner,
      jobname = "ZimportSRA",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        vars = list(
          files_in = inputFiles,
          files_out = outputFiles
        ),
        steps = list(
          if(!overwrite) JobSkipIfFileExists(JobVar("files_out")),
          JobBascetCommand(bascetInstance, list(
            "import-sra",
            JobArg("--sralist", JobVar("files_in"), sep = " "),
            if(!is.null(runinfo)) JobArg("--runinfo", runinfo, sep = " "),
            JobArg("--out", JobVar("files_out"), sep = " "),
            JobArg("--temp", JobEnv("BASCET_TEMPDIR"), sep = " "),
            JobArg("--threads", threads, sep = " "),
            JobArg("--runs-ahead", runsAhead, sep = " "),
            JobMaybeArg("--sra-workers", sraWorkers, sep = " "),
            JobArg("--prefetch", prefetch, sep = " "),
            JobArg("--fasterq-dump", fasterqDump, sep = " "),
            if(keepTemp) JobArg("--keep-temp")
          ))
        )
      ),
      arraysize = length(inputFiles)
    )
  } else {
    new_no_job()
  }
}

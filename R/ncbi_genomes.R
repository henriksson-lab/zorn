ncbiAssemblySummaryUrl <- function(db, group) {
  paste0("https://ftp.ncbi.nlm.nih.gov/genomes/", db, "/", group, "/assembly_summary.txt")
}

defaultNcbiGenomeCacheDir <- function() {
  file.path(defaultZornCacheDir(), "ncbi")
}

normalizeNcbiAssemblySummary <- function(path) {
  header <- readLines(path, n = 100, warn = FALSE)
  if(length(header) == 0) {
    stop("NCBI assembly summary is empty")
  }
  header_idx <- grep("^#\\s*assembly_accession\\b", header)
  if(length(header_idx) == 0) {
    header_idx <- grep("^#assembly_accession\\b", header)
  }
  column_line <- header[header_idx]
  if(length(column_line) == 0) {
    stop("Could not find NCBI assembly_summary header line")
  }
  columns <- strsplit(sub("^#\\s*", "", column_line[[1]]), "\t", fixed = FALSE)[[1]]
  columns <- make.names(columns, unique = TRUE)
  columns <- gsub("\\.", "_", columns)
  data.table::fread(
    path,
    sep = "\t",
    skip = header_idx[[1]] - 1L,
    header = TRUE,
    col.names = columns,
    data.table = FALSE,
    quote = "",
    na.strings = c("", "na", "NA")
  )
}

###############################################
#' Download and cache an NCBI genome assembly summary
#'
#' Downloads the NCBI `assembly_summary.txt` metadata file for a RefSeq/GenBank
#' group. By default this caches RefSeq bacteria metadata under
#' `file.path(tools::R_user_dir("Zorn", "cache"), "ncbi", "refseq",
#' "bacteria", "assembly_summary.txt")`. This is platform-specific and portable
#' across Linux, macOS, and Windows. If the cached file already exists, it is
#' reused unless `overwrite = TRUE`.
#'
#' @param db NCBI source database, usually `"refseq"` or `"genbank"`.
#' @param group NCBI genome group, default `"bacteria"`.
#' @param cacheDir Cache root directory. Defaults to the Zorn cache under the
#'   user cache directory.
#' @param overwrite Re-download when the cached file already exists.
#' @param dest Optional explicit output file path. When provided, this path is
#'   used instead of `file.path(cacheDir, db, group, "assembly_summary.txt")`.
#'
#' @return Path to the cached `assembly_summary.txt`.
#' @export
BascetDownloadNcbiGenomeMetadata <- function(
    db = "refseq",
    group = "bacteria",
    cacheDir = defaultNcbiGenomeCacheDir(),
    overwrite = FALSE,
    dest = NULL
) {
  stopifnot(is.character(db), length(db) == 1, nzchar(db))
  stopifnot(is.character(group), length(group) == 1, nzchar(group))
  stopifnot(is.character(cacheDir), length(cacheDir) == 1, nzchar(cacheDir))
  stopifnot(is.logical(overwrite), length(overwrite) == 1)
  stopifnot(is.null(dest) || (is.character(dest) && length(dest) == 1 && nzchar(dest)))

  if(is.null(dest)) {
    dest <- file.path(cacheDir, db, group, "assembly_summary.txt")
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if(file.exists(dest) && !overwrite) {
    return(normalizePath(dest, winslash = "/", mustWork = TRUE))
  }

  url <- ncbiAssemblySummaryUrl(db, group)
  tmp <- paste0(dest, ".tmp")
  utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
  file.rename(tmp, dest)
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

###############################################
#' Read cached NCBI genome assembly metadata
#'
#' @param path Path to `assembly_summary.txt`.
#'
#' @return A data frame with normalized column names.
#' @export
BascetReadNcbiGenomeMetadata <- function(path) {
  stopifnot(is.character(path), length(path) == 1, file.exists(path))
  normalizeNcbiAssemblySummary(path)
}

###############################################
#' Write sharded NCBI genome download inputs
#'
#' `BascetPrepareNcbiGenomeFetchLists()` filters/normalizes an NCBI assembly
#' metadata table and writes one TSV input shard per Bascet download job. The
#' generated shards contain `cell_id`, `assembly_accession`, and `ftp_path`.
#'
#' @param assemblies Assembly metadata data frame, or path to `assembly_summary.txt`.
#' @param bascetRoot The root folder where all Bascets are stored.
#' @param outputName Name of output shard prefix.
#' @param numShards Number of TSV shards to create. If NULL, derived from `genomesPerShard`.
#' @param genomesPerShard Target number of genomes per shard when `numShards` is NULL.
#' @param latest Keep only rows where `version_status == "latest"`.
#' @param excludeFromRefseq Keep only rows where `excluded_from_refseq` is empty.
#' @param assemblyLevel Optional assembly levels to keep.
#' @param refseqCategory Optional RefSeq categories to keep.
#' @param overwrite Overwrite existing output files.
#'
#' @return A data frame describing genome-to-shard assignment.
#' @export
BascetPrepareNcbiGenomeFetchLists <- function(
    assemblies,
    bascetRoot,
    outputName = "tofetch",
    numShards = NULL,
    genomesPerShard = 1000,
    latest = TRUE,
    excludeFromRefseq = TRUE,
    assemblyLevel = NULL,
    refseqCategory = NULL,
    overwrite = FALSE
) {
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.null(numShards) || is.positive.integer(numShards))
  stopifnot(is.positive.integer(genomesPerShard))
  stopifnot(is.logical(latest), length(latest) == 1)
  stopifnot(is.logical(excludeFromRefseq), length(excludeFromRefseq) == 1)
  stopifnot(is.null(assemblyLevel) || is.character(assemblyLevel))
  stopifnot(is.null(refseqCategory) || is.character(refseqCategory))
  stopifnot(is.logical(overwrite), length(overwrite) == 1)

  if(is.character(assemblies) && length(assemblies) == 1) {
    assemblies <- BascetReadNcbiGenomeMetadata(assemblies)
  }
  stopifnot(is.data.frame(assemblies))
  required <- c("assembly_accession", "ftp_path")
  missing <- setdiff(required, names(assemblies))
  if(length(missing) > 0) {
    stop(paste("Missing required assembly metadata columns:", paste(missing, collapse = ", ")))
  }

  keep <- !is.na(assemblies$assembly_accession) &
    nzchar(assemblies$assembly_accession) &
    !is.na(assemblies$ftp_path) &
    nzchar(assemblies$ftp_path)
  if(latest && "version_status" %in% names(assemblies)) {
    keep <- keep & assemblies$version_status == "latest"
  }
  if(excludeFromRefseq && "excluded_from_refseq" %in% names(assemblies)) {
    keep <- keep & (is.na(assemblies$excluded_from_refseq) | assemblies$excluded_from_refseq == "")
  }
  if(!is.null(assemblyLevel)) {
    stopifnot("assembly_level" %in% names(assemblies))
    keep <- keep & assemblies$assembly_level %in% assemblyLevel
  }
  if(!is.null(refseqCategory)) {
    stopifnot("refseq_category" %in% names(assemblies))
    keep <- keep & assemblies$refseq_category %in% refseqCategory
  }

  assemblies <- assemblies[keep, , drop = FALSE]
  assemblies <- assemblies[order(assemblies$assembly_accession), , drop = FALSE]
  rownames(assemblies) <- NULL
  if(nrow(assemblies) == 0) {
    stop("No NCBI genome assemblies remain after filtering")
  }
  if(is.null(numShards)) {
    numShards <- ceiling(nrow(assemblies) / genomesPerShard)
  }

  shard <- floor(((seq_len(nrow(assemblies)) - 1) * numShards) / nrow(assemblies)) + 1L
  outFiles <- makeOutputShardNames(bascetRoot, outputName, "ncbilist", numShards)
  manifestFile <- file.path(bascetRoot, paste0(outputName, ".assemblies.csv"))
  existing <- c(outFiles, manifestFile)
  if(!overwrite) {
    if(all(file.exists(existing))) {
      print("All NCBI genome fetch-list files already exist; reusing them. To change this behaviour, set overwrite=TRUE")
      return(invisible(utils::read.csv(manifestFile, stringsAsFactors = FALSE, check.names = FALSE)))
    } else if(any(file.exists(existing))) {
      print("Some NCBI genome fetch-list files exist, but not all(!). Fetch lists will be regenerated")
    }
  }

  for(i in seq_len(numShards)) {
    shard_table <- data.frame(
      cell_id = assemblies$assembly_accession[shard == i],
      assembly_accession = assemblies$assembly_accession[shard == i],
      ftp_path = assemblies$ftp_path[shard == i],
      stringsAsFactors = FALSE
    )
    utils::write.table(
      shard_table,
      file = outFiles[i],
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  assemblies$zorn_shard <- shard
  assemblies$zorn_fetch_tsv <- outFiles[shard]
  utils::write.csv(assemblies, file = manifestFile, row.names = FALSE)
  invisible(assemblies)
}

###############################################
#' Download sharded NCBI genome inputs into Bascet files
#'
#' Runs `bascet ncbi-genome-download` once per `.ncbilist` shard produced by
#' `BascetPrepareNcbiGenomeFetchLists()`.
#'
#' @param bascetRoot The root folder where all Bascets are stored.
#' @param inputName Name of input `.ncbilist` shard prefix.
#' @param outputName Name of output shard prefix.
#' @param outFormat Output shard extension. Use `"zip"` to store
#'   `{assembly}/contigs.fa` entries, or `"tirp.gz"` to store contigs as R1 reads
#'   with dummy `F` quality scores and empty R2/UMI.
#' @param threads Parallel genome workers per Bascet job. Defaults to runner ncpu.
#' @param downloadStartsPerSecond Global per-job NCBI download start rate.
#' @param queueSize Optional completed-fragment queue size.
#' @param maxRetries Download retry count.
#' @param keepTemp Keep per-job temporary download and fragment files.
#' @param overwrite Whether to submit jobs when output files already exist.
#' @param runner Bascet runner.
#' @param bascetInstance Bascet instance.
#'
#' @return A zorn job object.
#' @export
BascetDownloadNcbiGenomes <- function(
    bascetRoot,
    inputName = "tofetch",
    outputName = "contigs",
    outFormat = "zip",
    threads = NULL,
    downloadStartsPerSecond = 2,
    queueSize = NULL,
    maxRetries = 5,
    keepTemp = FALSE,
    overwrite = FALSE,
    runner = GetDefaultBascetRunner(),
    bascetInstance = GetDefaultBascetInstance()
) {
  stopifnot(dir.exists(bascetRoot))
  bascetRoot <- normalizeBascetRoot(bascetRoot)
  stopifnot(is.valid.shardname(inputName))
  stopifnot(is.valid.shardname(outputName))
  stopifnot(is.character(outFormat), length(outFormat) == 1)
  stopifnot(outFormat %in% c("zip", "tirp.gz", "tirp"))
  stopifnot(is.null(threads) || is.valid.threadcount(threads))
  stopifnot(is.numeric(downloadStartsPerSecond), length(downloadStartsPerSecond) == 1,
            !is.na(downloadStartsPerSecond), is.finite(downloadStartsPerSecond),
            downloadStartsPerSecond > 0)
  stopifnot(is.null(queueSize) || is.positive.integer(queueSize))
  stopifnot(is.positive.integer(maxRetries))
  stopifnot(is.logical(keepTemp), length(keepTemp) == 1)
  stopifnot(is.logical(overwrite), length(overwrite) == 1)
  stopifnot(is.runner(runner))
  stopifnot(is.bascet.instance(bascetInstance))

  inputShards <- stringr::str_sort(
    list.files(
      bascetRoot,
      pattern = paste0("^", inputName, "\\.[0-9]+\\.ncbilist$")
    ),
    numeric = TRUE
  )
  if(length(inputShards) == 0) {
    stop(paste0("No NCBI genome list shards found for inputName=", inputName))
  }
  inputFiles <- file.path(bascetRoot, inputShards)
  outputFiles <- makeOutputShardNames(bascetRoot, outputName, outFormat, length(inputFiles))

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
      jobname = "ZncbiGenomes",
      bascetInstance = bascetInstance,
      cmd = JobScript(
        vars = list(
          files_in = inputFiles,
          files_out = outputFiles
        ),
        steps = list(
          if(!overwrite) JobSkipIfFileExists(JobVar("files_out")),
          JobBascetCommand(bascetInstance, list(
            "ncbi-genome-download",
            JobArg("--input", JobVar("files_in"), sep = " "),
            JobArg("--out", JobVar("files_out"), sep = " "),
            JobArg("--temp", JobEnv("BASCET_TEMPDIR"), sep = " "),
            JobArg("--threads", threads, sep = " "),
            JobMaybeArg("--queue-size", queueSize, sep = " "),
            JobArg("--download-starts-per-second", downloadStartsPerSecond, sep = " "),
            JobArg("--max-retries", maxRetries, sep = " "),
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

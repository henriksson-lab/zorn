################################################################################
################ Structured job scripts ########################################
################################################################################

###############################################
#' Check that parameter is a valid environment variable name
#'
#' @param x A string representing an environment variable name.
#'
#' @return TRUE if `x` is a scalar shell-compatible environment variable name,
#'   otherwise FALSE.
#' @noRd
is.valid.env.variable <- function(x) {
  is.character(x) && length(x) == 1 && stringr::str_detect(x, "^[A-Za-z_][A-Za-z0-9_]*$")
}

###############################################
#' Check that an object is a JobScript
#'
#' @param x Object to test.
#'
#' @return TRUE if `x` inherits from `JobScript`, otherwise FALSE.
#' @noRd
is.jobscript <- function(x) {
  inherits(x, "JobScript")
}

###############################################
#' Build a structured job script
#'
#' Runners consume this through their target backend. SLURM renders to bash;
#' local execution interprets it directly.
#'
#' @param vars Named list of per-task values
#' @param steps List of job steps
#'
#' @return A JobScript object
#' @noRd
JobScript <- function(vars = list(), steps = list()) {
  stopifnot(is.list(vars), is.list(steps))
  if(!is.null(names(vars))) {
    for(nm in names(vars)) {
      stopifnot(is.valid.env.variable(nm))
    }
  }

  structure(
    list(
      vars = vars,
      steps = Filter(Negate(is.null), steps)
    ),
    class = "JobScript"
  )
}

###############################################
#' Reference a per-task JobScript variable
#'
#' `JobVar` values are resolved as `script$vars[[name]][task_id + 1]` by the
#' local runner and as `${name[$TASK_ID]}` by the bash renderer.
#'
#' @param name Variable name.
#'
#' @return A `JobVar` placeholder.
#' @noRd
JobVar <- function(name) {
  stopifnot(is.valid.env.variable(name))
  structure(list(name = name), class = "JobVar")
}

###############################################
#' Reference an environment variable
#'
#' `JobEnv` values are resolved from `Sys.getenv()` locally and rendered as
#' `$NAME` for bash.
#'
#' @param name Environment variable name.
#'
#' @return A `JobEnv` placeholder.
#' @noRd
JobEnv <- function(name) {
  stopifnot(is.valid.env.variable(name))
  structure(list(name = name), class = "JobEnv")
}

###############################################
#' Read a file and join lines with commas
#'
#' This is used for command arguments that expect a comma-separated list loaded
#' from a file generated earlier in the same job.
#'
#' @param path File path value. May be a literal path or another job value such
#'   as `JobVar`.
#'
#' @return A `JobFileLinesCsv` placeholder.
#' @noRd
JobFileLinesCsv <- function(path) {
  structure(list(path = path), class = "JobFileLinesCsv")
}

###############################################
#' Build a command argument
#'
#' `sep = "="` renders/resolves as one argv element like `--flag=value`.
#' `sep = " "` renders/resolves as two argv elements like `--flag value`.
#'
#' @param flag Flag name, or NULL for positional arguments.
#' @param value Argument value. May be a literal value or a job placeholder.
#' @param sep Separator between flag and value.
#'
#' @return A `JobArg` object.
#' @noRd
JobArg <- function(flag = NULL, value = NULL, sep = c("=", " ")) {
  sep <- match.arg(sep)
  structure(
    list(flag = flag, value = value, sep = sep),
    class = "JobArg"
  )
}

###############################################
#' Build an optional command argument
#'
#' Returns NULL when `value` is NULL so callers can include optional arguments
#' directly inside `list(...)` and let `JobCommand`/`JobBascetCommand` filter
#' them out.
#'
#' @param flag Flag name, or NULL for positional arguments.
#' @param value Argument value. NULL suppresses the argument.
#' @param format Function applied to non-NULL values before constructing the
#'   argument.
#' @param sep Separator between flag and value.
#'
#' @return A `JobArg` object, or NULL.
#' @noRd
JobMaybeArg <- function(flag = NULL, value = NULL, format = identity, sep = c("=", " ")) {
  if(is.null(value)) {
    return(NULL)
  }
  JobArg(flag, format(value), sep = match.arg(sep))
}

###############################################
#' Skip current task if a file exists
#'
#' Local execution returns a skip sentinel. Bash rendering emits an early
#' `exit 0` block for the task.
#'
#' @param path File path value. May be a literal path or a job placeholder.
#'
#' @return A `JobSkipIfFileExists` step.
#' @noRd
JobSkipIfFileExists <- function(path) {
  structure(list(path = path), class = "JobSkipIfFileExists")
}

###############################################
#' Make temporary per-task files from R-side list content
#'
#' Each task gets one generated file. Local execution creates the current task's
#' file lazily; bash rendering creates all files before script submission and
#' exposes them as a per-task array.
#'
#' @param name Job variable that will hold per-task file paths.
#' @param list_content List of file contents, one element per task.
#'
#' @return A `JobFiles` step.
#' @noRd
JobFiles <- function(name, list_content) {
  stopifnot(is.valid.env.variable(name), is.list(list_content))
  structure(
    list(name = name, list_content = list_content),
    class = "JobFiles"
  )
}

###############################################
#' Make one temporary file from R-side content
#'
#' The generated path is stored in an environment variable for later use through
#' `JobEnv(name)`.
#'
#' @param name Environment variable that will hold the file path.
#' @param lines File contents.
#'
#' @return A `JobFile` step.
#' @noRd
JobFile <- function(name, lines) {
  stopifnot(is.valid.env.variable(name))
  structure(list(name = name, lines = lines), class = "JobFile")
}

###############################################
#' Set an environment variable in the job
#'
#' Local execution restores the previous environment value after the task.
#' Bash rendering emits an `export` line.
#'
#' @param name Environment variable name.
#' @param value Value to assign. May be a literal value or a job placeholder.
#'
#' @return A `JobSetEnv` step.
#' @noRd
JobSetEnv <- function(name, value) {
  stopifnot(is.valid.env.variable(name))
  structure(list(name = name, value = value), class = "JobSetEnv")
}

###############################################
#' Ensure a directory exists
#'
#' Local execution calls `dir.create(..., recursive = TRUE)`. Bash rendering
#' emits `mkdir -p`.
#'
#' @param path Directory path value. May be a literal path or a job placeholder.
#'
#' @return A `JobEnsureDir` step.
#' @noRd
JobEnsureDir <- function(path) {
  structure(list(path = path), class = "JobEnsureDir")
}

###############################################
#' Build an external command step
#'
#' `JobCommand` is for non-bascet tools. Prefer more specific typed steps when
#' the operation is known, such as `JobEnsureDir` for directory creation.
#'
#' @param command Command name or executable path.
#' @param args Command arguments. Elements may be literals, `JobArg` objects, or
#'   job placeholders.
#' @param prepend Command prefix, commonly a simple container launcher. Local
#'   execution splits this on whitespace; complex shell syntax is not supported.
#'
#' @return A `JobCommand` step.
#' @noRd
JobCommand <- function(command, args = list(), prepend = "") {
  structure(
    list(command = command, args = Filter(Negate(is.null), args), prepend = prepend),
    class = "JobCommand"
  )
}

###############################################
#' Build a bascet command step
#'
#' Adds bascet executable path and logging flags from `bascetInstance`, then
#' appends `args`.
#'
#' @param bascetInstance Bascet instance.
#' @param args Bascet command arguments. Elements may be literals, `JobArg`
#'   objects, or job placeholders.
#'
#' @return A `JobBascetCommand` step.
#' @noRd
JobBascetCommand <- function(bascetInstance, args = list()) {
  stopifnot(is.bascet.instance(bascetInstance))
  structure(
    list(bascetInstance = bascetInstance, args = Filter(Negate(is.null), args)),
    class = "JobBascetCommand"
  )
}

###############################################
#' Build a pipeline step
#'
#' Pipeline steps connect stdout from one command to stdin of the next. Local
#' execution currently captures intermediate output in memory.
#'
#' @param ... Command-like steps to connect in order.
#'
#' @return A `JobPipeline` step.
#' @noRd
JobPipeline <- function(...) {
  steps <- list(...)
  stopifnot(length(steps) > 0)
  structure(list(steps = steps), class = "JobPipeline")
}

###############################################
#' Redirect step output
#'
#' At present only stdout redirection is modeled.
#'
#' @param step Step whose output should be redirected.
#' @param stdout Output path value, or NULL for no redirection.
#'
#' @return A `JobRedirect` step.
#' @noRd
JobRedirect <- function(step, stdout = NULL) {
  structure(list(step = step, stdout = stdout), class = "JobRedirect")
}

###############################################
#' Render a JobScript to bash lines
#'
#' This is the backend used by SLURM. Local execution does not use this except
#' for optional debug display.
#'
#' @param script JobScript object.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector of bash script lines.
#' @noRd
renderJobScriptBash <- function(script, task_id = "$TASK_ID") {
  stopifnot(is.jobscript(script))
  c(
    renderJobVarsBash(script$vars),
    unlist(lapply(script$steps, renderJobStepBash, task_id = task_id), use.names = FALSE)
  )
}

###############################################
#' Render JobScript variables as bash arrays
#'
#' @param vars Named list of per-task values.
#'
#' @return Character vector of bash array assignments.
#' @noRd
renderJobVarsBash <- function(vars) {
  if(length(vars) == 0) {
    return(character())
  }

  unlist(Map(function(name, vals) {
    paste0(
      name,
      "=(",
      stringr::str_flatten(shQuote(as.character(vals)), collapse = " "),
      ")"
    )
  }, names(vars), vars), use.names = FALSE)
}

###############################################
#' Render one job value for bash
#'
#' @param x Literal value or job placeholder.
#' @param task_id Bash expression for the current task id.
#' @param quote_literal Whether literal values should be shell-quoted.
#'
#' @return Character scalar suitable for bash output.
#' @noRd
renderJobValueBash <- function(x, task_id = "$TASK_ID", quote_literal = TRUE) {
  if(is.null(x)) {
    return(character())
  }
  if(inherits(x, "JobVar")) {
    return(paste0("${", x$name, "[", task_id, "]}"))
  }
  if(inherits(x, "JobEnv")) {
    return(paste0("$", x$name))
  }
  if(inherits(x, "JobFileLinesCsv")) {
    path <- renderJobValueBash(x$path, task_id = task_id, quote_literal = FALSE)
    return(paste0("$(cat ", path, " | paste --serial -d, - -)"))
  }

  val <- as.character(x)
  if(quote_literal) {
    shQuote(val)
  } else {
    val
  }
}

###############################################
#' Render one command argument for bash
#'
#' @param arg Literal value, job placeholder, or `JobArg`.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector of one or more bash argv fragments.
#' @noRd
renderJobArgBash <- function(arg, task_id = "$TASK_ID") {
  if(inherits(arg, "JobArg")) {
    if(is.null(arg$flag)) {
      return(renderJobValueBash(arg$value, task_id = task_id))
    }
    if(is.null(arg$value)) {
      return(arg$flag)
    }

    value <- renderJobValueBash(arg$value, task_id = task_id)
    if(arg$sep == "=") {
      paste0(arg$flag, "=", value)
    } else {
      c(arg$flag, value)
    }
  } else {
    renderJobValueBash(arg, task_id = task_id)
  }
}

###############################################
#' Render command arguments for bash
#'
#' @param args List of literal values, job placeholders, or `JobArg` objects.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector of bash argv fragments.
#' @noRd
renderJobArgsBash <- function(args, task_id = "$TASK_ID") {
  unlist(lapply(args, renderJobArgBash, task_id = task_id), use.names = FALSE)
}

###############################################
#' Render a JobScript step for bash
#'
#' Generic used by step-specific render methods.
#'
#' @param step JobScript step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector of bash script lines.
#' @noRd
renderJobStepBash <- function(step, task_id = "$TASK_ID") {
  UseMethod("renderJobStepBash")
}

###############################################
#' Render a skip-if-file-exists step for bash
#'
#' @param step `JobSkipIfFileExists` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector containing an early-exit shell block.
#' @noRd
renderJobStepBash.JobSkipIfFileExists <- function(step, task_id = "$TASK_ID") {
  path <- renderJobValueBash(step$path, task_id = task_id, quote_literal = FALSE)
  c(
    paste0("if [ -f ", path, " ]; then"),
    "  echo \"Skipping job as the output exists already\"",
    "  exit 0",
    "fi"
  )
}

###############################################
#' Render per-task temporary files for bash
#'
#' @param step `JobFiles` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector creating a bash array of generated file paths and a
#'   cleanup trap.
#' @noRd
renderJobStepBash.JobFiles <- function(step, task_id = "$TASK_ID") {
  # This intentionally mirrors the old file expander behavior. Files are created
  # from R before job submission, then exposed as per-task shell arrays.
  ts <- as.character(format(as.numeric(Sys.time()) * 1000, scientific = FALSE))
  path_tmp <- file.path(".jobdata")
  dir.create(path_tmp, showWarnings = FALSE)

  all_files <- character(length(step$list_content))
  for(i in seq_along(step$list_content)) {
    one_tmp <- file.path(path_tmp, paste0(ts, "-", i))
    if(file.exists(one_tmp)) {
      stop(paste0("Attempting to create file ", one_tmp, ", but it already exists"))
    }
    writeLines(con = one_tmp, text = step$list_content[[i]])
    all_files[i] <- one_tmp
  }

  c(
    renderJobVarsBash(setNames(list(all_files), step$name)),
    paste0("trap \"rm -rf ${", step$name, "[", task_id, "]}\" EXIT")
  )
}

###############################################
#' Render a single temporary file for bash
#'
#' @param step `JobFile` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character vector that creates a temp file and writes content to it.
#' @noRd
renderJobStepBash.JobFile <- function(step, task_id = "$TASK_ID") {
  c(
    paste0(step$name, "=$(mktemp)"),
    paste0(
      "echo -e ",
      shQuote(stringr::str_flatten(step$lines, collapse = "\\n")),
      " > ${", step$name, "}"
    )
  )
}

###############################################
#' Render an environment assignment for bash
#'
#' @param step `JobSetEnv` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character scalar `export` command.
#' @noRd
renderJobStepBash.JobSetEnv <- function(step, task_id = "$TASK_ID") {
  paste0(
    "export ",
    step$name,
    "=",
    renderJobValueBash(step$value, task_id = task_id)
  )
}

###############################################
#' Render a directory creation step for bash
#'
#' @param step `JobEnsureDir` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character scalar `mkdir -p` command.
#' @noRd
renderJobStepBash.JobEnsureDir <- function(step, task_id = "$TASK_ID") {
  paste("mkdir -p", renderJobValueBash(step$path, task_id = task_id))
}

###############################################
#' Render an external command for bash
#'
#' @param step `JobCommand` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character scalar command line.
#' @noRd
renderJobStepBash.JobCommand <- function(step, task_id = "$TASK_ID") {
  stringr::str_flatten(
    c(
      if(nzchar(step$prepend)) step$prepend,
      step$command,
      renderJobArgsBash(step$args, task_id = task_id)
    ),
    collapse = " "
  )
}

###############################################
#' Render a bascet command for bash
#'
#' @param step `JobBascetCommand` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character scalar command line.
#' @noRd
renderJobStepBash.JobBascetCommand <- function(step, task_id = "$TASK_ID") {
  bi <- step$bascetInstance
  stringr::str_flatten(
    c(
      if(nzchar(bi@prependCmd)) bi@prependCmd,
      shQuote(bi@bin),
      if(bi@logToFile) "--log-mode=$BASCET_LOGFILE",
      paste0("--log-level=", bi@logLevel),
      renderJobArgsBash(step$args, task_id = task_id)
    ),
    collapse = " "
  )
}

###############################################
#' Render a pipeline for bash
#'
#' @param step `JobPipeline` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character scalar pipeline command.
#' @noRd
renderJobStepBash.JobPipeline <- function(step, task_id = "$TASK_ID") {
  stringr::str_flatten(
    vapply(step$steps, renderJobStepBash, character(1), task_id = task_id),
    collapse = " | "
  )
}

###############################################
#' Render output redirection for bash
#'
#' @param step `JobRedirect` step.
#' @param task_id Bash expression for the current task id.
#'
#' @return Character scalar command line, with redirection if requested.
#' @noRd
renderJobStepBash.JobRedirect <- function(step, task_id = "$TASK_ID") {
  rendered <- renderJobStepBash(step$step, task_id = task_id)
  if(!is.null(step$stdout)) {
    stdout <- renderJobValueBash(step$stdout, task_id = task_id, quote_literal = FALSE)
    rendered <- paste(rendered, ">", stdout)
  }
  rendered
}


################################################################################
################ Direct local execution ########################################
################################################################################

###############################################
#' Run one task from a JobScript directly from R
#' @param script JobScript object
#' @param task_id Zero-based task id
#' @param env Named list of environment variables for this task
#' @return TRUE if the task ran, FALSE if it was skipped
#' @noRd
runJobScriptLocal <- function(script, task_id, env = list()) {
  stopifnot(is.jobscript(script))

  ctx <- new.env(parent = emptyenv())
  ctx$script <- script
  ctx$task_id <- task_id
  ctx$task_index <- task_id + 1L
  ctx$vars <- script$vars
  ctx$cleanup <- character()
  ctx$env_names <- character()
  ctx$old_env <- list()

  on.exit({
    for(path in ctx$cleanup) {
      if(file.exists(path)) {
        unlink(path, recursive = TRUE, force = TRUE)
      }
    }
    restoreJobEnvLocal(ctx)
  }, add = TRUE)

  for(nm in names(env)) {
    setJobEnvLocal(ctx, nm, env[[nm]])
  }

  for(step in script$steps) {
    ret <- runJobStepLocal(step, ctx)
    if(identical(ret, "skip")) {
      return(FALSE)
    }
  }

  TRUE
}

###############################################
#' Set an environment variable for local task execution
#'
#' Records the previous value in the task context so it can be restored after
#' the task completes.
#'
#' @param ctx Local execution context.
#' @param name Environment variable name.
#' @param value Value to assign. NULL or NA unsets the variable.
#'
#' @return Invisibly returns NULL.
#' @noRd
setJobEnvLocal <- function(ctx, name, value) {
  stopifnot(is.valid.env.variable(name))
  if(!(name %in% ctx$env_names)) {
    ctx$env_names <- c(ctx$env_names, name)
    old <- Sys.getenv(name, unset = NA_character_)
    ctx$old_env[[name]] <- old
  }
  if(is.na(value) || is.null(value)) {
    Sys.unsetenv(name)
  } else {
    do.call(Sys.setenv, as.list(stats::setNames(as.character(value), name)))
  }
}

###############################################
#' Restore environment variables changed by local execution
#'
#' @param ctx Local execution context.
#'
#' @return Invisibly returns NULL.
#' @noRd
restoreJobEnvLocal <- function(ctx) {
  for(name in rev(ctx$env_names)) {
    old <- ctx$old_env[[name]]
    if(length(old) == 0 || is.na(old)) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(old, name)))
    }
  }
}

###############################################
#' Resolve a job value for local execution
#'
#' Converts placeholders such as `JobVar`, `JobEnv`, and `JobFileLinesCsv` into
#' concrete values for the current task.
#'
#' @param x Literal value or job placeholder.
#' @param ctx Local execution context.
#'
#' @return Character value, except list-backed `JobVar` values may return the
#'   original list element.
#' @noRd
resolveJobValueLocal <- function(x, ctx) {
  if(is.null(x)) {
    return(character())
  }
  if(inherits(x, "JobVar")) {
    vals <- ctx$vars[[x$name]]
    if(is.null(vals)) {
      stop(paste0("Unknown JobScript variable: ", x$name))
    }
    if(is.list(vals) && !is.data.frame(vals)) {
      return(vals[[ctx$task_index]])
    }
    if(length(vals) == 1) {
      return(as.character(vals[[1]]))
    }
    return(as.character(vals[[ctx$task_index]]))
  }
  if(inherits(x, "JobEnv")) {
    return(Sys.getenv(x$name))
  }
  if(inherits(x, "JobFileLinesCsv")) {
    path <- resolveJobValueLocal(x$path, ctx)
    if(!file.exists(path)) {
      stop(paste0("Cannot read missing file for JobFileLinesCsv: ", path))
    }
    return(stringr::str_flatten(readLines(path, warn = FALSE), collapse = ","))
  }

  as.character(x)
}

###############################################
#' Resolve one command argument for local execution
#'
#' @param arg Literal value, job placeholder, or `JobArg`.
#' @param ctx Local execution context.
#'
#' @return Character vector of one or more argv elements.
#' @noRd
resolveJobArgLocal <- function(arg, ctx) {
  if(inherits(arg, "JobArg")) {
    if(is.null(arg$flag)) {
      return(resolveJobValueLocal(arg$value, ctx))
    }
    if(is.null(arg$value)) {
      return(arg$flag)
    }
    value <- resolveJobValueLocal(arg$value, ctx)
    if(arg$sep == "=") {
      paste0(arg$flag, "=", value)
    } else {
      c(arg$flag, value)
    }
  } else {
    resolveJobValueLocal(arg, ctx)
  }
}

###############################################
#' Resolve command arguments for local execution
#'
#' @param args List of literal values, job placeholders, or `JobArg` objects.
#' @param ctx Local execution context.
#'
#' @return Character vector of argv elements.
#' @noRd
resolveJobArgsLocal <- function(args, ctx) {
  unlist(lapply(args, resolveJobArgLocal, ctx = ctx), use.names = FALSE)
}

###############################################
#' Split a command prefix for local execution
#'
#' This supports simple prefixes such as `singularity run image.sif`. Complex
#' shell syntax is not modeled by local direct execution.
#'
#' @param prepend Prefix string.
#'
#' @return Character vector of prefix argv elements.
#' @noRd
splitPrependLocal <- function(prepend) {
  if(is.null(prepend) || !nzchar(prepend)) {
    return(character())
  }
  strsplit(prepend, "[[:space:]]+")[[1]]
}

###############################################
#' Build executable and argv for local execution
#'
#' Applies `prepend` by turning the prefix command into the executable and
#' appending the requested command plus its arguments.
#'
#' @param command Command name or executable path.
#' @param args Command arguments.
#' @param prepend Command prefix.
#' @param ctx Local execution context.
#'
#' @return List with `command` and `args` elements.
#' @noRd
buildJobCommandLocal <- function(command, args, prepend, ctx) {
  resolved_args <- resolveJobArgsLocal(args, ctx)
  prepend_parts <- splitPrependLocal(prepend)

  if(length(prepend_parts) == 0) {
    list(command = command, args = resolved_args)
  } else {
    list(
      command = prepend_parts[[1]],
      args = c(prepend_parts[-1], command, resolved_args)
    )
  }
}

###############################################
#' Run a command with `system2`
#'
#' Centralizes status checking for direct local execution. When `capture` is
#' TRUE, stdout is returned for use by pipeline and redirect steps.
#'
#' @param command Command name or executable path.
#' @param args Character vector of argv elements.
#' @param input Optional stdin content.
#' @param capture Whether to capture stdout.
#' @param stdout `system2` stdout setting when not capturing.
#' @param stderr `system2` stderr setting.
#'
#' @return Captured stdout when `capture` is TRUE; otherwise invisible
#'   character vector.
#' @noRd
runSystem2Local <- function(command, args, input = NULL, capture = FALSE, stdout = "", stderr = "") {
  if(capture) {
    out <- tryCatch(
      system2(command, args = args, input = input, stdout = TRUE, stderr = stderr),
      warning = function(w) {
        value <- conditionMessage(w)
        attr(value, "status") <- 1L
        value
      }
    )
    status <- attr(out, "status")
    if(is.null(status)) {
      status <- 0L
    }
    if(status != 0L) {
      stop(paste0("Command failed with exit status ", status, ": ", command))
    }
    return(out)
  }

  status <- system2(command, args = args, input = input, stdout = stdout, stderr = stderr)
  if(status != 0L) {
    stop(paste0("Command failed with exit status ", status, ": ", command))
  }
  invisible(character())
}

###############################################
#' Run one JobScript step locally
#'
#' Generic used by step-specific local execution methods.
#'
#' @param step JobScript step.
#' @param ctx Local execution context.
#' @param input Optional stdin content for command-like steps.
#' @param capture Whether command-like output should be captured.
#'
#' @return Step-specific result. Skip steps may return the sentinel `"skip"`.
#' @noRd
runJobStepLocal <- function(step, ctx, input = NULL, capture = FALSE) {
  UseMethod("runJobStepLocal")
}

###############################################
#' Run skip-if-file-exists locally
#'
#' @param step `JobSkipIfFileExists` step.
#' @param ctx Local execution context.
#' @param input Ignored.
#' @param capture Ignored.
#'
#' @return `"skip"` when the file exists; otherwise invisible NULL.
#' @noRd
runJobStepLocal.JobSkipIfFileExists <- function(step, ctx, input = NULL, capture = FALSE) {
  path <- resolveJobValueLocal(step$path, ctx)
  if(file.exists(path)) {
    message("Skipping job as the output exists already")
    return("skip")
  }
  invisible(NULL)
}

###############################################
#' Run per-task temporary file creation locally
#'
#' Creates only the current task's file and stores the path in the local context
#' variables under `step$name`.
#'
#' @param step `JobFiles` step.
#' @param ctx Local execution context.
#' @param input Ignored.
#' @param capture Ignored.
#'
#' @return Invisibly returns NULL.
#' @noRd
runJobStepLocal.JobFiles <- function(step, ctx, input = NULL, capture = FALSE) {
  path_tmp <- file.path(".jobdata")
  dir.create(path_tmp, showWarnings = FALSE)
  one_tmp <- tempfile(pattern = paste0(step$name, "."), tmpdir = path_tmp)
  writeLines(con = one_tmp, text = step$list_content[[ctx$task_index]])
  ctx$vars[[step$name]] <- one_tmp
  ctx$cleanup <- c(ctx$cleanup, one_tmp)
  invisible(NULL)
}

###############################################
#' Run single temporary file creation locally
#'
#' Creates a temp file, writes `step$lines`, and stores the path in the named
#' environment variable.
#'
#' @param step `JobFile` step.
#' @param ctx Local execution context.
#' @param input Ignored.
#' @param capture Ignored.
#'
#' @return Invisibly returns NULL.
#' @noRd
runJobStepLocal.JobFile <- function(step, ctx, input = NULL, capture = FALSE) {
  one_tmp <- tempfile(pattern = paste0(step$name, "."))
  writeLines(con = one_tmp, text = step$lines)
  setJobEnvLocal(ctx, step$name, one_tmp)
  ctx$cleanup <- c(ctx$cleanup, one_tmp)
  invisible(NULL)
}

###############################################
#' Run environment assignment locally
#'
#' @param step `JobSetEnv` step.
#' @param ctx Local execution context.
#' @param input Ignored.
#' @param capture Ignored.
#'
#' @return Invisibly returns NULL.
#' @noRd
runJobStepLocal.JobSetEnv <- function(step, ctx, input = NULL, capture = FALSE) {
  setJobEnvLocal(ctx, step$name, resolveJobValueLocal(step$value, ctx))
  invisible(NULL)
}

###############################################
#' Run directory creation locally
#'
#' @param step `JobEnsureDir` step.
#' @param ctx Local execution context.
#' @param input Ignored.
#' @param capture Ignored.
#'
#' @return Invisibly returns NULL.
#' @noRd
runJobStepLocal.JobEnsureDir <- function(step, ctx, input = NULL, capture = FALSE) {
  dir.create(resolveJobValueLocal(step$path, ctx), recursive = TRUE, showWarnings = FALSE)
  invisible(NULL)
}

###############################################
#' Run an external command locally
#'
#' @param step `JobCommand` step.
#' @param ctx Local execution context.
#' @param input Optional stdin content.
#' @param capture Whether stdout should be captured.
#'
#' @return Captured stdout when requested; otherwise invisible character vector.
#' @noRd
runJobStepLocal.JobCommand <- function(step, ctx, input = NULL, capture = FALSE) {
  cmd <- buildJobCommandLocal(step$command, step$args, step$prepend, ctx)
  runSystem2Local(cmd$command, cmd$args, input = input, capture = capture)
}

###############################################
#' Run a bascet command locally
#'
#' Adds bascet logging flags and resolves all command arguments before invoking
#' the bascet executable.
#'
#' @param step `JobBascetCommand` step.
#' @param ctx Local execution context.
#' @param input Optional stdin content.
#' @param capture Whether stdout should be captured.
#'
#' @return Captured stdout when requested; otherwise invisible character vector.
#' @noRd
runJobStepLocal.JobBascetCommand <- function(step, ctx, input = NULL, capture = FALSE) {
  bi <- step$bascetInstance
  args <- c(
    if(bi@logToFile) paste0("--log-mode=", Sys.getenv("BASCET_LOGFILE")),
    paste0("--log-level=", bi@logLevel),
    resolveJobArgsLocal(step$args, ctx)
  )
  cmd <- buildJobCommandLocal(bi@bin, args, bi@prependCmd, ctx)
  runSystem2Local(cmd$command, cmd$args, input = input, capture = capture)
}

###############################################
#' Run a pipeline locally
#'
#' Executes each command-like step in order, passing captured stdout as stdin to
#' the next step. Intermediate output is stored in memory.
#'
#' @param step `JobPipeline` step.
#' @param ctx Local execution context.
#' @param input Optional stdin content for the first step.
#' @param capture Whether final stdout should be returned.
#'
#' @return Captured final stdout when requested; otherwise invisible NULL.
#' @noRd
runJobStepLocal.JobPipeline <- function(step, ctx, input = NULL, capture = FALSE) {
  cur_input <- input
  for(i in seq_along(step$steps)) {
    cur_input <- runJobStepLocal(step$steps[[i]], ctx, input = cur_input, capture = TRUE)
  }
  if(capture) {
    cur_input
  } else {
    if(length(cur_input) > 0) {
      writeLines(cur_input)
    }
    invisible(NULL)
  }
}

###############################################
#' Run output redirection locally
#'
#' Captures the wrapped step output and writes it to the resolved stdout path.
#'
#' @param step `JobRedirect` step.
#' @param ctx Local execution context.
#' @param input Optional stdin content.
#' @param capture Whether output should be captured when no redirection is
#'   configured.
#'
#' @return Invisibly returns NULL when redirecting; otherwise returns the wrapped
#'   step result.
#' @noRd
runJobStepLocal.JobRedirect <- function(step, ctx, input = NULL, capture = FALSE) {
  if(!is.null(step$stdout)) {
    out <- runJobStepLocal(step$step, ctx, input = input, capture = TRUE)
    writeLines(out, con = resolveJobValueLocal(step$stdout, ctx))
    return(invisible(NULL))
  }
  runJobStepLocal(step$step, ctx, input = input, capture = capture)
}

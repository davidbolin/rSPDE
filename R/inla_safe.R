#' Load INLA safely for examples and tests
#'
#' Checks that a sufficiently recent INLA package and matching executable are
#' available. In non-interactive use it also limits INLA to one thread for
#' repeatability. An optional Cgeneric symbol can be required when an rSPDE
#' feature depends on a model that may not yet be included in every INLA build.
#'
#' @param multicore Logical; if `TRUE`, leave INLA's `num.threads` option
#'   unchanged. If `FALSE`, set it to `"1:1:1"`. The default allows multicore
#'   only in interactive sessions outside testthat.
#' @param quietly Logical; suppress diagnostic messages when `TRUE`.
#' @param minimum_version Minimum acceptable INLA version. This should match
#'   the requirement in `DESCRIPTION`.
#' @param required_symbol Optional name of a required INLA Cgeneric symbol.
#'
#' @return `TRUE` when INLA is available and safe to use; otherwise `FALSE`.
#' @export
#' @keywords internal
#'
#' @examples
#' if (rspde_safe_inla()) {
#'   # Run INLA-dependent calculations.
#' }
rspde_safe_inla <- function(multicore = NULL,
                            quietly = FALSE,
                            minimum_version = "24.12.01",
                            required_symbol = NULL) {
  inla_version <- tryCatch(
    utils::packageVersion("INLA"),
    error = function(e) NULL
  )
  if (is.null(inla_version)) {
    if (!quietly) message("Package 'INLA' is not installed.")
    return(FALSE)
  }
  if (inla_version < numeric_version(minimum_version)) {
    if (!quietly) {
      message(
        "Installed INLA version is ", inla_version,
        "; version >= ", minimum_version, " is required."
      )
    }
    return(FALSE)
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    if (!quietly) message("Package 'INLA' could not be loaded safely.")
    return(FALSE)
  }

  inla_call <- tryCatch(
    INLA::inla.getOption("inla.call"),
    error = function(e) NULL
  )
  if (!is.character(inla_call) || length(inla_call) < 1L) {
    if (!quietly) {
      message("inla.getOption('inla.call') failed; INLA is not installed correctly.")
    }
    return(FALSE)
  }
  run_wrappers <- grepl("\\.run$", inla_call)
  if (all(run_wrappers)) {
    binaries <- sub("\\.run$", "", inla_call)
    if (!all(file.exists(binaries))) {
      if (!quietly) {
        message(
          "The INLA executable was not found; INLA may have a platform mismatch."
        )
      }
      return(FALSE)
    }
  }

  if (!is.null(required_symbol) &&
      !rspde_cgeneric_symbol_available(required_symbol)) {
    if (!quietly) {
      message(
        "The required Cgeneric symbol '", required_symbol,
        "' is unavailable in both INLA and the local rSPDE shared library."
      )
    }
    return(FALSE)
  }

  if (is.null(multicore)) {
    multicore <- interactive() &&
      !identical(Sys.getenv("TESTTHAT"), "true")
  }
  if (!multicore) {
    threads <- tryCatch(
      INLA::inla.getOption("num.threads"),
      error = function(e) NULL
    )
    if (is.null(threads)) {
      if (!quietly) message("Could not read INLA's num.threads option.")
      return(FALSE)
    }
    if (!identical(threads, "1:1:1")) {
      if (!quietly) {
        message(
          "Changing INLA option num.threads from '", threads,
          "' to '1:1:1'."
        )
      }
      INLA::inla.setOption(num.threads = "1:1:1")
    }
  }
  TRUE
}

#' @rdname rspde_safe_inla
#' @param envir Environment in which test cleanup handlers are registered.
#' @return `local_rspde_safe_inla()` is called for its testthat skip and local
#'   option-setting side effects.
#' @export
local_rspde_safe_inla <- function(multicore = FALSE,
                                  quietly = TRUE,
                                  required_symbol = NULL,
                                  envir = parent.frame()) {
  available <- rspde_safe_inla(
    multicore = multicore,
    quietly = quietly,
    required_symbol = required_symbol
  )
  skip_message <- if (is.null(required_symbol)) {
    "INLA is not available safely."
  } else {
    paste0("INLA Cgeneric model '", required_symbol, "' is not available.")
  }
  testthat::skip_if_not(available, skip_message)

  if (requireNamespace("INLA", quietly = TRUE)) {
    option_names <- c(
      "num.threads", "inla.timeout", "fmesher.evolution.warn",
      "fmesher.evolution.verbosity"
    )
    old_options <- lapply(option_names, function(name) {
      tryCatch(INLA::inla.getOption(name), error = function(e) e)
    })
    if (any(vapply(old_options, inherits, logical(1), what = "error"))) {
      testthat::skip("INLA options could not be read safely.")
    }
    names(old_options) <- option_names
    withr::defer(do.call(INLA::inla.setOption, old_options), envir = envir)
  }

  INLA::inla.setOption(
    inla.timeout = 60,
    fmesher.evolution.warn = TRUE,
    fmesher.evolution.verbosity = "stop"
  )
  invisible(TRUE)
}

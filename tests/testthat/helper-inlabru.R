# inlabru::bru() stops when inlabru cannot use the installed INLA. In
# particular, inlabru <= 2.15.0 expects INLA's 'inla.call' option to hold the
# path to the INLA executable, but INLA >= 26.11 leaves that option unset. In
# that case, point 'inla.call' to the INLA executable until the calling test
# ends. Skip the test if inlabru still cannot use INLA.
skip_if_inlabru_cannot_use_inla <- function(envir = parent.frame()) {
  testthat::skip_if_not_installed("inlabru")
  inlabru_can_use_inla <- function() {
    tryCatch(
      isTRUE(inlabru::bru_safe_inla(multicore = TRUE, quietly = TRUE)),
      error = function(e) FALSE
    )
  }
  if (!inlabru_can_use_inla() &&
      requireNamespace("INLA", quietly = TRUE) &&
      is.null(INLA::inla.getOption("inla.call"))) {
    inla_call <- rspde_inla_call()
    if (!is.null(inla_call)) {
      INLA::inla.setOption(inla.call = inla_call)
      withr::defer(INLA::inla.setOption(inla.call = NULL), envir = envir)
    }
  }
  testthat::skip_if_not(
    inlabru_can_use_inla(),
    "inlabru cannot use the installed INLA."
  )
}

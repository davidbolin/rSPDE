# inlabru::bru() stops when inlabru cannot use the installed INLA. This
# happens, for example, when inlabru's own INLA check does not support a newer
# INLA version. Skip tests that call bru() in that case.
skip_if_inlabru_cannot_use_inla <- function() {
  testthat::skip_if_not_installed("inlabru")
  ok <- tryCatch(
    isTRUE(inlabru::bru_safe_inla(multicore = TRUE, quietly = TRUE)),
    error = function(e) FALSE
  )
  testthat::skip_if_not(ok, "inlabru cannot use the installed INLA.")
}

# rspde_lme() warns when an optimisation method cannot produce a numerically
# positive-definite Hessian and it falls back to another method or result.
# In the tests this is an optimisation nuisance, not a failure, so these
# helpers muffle only that warning. Every other warning is still reported.
muffle_hessian_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("positive-definite Hessian", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

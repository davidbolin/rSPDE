if (requireNamespace("inlabru", quietly = TRUE)) {
  inlabru::bru_options_set_local(
    "bru_compat_pre_2_14_enable" = FALSE,
    .envir = testthat::teardown_env()
  )
}

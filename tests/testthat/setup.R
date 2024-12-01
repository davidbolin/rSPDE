if (requireNamespace("INLA", quietly = TRUE)) {
  INLA::inla.setOption(fmesher.evolution.warn = TRUE)
  INLA::inla.setOption(fmesher.evolution.verbosity = "stop")
}

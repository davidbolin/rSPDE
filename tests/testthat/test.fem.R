context("FEM calculations")

test_that("FEM matrix construction", {
    testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

  x <- c(0,1,1.75,2)
  fem <- rSPDE.fem1d(x)
  mesh <- fmesher::fm_mesh_1d(x)
  fem_inla <- fmesher::fm_fem(mesh)

  # Pr multiplication
  expect_equal(fem$C, fem_inla$c1, tolerance = 1e-10)
  expect_equal(fem$G, fem_inla$g1, tolerance = 1e-10)

    INLA::inla.setOption(num.threads = old_threads)
})

test_that("A matrix construction", {
    testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }

  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")

  x <- c(0,1,1.75,2)
  fem <- rSPDE.fem1d(x)
  loc <- c(0.6,1.1,2)
  A <- rSPDE.A1d(x,loc)
  mesh <- fmesher::fm_mesh_1d(x)
  A_inla <- fmesher::fm_basis(mesh,loc)

  # Pr multiplication
  expect_equal(A, A_inla, tolerance = 1e-10)

    INLA::inla.setOption(num.threads = old_threads)
})

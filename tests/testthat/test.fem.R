context("FEM calculations")

test_that("FEM matrix construction", {
  x <- c(0,1,1.75,2)
  fem <- rSPDE.fem1d(x)
  mesh <- INLA::inla.mesh.1d(x)
  fem_inla <- INLA::inla.mesh.fem(mesh)
  
  # Pr multiplication
  expect_equal(fem$C, fem_inla$c1, tolerance = 1e-10)
  expect_equal(fem$G, fem_inla$g1, tolerance = 1e-10)
})

test_that("A matrix construction", {
  x <- c(0,1,1.75,2)
  fem <- rSPDE.fem1d(x)
  loc <- c(0.6,1.1,2)
  A <- rSPDE.A1d(x,loc)
  mesh <- INLA::inla.mesh.1d(x)
  A_inla <- INLA::inla.mesh.1d.A(mesh,loc)
  
  # Pr multiplication
  expect_equal(A, A_inla, tolerance = 1e-10)
})
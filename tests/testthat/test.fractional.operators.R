context("factional.operators")

test_that("Operator construction for fractional stationary Matern", {
  
  x <- seq(from = 0, to = 1, length.out = 51)
  fem <- rSPDE.fem1d(x)
  
  d <- 1
  nu <- 0.8
  sigma <- 0.5
  kappa <- 20
  
  op1 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                              G = fem$G, C = fem$C, d = 1,
                          type="operator")
  
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(d /2) * gamma(nu + d/2)))
  beta <- (nu + d/2)/2
  
  op2 <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu,
                              G = fem$G, C = fem$C, d = d)
  
  L <- fem$G + kappa^2*fem$C
  op3 <- fractional.operators(L = L, scale.factor = kappa^2, tau = tau, 
                              beta = beta, C = fem$C)
  v <- t(rSPDE.A1d(x,0.5))
  c1 <- as.vector(Sigma.mult(op1,v))
  c2 <- as.vector(Sigma.mult(op2,v))
  c3 <- as.vector(Sigma.mult(op3,v))
  c0 <- as.vector(matern.covariance(abs(x-0.5),kappa = kappa, nu = nu, sigma = sigma))
  
  expect_equal(c1, c2, tolerance = 1e-10)
  expect_equal(c2, c3, tolerance = 1e-10)
  expect_equal(c3, c0, tolerance = 0.02)
  
})

test_that("Operator construction for non-fractional stationary Matern", {
  
  x <- seq(from = 0, to = 1, length.out = 51)
  fem <- rSPDE.fem1d(x)
  
  d <- 1
  nu <- 1.5
  sigma <- 0.5
  kappa <- 20
  
  op1 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                          G = fem$G, C = fem$C, d = 1,
                          type="operator")
  
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(d /2) * gamma(nu + d/2)))
  beta <- (nu + d/2)/2
  
  op2 <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu,
                               G = fem$G, C = fem$C, d = d)
  
  L <- fem$G + kappa^2*fem$C
  op3 <- fractional.operators(L = L, scale.factor = kappa^2, tau = tau, 
                              beta = beta, C = fem$C)
  v <- t(rSPDE.A1d(x,0.5))
  c1 <- as.vector(Sigma.mult(op1,v))
  c2 <- as.vector(Sigma.mult(op2,v))
  c3 <- as.vector(Sigma.mult(op3,v))
  c0 <- as.vector(matern.covariance(abs(x-0.5),kappa = kappa, nu = nu, sigma = sigma))
  
  expect_equal(c1, c2, tolerance = 1e-10)
  expect_equal(c2, c3, tolerance = 1e-10)
  expect_equal(c3, c0, tolerance = 0.02)
  
})


test_that("Operator construction for fractional stationary Matern with beta>1", {
  
  x <- seq(from = 0, to = 1, length.out = 51)
  fem <- rSPDE.fem1d(x)
  
  d <- 1
  nu <- 2
  sigma <- 0.5
  kappa <- 20
  
  op1 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                          G = fem$G, C = fem$C, d = 1,
                          type="operator")
  
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(d /2) * gamma(nu + d/2)))
  beta <- (nu + d/2)/2
  
  op2 <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu,
                               G = fem$G, C = fem$C, d = d)
  
  L <- fem$G + kappa^2*fem$C
  op3 <- fractional.operators(L = L, scale.factor = kappa^2, tau = tau, 
                              beta = beta, C = fem$C)
  v <- t(rSPDE.A1d(x,0.5))
  c1 <- as.vector(Sigma.mult(op1,v))
  c2 <- as.vector(Sigma.mult(op2,v))
  c3 <- as.vector(Sigma.mult(op3,v))
  c0 <- as.vector(matern.covariance(abs(x-0.5),kappa = kappa, nu = nu, sigma = sigma))
  
  expect_equal(c1, c2, tolerance = 1e-10)
  expect_equal(c2, c3, tolerance = 1e-10)
  expect_equal(c3, c0, tolerance = 0.02)
  
})


test_that("Operator construction for non-stationary Matern", {
  
  x <- seq(from = 0, to = 1, length.out = 51)
  fem <- rSPDE.fem1d(x)
  
  d <- 1
  nu <- 0.8
  kappa <-  10*(1+2*x^2)
  tau <-  0.1*(1 - 0.7*x^2)
  op1 <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu, 
                              G = fem$G, C = fem$C, d = d, m=1)
  
  beta <- (nu + d/2)/2
  
  L <- fem$G + fem$C%*%Matrix::Diagonal(dim(fem$C)[1],kappa^2)
  op2 <- fractional.operators(L = L, scale.factor = min(kappa)^2, tau = tau, 
                              beta = beta, C = fem$C)
  v <- t(rSPDE.A1d(x,0.5))
  c1 <- as.vector(Sigma.mult(op1,v))
  c2 <- as.vector(Sigma.mult(op2,v))
  expect_equal(c1, c2, tolerance = 1e-10)
  
})

## Weighted-L2 rational coefficients in the INLA interface: the precision the
## cgeneric model assembles, the block count that rspde.make.A() and
## rspde.make.index() have to follow, and the coefficient table the model is
## given. Everything else is in test.rational.wl2.R.
##
## The cgeneric tests need the compiled shared library, which only exists in an
## installed package:
##     RSPDE_COMPILE=1 R CMD INSTALL .
##     NOT_CRAN=true Rscript -e 'library(rSPDE); testthat::test_dir("tests/testthat")'
## Under devtools::load_all() the library is not on disk and they skip.

test_that("the precision builders agree for the weighted-L2 type", {
  ## rspde.matern.precision() assembles the blocks from Malpha and Malpha2,
  ## while the blocks stored in the model are assembled from L and C^-1; the
  ## two forms must agree. rspde.matern.precision.opt() is the vectorised
  ## version of the same expression.
  x <- seq(from = 0, to = 1, length.out = 61)
  n <- length(x)
  for (nu in c(0.25, 1, 1.9)) {
    op <- matern.operators(
      loc_mesh = x, nu = nu, range = 0.2, sigma = 1, d = 1, m = 2,
      parameterization = "matern", type_rational_approximation = "wl2"
    )
    fm <- op$fem_mesh_matrices
    Ci <- Matrix::Diagonal(n, 1 / diag(fm$c0))
    fm$g2 <- fm$g1 %*% Ci %*% fm$g1
    fm$g3 <- fm$g2 %*% Ci %*% fm$g1
    fm$g4 <- fm$g3 %*% Ci %*% fm$g1
    Q <- rspde.matern.precision(
      kappa = op$kappa, nu = nu, tau = op$tau, rspde.order = 2, dim = 1,
      fem_mesh_matrices = fm, type_rational_approx = "wl2",
      wl2_table = op$wl2_table
    )
    expect_equal(max(abs(Q - op$Q)) / max(abs(op$Q)), 0, tolerance = 1e-12)

    if (nu >= 1) {
      ## The sparsity grows with floor(alpha), so the pattern is the one of the
      ## highest power the blocks use.
      gr <- as.matrix(fm[[paste0("g", floor(nu + 0.5) + 1)]]) != 0
      pick <- function(M) as.matrix(M)[gr]
      qv <- rspde.matern.precision.opt(
        kappa = op$kappa, nu = nu, tau = op$tau, rspde.order = 2, dim = 1,
        fem_matrices = list(
          C = pick(fm$c0), G = pick(fm$g1), G_2 = pick(fm$g2),
          G_3 = pick(fm$g3), G_4 = pick(fm$g4)
        ),
        graph = NULL, sharp = TRUE, type_rational_approx = "wl2",
        wl2_table = op$wl2_table
      )
      ref <- c(
        pick(op$Q[seq_len(n), seq_len(n)]),
        pick(op$Q[n + seq_len(n), n + seq_len(n)])
      )
      expect_equal(max(abs(qv - ref)) / max(abs(ref)), 0, tolerance = 1e-12)
    }
  }
})

test_that("the INLA interface builds the weighted-L2 precision", {
  ## The cgeneric model assembles the blocks in C from the finite element
  ## matrices; precision() builds the same model in R through
  ## matern.operators(). The two must agree, for every class.
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_stat_frac_wl2_model"
  )
  set.seed(1)
  loc <- cbind(runif(60), runif(60))
  mesh <- fmesher::fm_mesh_2d_inla(
    loc = loc, max.edge = c(0.25, 0.5), cutoff = 0.08
  )
  for (nu in c(0.4, 0.9, 1.4)) {
    for (m in c(1, 2, 3)) {
      mod <- rspde.matern(
        mesh = mesh, nu = nu, rspde.order = m, parameterization = "matern",
        type.rational.approx = "wl2", shared_lib = "rSPDE",
        start.theta = c(log(1), log(0.4))
      )
      qq <- INLA::inla.cgeneric.q(mod)
      ## m blocks, not m + 1: there is no constant term
      expect_equal(nrow(qq$Q), m * mesh$n)
      ## precision() uses nu + 1e-10, so the agreement stops at about 1e-9
      expect_equal(
        max(abs(qq$Q - precision(mod))) / max(abs(precision(mod))), 0,
        tolerance = 1e-7
      )
    }
  }
})

test_that("the A matrix and the index have the weighted-L2 block count", {
  set.seed(1)
  loc <- cbind(runif(40), runif(40))
  mesh <- fmesher::fm_mesh_2d_inla(
    loc = loc, max.edge = c(0.3, 0.6), cutoff = 0.1
  )
  for (m in 1:3) {
    for (ty in c("brasil", "wl2")) {
      nb <- if (ty == "wl2") m else m + 1
      A <- rspde.make.A(
        mesh = mesh, loc = loc, rspde.order = m, nu = 0.9,
        type.rational.approx = ty
      )
      idx <- rspde.make.index(
        name = "f", mesh = mesh, rspde.order = m, nu = 0.9,
        type.rational.approx = ty
      )
      expect_equal(ncol(A), nb * mesh$n)
      expect_equal(length(idx$f), nb * mesh$n)
    }
  }
})

test_that("the INLA interface refuses the weighted-L2 type where it has none", {
  set.seed(1)
  mesh <- fmesher::fm_mesh_2d_inla(
    loc = cbind(runif(40), runif(40)), max.edge = c(0.3, 0.6), cutoff = 0.1
  )
  ## the parsimonious model is a different approximation
  expect_error(
    rspde.matern(
      mesh = mesh, nu = 0.9, rspde.order = 0, type.rational.approx = "wl2"
    ),
    "rspde.order >= 1"
  )
  ## floor(alpha) beyond 2 is out of scope for the classes
  expect_error(
    rspde.matern(
      mesh = mesh, rspde.order = 2, nu.upper.bound = 3,
      type.rational.approx = "wl2"
    ),
    "0, 1 or 2"
  )
})

test_that("the INLA interface estimates nu with the weighted-L2 type", {
  ## nu is estimated, so the coefficients reach the cgeneric model as a table
  ## with one block of 999 rows per floor(alpha). The lookup must pick the
  ## right block, so nu is checked on both sides of the band boundary.
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_stat_general_wl2_model"
  )
  set.seed(1)
  loc <- cbind(runif(60), runif(60))
  mesh <- fmesher::fm_mesh_2d_inla(
    loc = loc, max.edge = c(0.25, 0.5), cutoff = 0.08
  )
  for (m in c(1, 2, 3)) {
    for (nu in c(0.3, 0.9, 1.4, 1.9)) {
      mod <- rspde.matern(
        mesh = mesh, rspde.order = m, nu.upper.bound = 2,
        parameterization = "matern", type.rational.approx = "wl2",
        shared_lib = "rSPDE", start.theta = c(log(1), log(0.4)),
        start.nu = nu
      )
      qq <- INLA::inla.cgeneric.q(mod)
      expect_equal(nrow(qq$Q), m * mesh$n)
      expect_equal(
        max(abs(qq$Q - precision(mod))) / max(abs(precision(mod))), 0,
        tolerance = 1e-7
      )
    }
  }
})

test_that("the cgeneric table has one block per floor(alpha)", {
  tb <- rSPDE:::wl2_cgeneric_table(d = 2, m = 2, nu_upper_bound = 2)
  expect_equal(dim(tb), c(2 * 999, 2 * 2 + 2))
  expect_equal(attr(tb, "m_alpha_min"), 1)
  expect_equal(
    colnames(tb), c("alpha", "r1", "r2", "p1", "p2", "p0")
  )
  ## the rows are the two bands, in order
  expect_equal(unname(tb[c(1, 999, 1000), "alpha"]), c(1.001, 1.999, 2.001))
  ## and each row is a valid model
  expect_true(all(tb[, c("r1", "r2")] > 0))
  expect_true(all(tb[, c("p1", "p2")] < 1))
  expect_true(all(tb[, "p0"] >= 0 & tb[, "p0"] < 1))
  ## a row agrees with a direct interpolation
  cf <- rSPDE:::wl2_interp_coefficients(
    rSPDE:::wl2_coefficient_table(d = 2, m = 2, m_alpha = 2), 2.5
  )
  expect_equal(unname(tb[1000 + 499, -1]), c(cf$r, cf$p, cf$p0))
  ## one band only when the prior on nu cannot reach the next
  tb1 <- rSPDE:::wl2_cgeneric_table(d = 1, m = 1, nu_upper_bound = 0.4)
  expect_equal(nrow(tb1), 999)
  expect_equal(attr(tb1, "m_alpha_min"), 0)
})

test_that("the exact one-dimensional INLA models take the weighted-L2 type", {
  ## rspde.matern1d builds the exact Markov representation in C. With the
  ## weighted-L2 classes there is no constant term, so the field is
  ## rspde.order blocks of floor(alpha) + 1 entries per location rather than
  ## that plus a k-block, and the pole blocks carry the shared shift.
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_1d_general_wl2_model"
  )
  set.seed(1)
  loc <- sort(runif(25))
  sigma <- 1.2
  range <- 0.3
  for (m in 1:3) {
    for (nu in c(0.3, 0.9, 1.4)) {
      kappa <- sqrt(8 * nu) / range
      ## fixed nu: the index mapping is the identity, so Q is directly the
      ## matrix the R construction builds
      mod <- rspde.matern1d(
        loc = loc, rspde.order = m, nu = nu, parameterization = "matern",
        type.rational.approx = "wl2", shared_lib = "rSPDE",
        start.theta = c(log(sigma), log(range))
      )
      qq <- INLA::inla.cgeneric.q(mod)
      ref <- rSPDE:::matern.rational.precision(
        loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma,
        type_rational = "wl2"
      )$Q
      ## m blocks of floor(alpha) + 1 per location, and no k-block
      expect_equal(nrow(qq$Q), m * (floor(nu + 0.5) + 1) * length(loc))
      expect_equal(nrow(qq$Q), nrow(ref))
      expect_equal(max(abs(qq$Q - ref)) / max(abs(ref)), 0, tolerance = 1e-8)
    }
  }
})

test_that("the one-dimensional INLA model pads the graph for an estimated nu", {
  ## With nu estimated the graph is sized for nu.upper.bound and the entries
  ## the current nu does not use are left isolated with a unit diagonal. The
  ## active block must still be the model the R code builds.
  local_rspde_safe_inla(
    required_symbol = "inla_cgeneric_rspde_1d_general_wl2_model"
  )
  set.seed(1)
  loc <- sort(runif(25))
  sigma <- 1.2
  range <- 0.3
  m <- 2
  for (nu in c(0.3, 0.9, 1.4)) {
    kappa <- sqrt(8 * nu) / range
    mod <- rspde.matern1d(
      loc = loc, rspde.order = m, nu.upper.bound = 2,
      parameterization = "matern", type.rational.approx = "wl2",
      shared_lib = "rSPDE", start.theta = c(log(sigma), log(range)),
      start.nu = nu
    )
    Qc <- INLA::inla.cgeneric.q(mod)$Q
    ## the graph is sized for the upper bound throughout
    expect_equal(nrow(Qc), m * (floor(2 + 0.5) + 1) * length(loc))
    dg <- diag(Qc)
    offd <- (Matrix::rowSums(abs(Qc)) - abs(dg)) > 1e-12
    act <- which(offd | abs(dg - 1) > 1e-12)
    ref <- rSPDE:::matern.rational.precision(
      loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma,
      type_rational = "wl2"
    )$Q
    expect_equal(length(act), nrow(ref))
    expect_equal(
      max(abs(Qc[act, act] - ref)) / max(abs(ref)), 0, tolerance = 1e-8
    )
  }
})

test_that("the one-dimensional INLA interface refuses what it cannot do", {
  set.seed(1)
  loc <- sort(runif(20))
  ## an integer alpha needs no rational approximation
  expect_error(
    rspde.matern1d(
      loc = loc, nu = 0.5, rspde.order = 2, type.rational.approx = "wl2"
    ),
    "non-integer alpha"
  )
  ## floor(alpha) beyond 2 is out of scope for the classes
  expect_error(
    rspde.matern1d(
      loc = loc, nu.upper.bound = 3, rspde.order = 2,
      type.rational.approx = "wl2"
    ),
    "0, 1 or 2"
  )
})

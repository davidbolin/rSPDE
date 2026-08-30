#' Hybrid Whittle-Matern SPDE model for INLA / inlabru (alpha = 2)
#'
#' Creates an INLA cgeneric latent model for the hybrid Whittle-Matern
#' SPDE
#' \deqn{Y(s) = \beta_X^T L^{-1} X(s) + \tau^{-1} L^{-1} W(s),}
#' where \eqn{L = \kappa^2 - \Delta} and the deterministic mean
#' component \eqn{\beta_X^T L^{-1} X} is part of the latent model
#' itself: the cgeneric stores the centred field
#' \eqn{u(s) = \tau^{-1} L^{-1} W(s)} and returns
#' \eqn{\mu(s) = \beta_X^T L^{-1} X(s)} via the `mu` callback. The
#' coefficients `beta_X`, together with `tau` and `kappa`, are estimated
#' jointly as cgeneric hyperparameters -- this is the genuine fully
#' Bayesian fit, not a plug-in approximation.
#'
#' Currently the implementation is restricted to the case alpha = 2
#' (which is the main case of practical interest: in 2d this is the
#' Matern with `nu = 1`).
#'
#' @param mesh An `fmesher` mesh (1d or 2d).
#' @param X Numeric vector or matrix of covariate fields evaluated at
#' the mesh nodes; one column per covariate.
#' @param prior.tau A list with `mean` and `prec` (Gaussian prior on
#' `log(tau)`). Defaults: mean = 0, prec = 0.1.
#' @param prior.kappa A list with `mean` and `prec` (Gaussian prior on
#' `log(kappa)`). Defaults: mean derived from the mesh range, prec = 0.1.
#' @param prior.beta_x A list with `mean` (vector of length `p`) and
#' `prec` (vector of length `p`). Independent Gaussian priors on the
#' regression coefficients. Defaults: mean = 0, prec = 0.001
#' (uninformative).
#' @param start.ltau,start.lkappa Starting values for `log(tau)` and
#' `log(kappa)`. Default to the prior means.
#' @param start.beta_x Numeric vector of starting values for `beta_X`.
#' Default: zeros.
#' @param ... Additional arguments passed to the underlying model
#' constructor.
#' @param separate_kappa_mu Logical. If `TRUE`, the deterministic mean
#' operator uses its own range parameter `kappa_mu`, estimated as an
#' extra hyperparameter, rather than sharing `kappa` with the random
#' field. Default: `FALSE`.
#' @param prior.kappa_mu A list with `mean` and `prec` (Gaussian prior on
#' `log(kappa_mu)`). Only used when `separate_kappa_mu = TRUE`. Defaults
#' to the same prior as `prior.kappa`.
#' @param start.lkappa_mu Starting value for `log(kappa_mu)`. Only used
#' when `separate_kappa_mu = TRUE`. Defaults to the prior mean.
#' @param debug Passed to INLA.
#' @param shared_lib Which shared lib to use for the cgeneric. See
#' [rspde.matern()] for details.
#'
#' @return An INLA cgeneric latent model that can be used in INLA
#' formulae as `f(idx, model = <model>)`. 
#'
#' @details
#' The cgeneric hyperparameters are
#' \deqn{\theta = (\log\tau,\; \log\kappa,\; \beta_{X,1}, \ldots, \beta_{X,p}).}
#' The model stores the FEM `C` and `G` matrices and the covariate
#' matrix `X` and (re)computes
#' \eqn{Q = \tau^2\, L\, C^{-1}\, L} and
#' \eqn{\mu = L^{-1}(X \beta_X)} at every INLA call. The mean
#' computation is a single sparse Cholesky solve done with Eigen.
#'
#' @export
#' @seealso [hybrid.spde()], [rspde.matern()]
#' @examplesIf rspde_safe_inla(required_symbol = "inla_cgeneric_rspde_hybrid_alpha2_model")
#' library(INLA)
#' x <- seq(0, 1, length.out = 81)
#' mesh <- fmesher::fm_mesh_1d(x)
#' X <- cbind(sin(2 * pi * x), cos(2 * pi * x))
#' hyb <- rspde.hybrid.matern(mesh = mesh, X = X)
#' # ... use in formula: y ~ -1 + f(idx, model = hyb)
rspde.hybrid.matern <- function(mesh,
                                X,
                                prior.tau = NULL,
                                prior.kappa = NULL,
                                prior.beta_x = NULL,
                                start.ltau = NULL,
                                start.lkappa = NULL,
                                start.beta_x = NULL,
                                ...,
                                separate_kappa_mu = FALSE,
                                prior.kappa_mu = NULL,
                                start.lkappa_mu = NULL,
                                debug = FALSE,
                                shared_lib = "detect") {

  if (!inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d"))) {
    stop("'mesh' must be an fmesher mesh (fm_mesh_1d or fm_mesh_2d).")
  }
  d <- fmesher::fm_manifold_dim(mesh)

  X <- as.matrix(X)
  if (any(!is.finite(X))) {
    stop("'X' must contain only finite values.")
  }

  # FEM matrices: mass-lumped C, stiffness G.
  fem <- fm_fem(mesh)
  n_mesh <- nrow(fem$c0)
  if (nrow(X) != n_mesh) {
    stop(paste0("'X' must have nrow equal to the number of mesh nodes (",
                n_mesh, ")."))
  }
  p <- ncol(X)

  # Mass-lumped C (diagonal). We pass it as a numeric vector
  # (its diagonal) since the C++ aux helpers expect a dense
  # vector representation.
  C_diag <- as.numeric(Matrix::rowSums(fem$c0))
  G_mat  <- fem$g1

  # Convert G to a triplet sparse matrix for cgeneric.
  to_dgT <- function(A) {
    A <- as(A, "TsparseMatrix")
    methods::as(A, "dgTMatrix")
  }
  G_dgT <- to_dgT(G_mat)

  # Graph of Q = tau^2 L C^{-1} L. The non-zero pattern is the same as
  # L*L which equals the pattern of L + G_2 where G_2 = G C^{-1} G.
  # We build a symbolic Q using kappa=1 just to read out the sparsity.
  # We then apply the same `transpose_cgeneric` indexing convention as
  # rspde.matern: swap i and j (so that we list the matrix in
  # column-major form) and keep only the entries with i <= j after the
  # swap (i.e. the lower triangle of the original).
  C_lumped <- Matrix::Diagonal(n_mesh, x = C_diag)
  L_sym <- G_mat + C_lumped
  Cinv  <- Matrix::Diagonal(n_mesh, x = 1 / C_diag)
  Q_sym <- L_sym %*% Cinv %*% L_sym
  Q_sym <- Matrix::forceSymmetric(Q_sym)
  Q_graph <- transpose_cgeneric(Q_sym)
  graph_i <- Q_graph@i
  graph_j <- Q_graph@j
  M <- length(graph_i)

  # Default priors.
  if (is.null(prior.tau$mean))  prior.tau$mean  <- 0
  if (is.null(prior.tau$prec))  prior.tau$prec  <- 0.1
  if (is.null(prior.kappa$prec)) prior.kappa$prec <- 0.1
  if (is.null(prior.kappa$mean)) {
    # 20% of mesh range as a generic default for kappa.
    if (d == 1) {
      mesh_range <- diff(range(mesh$loc))
    } else {
      mesh_range <- max(diff(range(mesh$loc[, 1])),
                        diff(range(mesh$loc[, 2])))
    }
    prior.kappa$mean <- log(sqrt(8 * 0.5) / (0.2 * mesh_range))
  }
  if (is.null(prior.beta_x$mean)) prior.beta_x$mean <- rep(0, p)
  if (is.null(prior.beta_x$prec)) prior.beta_x$prec <- rep(0.001, p)
  if (length(prior.beta_x$mean) != p)
    stop("'prior.beta_x$mean' must have length ncol(X).")
  if (length(prior.beta_x$prec) != p)
    stop("'prior.beta_x$prec' must have length ncol(X).")

  # kappa_mu prior / start (only used when separate).
  separate_kappa_mu <- isTRUE(separate_kappa_mu)
  if (separate_kappa_mu) {
    if (is.null(prior.kappa_mu$mean)) prior.kappa_mu$mean <- prior.kappa$mean
    if (is.null(prior.kappa_mu$prec)) prior.kappa_mu$prec <- prior.kappa$prec
    if (is.null(start.lkappa_mu))     start.lkappa_mu     <- prior.kappa_mu$mean
  } else {
    # Provide dummy values so the cgeneric C asserts succeed; the C code
    # ignores them when separate_kappa_mu = 0.
    prior.kappa_mu       <- list(mean = 0, prec = 1)
    start.lkappa_mu      <- 0
  }

  # Starting values.
  if (is.null(start.ltau))   start.ltau   <- prior.tau$mean
  if (is.null(start.lkappa)) start.lkappa <- prior.kappa$mean
  if (is.null(start.beta_x)) start.beta_x <- rep(0, p)
  if (length(start.beta_x) != p)
    stop("'start.beta_x' must have length ncol(X).")

  if (separate_kappa_mu) {
    start_theta <- c(start.ltau, start.lkappa, start.lkappa_mu,
                     as.numeric(start.beta_x))
  } else {
    start_theta <- c(start.ltau, start.lkappa, as.numeric(start.beta_x))
  }
  theta.prior.mean <- c(prior.tau$mean, prior.kappa$mean)
  theta.prior.prec <- diag(c(prior.tau$prec, prior.kappa$prec))

  rspde_lib <- get_shared_library(shared_lib)

  model <- do.call(
    eval(parse(text = "INLA::inla.cgeneric.define")),
    list(
      model = "inla_cgeneric_rspde_hybrid_alpha2_model",
      shlib = rspde_lib,
      n     = as.integer(n_mesh),
      debug = debug,
      graph_opt_i = as.integer(graph_i),
      graph_opt_j = as.integer(graph_j),
      p           = as.integer(p),
      separate_kappa_mu = as.integer(separate_kappa_mu),
      C_diag      = as.double(C_diag),
      G_mat       = G_dgT,
      X_mat       = as.matrix(X),
      theta.prior.mean = as.double(theta.prior.mean),
      theta.prior.prec = as.matrix(theta.prior.prec),
      start.theta      = as.double(start_theta),
      beta_x.prior.mean = as.double(prior.beta_x$mean),
      beta_x.prior.prec = as.double(prior.beta_x$prec),
      kappa_mu.prior.mean = as.double(prior.kappa_mu$mean),
      kappa_mu.prior.prec = as.double(prior.kappa_mu$prec)
    )
  )
  model <- rspde_resolve_cgeneric_model(model)
  rspde_check_cgeneric_symbol(model)

  # Attributes used by the inlabru mapper for inla_rspde models so that
  # bru_get_mapper.inla_rspde returns a single-replicate mesh mapper.
  model$mesh        <- mesh
  model$rspde.order <- 0L
  model$integer.nu  <- TRUE
  model$est_nu      <- 0L
  model$n.spde      <- n_mesh
  model$alpha       <- 2
  model$nu          <- 2 - d / 2
  model$d           <- d
  model$p           <- p
  model$X           <- X
  model$separate_kappa_mu <- separate_kappa_mu
  model$cgeneric_type <- "hybrid_alpha2"
  class(model) <- c("inla_rspde_hybrid_alpha2", "inla_rspde",
                    class(model))

  return(model)
}


#' @export
#' @method print inla_rspde_hybrid_alpha2
print.inla_rspde_hybrid_alpha2 <- function(x, ...) {
  cat("Hybrid Whittle-Matern SPDE model for INLA / inlabru (alpha = 2)\n")
  cat("  Y = beta_X^T L_mu^{-1} X + tau^{-1} L^{-1} W\n")
  cat("  Number of covariate fields (p):", x$p, "\n")
  cat("  Number of mesh nodes:          ", x$n.spde, "\n")
  cat("  Dimension d:                   ", x$d, "\n")
  cat("  Smoothness nu = 2 - d/2:       ", format(x$nu), "\n")
  if (isTRUE(x$separate_kappa_mu)) {
    cat("  Mean operator: separate kappa_mu (its own hyperparameter)\n")
    cat("\nHyperparameters: log(tau), log(kappa), log(kappa_mu), ",
        "beta_x1, ..., beta_x", x$p, "\n", sep = "")
  } else {
    cat("  Mean operator: kappa_mu tied to kappa\n")
    cat("\nHyperparameters: log(tau), log(kappa), beta_x1, ..., beta_x",
        x$p, "\n", sep = "")
  }
  invisible(x)
}

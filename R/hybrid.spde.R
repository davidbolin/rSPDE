#' Hybrid non-stationary Whittle-Matern SPDE model with non-zero mean
#'
#' `hybrid.spde` constructs a finite element/rational SPDE approximation
#' of the model
#' \deqn{Y(s) = \beta L^{-\alpha/2} X(s) + \tau^{-1} L^{-\alpha/2} W(s),}
#' where \eqn{L = \kappa(s)^2 - \Delta} is the (possibly non-stationary)
#' Whittle-Matern differential operator, \eqn{\tau} is a positive scale
#' that affects only the random component, \eqn{X} is a (matrix of)
#' covariate field(s) supplied at the mesh nodes, \eqn{\beta} is a
#' vector of regression coefficients and \eqn{W} is Gaussian white
#' noise. 
#'
#' @param X Numeric vector or matrix containing the covariate field(s)
#' evaluated at the mesh nodes. If a matrix is supplied, each column is
#' treated as a separate covariate field. The number of rows must equal
#' the number of mesh nodes of the underlying FEM discretization.
#' @param beta_X Numeric vector of regression coefficients corresponding
#' to the columns of `X`. Defaults to a vector of zeros of length
#' `ncol(X)`.
#' @param kappa,tau,theta,B.tau,B.kappa,B.sigma,B.range,alpha,nu,parameterization,G,C,d,graph,mesh,range_mesh,loc_mesh,m,type,type_rational_approximation,check_stationarity
#' Arguments passed to [spde.matern.operators()]. The default
#' `check_stationarity = TRUE` allows the underlying covariance model
#' to be returned as a stationary [matern.operators()] object whenever
#' `tau` and `kappa` are constant after the log-linear regression on
#' `theta` is applied.
#'
#' @return `hybrid.spde` returns an object of class `"hybrid_spde"`
#' which also inherits from the class of the object returned by
#' [spde.matern.operators()] (`spde_matern_operator` plus either
#' `rSPDEobj` or `CBrSPDEobj` depending on `type`). The object stores
#' the covariate matrix `X`, the regression coefficients `beta_X`, the
#' cached operator object `op_mean` used to evaluate \eqn{L^{-\alpha/2}},
#' and the evaluated mean function `mu` at the mesh nodes.
#'
#' @seealso [spde.matern.operators()], [simulate.hybrid_spde()],
#' [predict.hybrid_spde()]
#' @export
#' @examples
#' set.seed(1)
#' x <- seq(from = 0, to = 1, length.out = 51)
#' kappa <- 10
#' tau <- 0.5
#' alpha <- 1.3
#' X <- cbind(sin(2 * pi * x), cos(2 * pi * x))
#' beta_X <- c(1.5, -0.5)
#' op <- hybrid.spde(
#'   X = X, beta_X = beta_X,
#'   kappa = kappa, tau = tau, alpha = alpha,
#'   loc_mesh = x, d = 1, parameterization = "spde"
#' )
#' Y <- simulate(op, nsim = 1)
hybrid.spde <- function(X = NULL,
                        beta_X = NULL,
                        kappa = NULL,
                        tau = NULL,
                        theta = NULL,
                        B.tau = NULL,
                        B.kappa = NULL,
                        B.sigma = NULL,
                        B.range = NULL,
                        alpha = NULL,
                        nu = NULL,
                        parameterization = c("spde", "matern"),
                        G = NULL,
                        C = NULL,
                        d = NULL,
                        graph = NULL,
                        mesh = NULL,
                        range_mesh = NULL,
                        loc_mesh = NULL,
                        m = 1,
                        type = c("covariance", "operator"),
                        type_rational_approximation = c(
                          "brasil",
                          "chebfun",
                          "chebfunLB"
                        ),
                        check_stationarity = TRUE) {

  if (is.null(X)) {
    stop("'X' must be supplied as a vector or matrix of covariate values at the mesh nodes.")
  }
  X <- as.matrix(X)
  if (any(!is.finite(X))) {
    stop("'X' must contain only finite values.")
  }
  if (is.null(beta_X)) {
    beta_X <- rep(0, ncol(X))
  }
  beta_X <- as.numeric(beta_X)
  if (length(beta_X) != ncol(X)) {
    stop("'beta_X' must have length equal to ncol(X).")
  }

  parameterization <- parameterization[[1]]
  type <- type[[1]]

  # Only forward the B.* matrices that match the chosen parameterization
  # so that spde.matern.operators' missing()-based input checks behave
  # the same way as when called directly by the user.
  spde_args <- list(
    kappa = kappa, tau = tau, theta = theta,
    alpha = alpha, nu = nu,
    parameterization = parameterization,
    G = G, C = C, d = d,
    graph = graph, mesh = mesh,
    range_mesh = range_mesh, loc_mesh = loc_mesh,
    m = m,
    type = type,
    type_rational_approximation = type_rational_approximation,
    check_stationarity = check_stationarity
  )
  if (parameterization == "spde") {
    if (!is.null(B.tau))   spde_args$B.tau <- B.tau
    if (!is.null(B.kappa)) spde_args$B.kappa <- B.kappa
  } else {
    if (!is.null(B.sigma)) spde_args$B.sigma <- B.sigma
    if (!is.null(B.range)) spde_args$B.range <- B.range
  }
  op_cov <- do.call(spde.matern.operators, spde_args)

  n_mesh <- dim(op_cov$C)[1]
  if (nrow(X) != n_mesh) {
    stop(paste0(
      "The number of rows of 'X' (", nrow(X),
      ") must equal the number of mesh nodes (", n_mesh, ")."
    ))
  }

  # spde.matern.operators with type = "operator" does not explicitly
  # store G in the output, but we need it to be able to rebuild op_mean
  # when the model parameters are updated. Recover G from the original
  # FEM matrices when possible.
  if (is.null(op_cov$G)) {
    if (!is.null(mesh)) {
      fem <- fm_fem(mesh)
      op_cov$G <- fem$g1
    } else if (!is.null(graph)) {
      op_cov$G <- graph$mesh$G
    } else if (!is.null(loc_mesh) && (is.null(d) || d == 1)) {
      mesh_1d <- fm_mesh_1d(loc_mesh)
      fem <- fm_fem(mesh_1d)
      op_cov$G <- fem$g1
    } else if (!is.null(G)) {
      op_cov$G <- G
    } else {
      stop("Could not infer the stiffness matrix G; please supply it explicitly.")
    }
  }

  op_cov$X <- X
  op_cov$beta_X <- beta_X

  op_cov$op_mean <- build_op_mean_hybrid(
    C = op_cov$C, G = op_cov$G,
    kappa = op_cov$kappa, alpha = op_cov$alpha, m = op_cov$m
  )

  op_cov$mu <- compute_hybrid_mean(op_cov)

  class(op_cov) <- c("hybrid_spde", class(op_cov))
  return(op_cov)
}


#' @name update.hybrid_spde
#' @title Update parameters of hybrid_spde objects
#' @description Updates the model parameters of a `hybrid_spde` object.
#' Delegates the update of the covariance part to the underlying
#' [spde.matern.operators()] update method and additionally rebuilds the
#' operator used to evaluate the mean and recomputes the stored mean
#' vector.
#' @param object A `hybrid_spde` object.
#' @param X If non-null, update the covariate matrix.
#' @param beta_X If non-null, update the regression coefficients.
#' @param ... Additional arguments passed to the underlying update
#' method of [spde.matern.operators()].
#' @return An object of class `"hybrid_spde"` with the updated
#' parameters.
#' @method update hybrid_spde
#' @seealso [hybrid.spde()]
#' @export
update.hybrid_spde <- function(object, X = NULL, beta_X = NULL, ...) {
  # Stash hybrid-specific fields before delegating to the underlying
  # update method.
  X_store <- if (is.null(X)) object$X else as.matrix(X)
  beta_store <- if (is.null(beta_X)) object$beta_X else as.numeric(beta_X)
  if (length(beta_store) != ncol(X_store)) {
    stop("'beta_X' must have length equal to ncol(X).")
  }
  G_store <- object$G

  class(object) <- setdiff(class(object), "hybrid_spde")
  object$X <- NULL
  object$beta_X <- NULL
  object$op_mean <- NULL
  object$mu <- NULL

  new_object <- update(object, ...)

  if (is.null(new_object$G)) {
    new_object$G <- G_store
  }
  new_object$X <- X_store
  new_object$beta_X <- beta_store
  new_object$op_mean <- build_op_mean_hybrid(
    C = new_object$C, G = new_object$G,
    kappa = new_object$kappa,
    alpha = new_object$alpha,
    m = new_object$m
  )
  new_object$mu <- compute_hybrid_mean(new_object)

  class(new_object) <- c("hybrid_spde", class(new_object))
  return(new_object)
}


#' Simulation of a hybrid Whittle-Matern SPDE model
#'
#' Samples the hybrid SPDE model
#' \deqn{L^{\alpha/2}(\tau Y(s)) = \beta X(s) + W(s)}
#' by first sampling the mean-zero part using the simulation method of
#' the underlying [spde.matern.operators()] approximation and then
#' adding the deterministic mean
#' \eqn{\mu = \tau^{-1} \beta L^{-\alpha/2} X}.
#'
#' @param object A `hybrid_spde` object returned by [hybrid.spde()].
#' @param nsim Number of samples to generate.
#' @param seed Optional integer used to initialise the random number
#' generator.
#' @param ... Additional arguments passed to the underlying simulate
#' method (e.g. updated parameter values).
#' @return A matrix with `nsim` columns and one row per mesh node.
#' @method simulate hybrid_spde
#' @seealso [hybrid.spde()]
#' @export
simulate.hybrid_spde <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "hybrid_spde")
  U <- simulate(obj_base, nsim = nsim, ...)
  U <- as.matrix(U)
  mu <- compute_hybrid_mean(object)
  if (length(mu) != nrow(U)) {
    stop("Internal error: the dimension of the mean does not match the simulated field.")
  }
  U <- U + matrix(mu, nrow = nrow(U), ncol = ncol(U))
  return(U)
}


#' Kriging prediction for a hybrid Whittle-Matern SPDE model
#'
#' Computes posterior predictions for the hybrid model
#' \eqn{Y_i = (A Y)_i + \epsilon_i}, where \eqn{Y} satisfies
#' \eqn{L^{\alpha/2}(\tau Y) = \beta X + W}. The posterior of the
#' centred latent field \eqn{Y - \mu} is computed using the underlying
#' SPDE prediction routine; the deterministic mean \eqn{\mu} is then
#' added back at the prediction locations.
#'
#' @param object A `hybrid_spde` object.
#' @param A Projection matrix linking the observations to the mesh nodes.
#' @param Aprd Projection matrix linking the prediction locations to the
#' mesh nodes.
#' @param Y Observations.
#' @param sigma.e Standard deviation of the measurement noise.
#' @param compute.variances If TRUE, the kriging variances are returned.
#' @param posterior_samples If TRUE, posterior samples are returned.
#' @param n_samples Number of posterior samples.
#' @param only_latent If TRUE, the posterior samples are not perturbed
#' by the measurement noise.
#' @param ... Additional arguments passed to the underlying predict
#' method.
#' @return A list with components `mean`, optionally `variance` and
#' `samples`.
#' @method predict hybrid_spde
#' @seealso [hybrid.spde()]
#' @export
predict.hybrid_spde <- function(object, A, Aprd, Y, sigma.e,
                                compute.variances = FALSE,
                                posterior_samples = FALSE,
                                n_samples = 100,
                                only_latent = FALSE,
                                ...) {
  Y <- as.matrix(Y)
  mu <- compute_hybrid_mean(object)
  mu_obs <- as.vector(A %*% mu)
  mu_prd <- as.vector(Aprd %*% mu)
  Y_centered <- Y - matrix(mu_obs, nrow = nrow(Y), ncol = ncol(Y))

  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "hybrid_spde")
  out <- predict(obj_base, A = A, Aprd = Aprd, Y = Y_centered,
                 sigma.e = sigma.e,
                 compute.variances = compute.variances,
                 posterior_samples = posterior_samples,
                 n_samples = n_samples,
                 only_latent = only_latent,
                 ...)
  out$mean <- as.matrix(out$mean) +
    matrix(mu_prd, nrow = length(mu_prd), ncol = ncol(as.matrix(out$mean)))
  if (!is.null(out$samples)) {
    if (is.list(out$samples)) {
      out$samples <- lapply(out$samples, function(s) {
        s + matrix(mu_prd, nrow = nrow(s), ncol = ncol(s))
      })
    } else {
      out$samples <- out$samples +
        matrix(mu_prd, nrow = nrow(out$samples), ncol = ncol(out$samples))
    }
  }
  return(out)
}


#' @export
#' @method make_A hybrid_spde
make_A.hybrid_spde <- function(object, loc, ...) {
  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "hybrid_spde")
  return(make_A(obj_base, loc = loc, ...))
}


#' @export
#' @method cov_function_mesh hybrid_spde
cov_function_mesh.hybrid_spde <- function(object, p, ...) {
  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "hybrid_spde")
  return(cov_function_mesh(obj_base, p = p, ...))
}


#' @export
#' @method covariance_mesh hybrid_spde
covariance_mesh.hybrid_spde <- function(object, ...) {
  obj_base <- object
  class(obj_base) <- setdiff(class(obj_base), "hybrid_spde")
  return(covariance_mesh(obj_base, ...))
}


#' @export
#' @method print hybrid_spde
print.hybrid_spde <- function(x, ...) {
  cat("Hybrid Whittle-Matern SPDE model with non-zero mean\n")
  cat("  Y = beta_X L^{-alpha/2} X + tau^{-1} L^{-alpha/2} W\n")
  cat("  Number of mesh nodes:        ", dim(x$C)[1], "\n", sep = "")
  cat("  Number of covariate fields:  ", ncol(as.matrix(x$X)), "\n", sep = "")
  cat("  beta_X:                      ",
      paste(format(x$beta_X), collapse = ", "), "\n", sep = "")
  cat("  alpha:                       ", format(x$alpha), "\n", sep = "")
  cat("  Rational approximation order:", x$m, "\n")
  if (inherits(x, "CBrSPDEobj")) {
    cat("  Covariance part: covariance-based rational approximation\n")
  } else {
    cat("  Covariance part: operator-based rational approximation\n")
  }
  invisible(x)
}

#' @noRd
#'
#' Build an internal fractional.operators object with tau = 1 that is
#' used to evaluate L^{-alpha/2} X for the deterministic mean of the
#' hybrid SPDE model.
build_op_mean_hybrid <- function(C, G, kappa, alpha, m) {
    n <- dim(C)[1]
    if (length(kappa) == 1) {
        L <- G + C * kappa^2
    } else {
        if (length(kappa) != n) {
            stop("kappa must have length 1 or equal to the number of mesh nodes")
        }
        kp <- Matrix::Diagonal(n, kappa^2)
        L <- G + C %*% kp
    }
    fractional.operators(
        L = L,
        beta = alpha / 2,
        C = C,
        scale.factor = min(kappa)^2,
        m = m,
        tau = 1
    )
}

#' @noRd
#'
#' Compute the deterministic mean beta L^{-alpha/2} X(s) of a hybrid SPDE model 
#' at the mesh nodes. 
#' Uses the (cached) `op_mean` field of the object if it is available,
#' otherwise rebuilds it from the model state.
compute_hybrid_mean <- function(object) {
    if (is.null(object$X) || is.null(object$beta_X)) {
        return(rep(0, dim(object$C)[1]))
    }
    X <- as.matrix(object$X)
    if (all(object$beta_X == 0)) {
        return(rep(0, nrow(X)))
    }
    op_mean <- object$op_mean
    if (is.null(op_mean)) {
        op_mean <- build_op_mean_hybrid(
            C = object$C, G = object$G,
            kappa = object$kappa, alpha = object$alpha, m = object$m
        )
    }
    # Convert X (function point values at the mesh nodes) to the FEM
    # weak-form representation by multiplying with the mass-lumped
    # diagonal of C, then apply the operator-based approximation of
    # L^{-alpha/2}.
    C_diag <- Matrix::diag(op_mean$C)
    CX <- C_diag * X
    LinvX <- as.matrix(Pr.mult(op_mean, Pl.solve(op_mean, CX)))
    mu <- as.vector(LinvX %*% object$beta_X)
    return(mu)
}

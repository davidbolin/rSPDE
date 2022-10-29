###########################################
## fractional.computations.R
##
##   Copyright (C) 2018, 2019, David Bolin
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Simulation of a fractional SPDE using a rational SPDE approximation
#'
#' @description The function samples a Gaussian random field based on a
#' pre-computed rational SPDE approximation.
#'
#' @param object The rational SPDE approximation, computed
#' using [fractional.operators()],
#' [matern.operators()], or [spde.matern.operators()].
#' @param nsim The number of simulations.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param ... Currently not used.
#'
#' @return A matrix with the `n` samples as columns.
#' @seealso [simulate.CBrSPDEobj()]
#' @export
#' @method simulate rSPDEobj
#'
#' @examples
#' # Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation
#' op <- matern.operators(
#'   kappa = kappa, sigma = sigma,
#'   nu = nu, G = fem$G, C = fem$C, d = 1
#' )
#'
#' # Sample the model and plot the result
#' Y <- simulate(op)
#' plot(x, Y, type = "l", ylab = "u(x)", xlab = "x")
#'
simulate.rSPDEobj <- function(object,
                              nsim = 1,
                              seed = NULL,
                              ...) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if (!inherits(object, "rSPDEobj")) {
    stop("input object is not of class rSPDEobj")
  }
  m <- dim(object$Q)[1]
  z <- rnorm(nsim * m)
  dim(z) <- c(m, nsim)
  x <- Qsqrt.solve(object, z)
  x <- Pr.mult(object, x)

  return(x)
}


#' @name update.CBrSPDEobj
#' @title Update parameters of CBrSPDEobj objects
#' @description Function to change the parameters of a CBrSPDEobj object
#' @param object The covariance-based rational SPDE approximation,
#' computed using [matern.operators()]
#' @param user_kappa If non-null, update the range parameter
#' of the covariance function.
#' @param user_sigma If non-null, update the standard deviation of
#' the covariance function.
#' @param user_nu If non-null, update the shape parameter of the
#' covariance function.
#' @param user_m If non-null, update the order of the rational
#' approximation, which needs to be a positive integer.
#' @param compute_higher_order Logical. Should the higher order
#' finite element matrices be computed?
#' @param return_block_list Logical. For `type = "covariance"`,
#' should the block parts of the precision matrix be returned
#' separately as a list?
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are "chebfun",
#' "brasil" or "chebfunLB".
#' @param ... Currently not used.
#' @return It returns an object of class "CBrSPDEobj. This object contains the
#' same quantities listed in the output of [matern.operators()].
#' @method update CBrSPDEobj
#' @seealso [simulate.CBrSPDEobj()], [matern.operators()]
#' @export
#' @examples
#' # Compute the covariance-based rational approximation of a
#' # Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#' op_cov
#'
#' # Update the range parameter of the model:
#' op_cov <- update(op_cov, user_kappa = 20)
#' op_cov
#'
update.CBrSPDEobj <- function(object, user_nu = NULL,
                              user_kappa = NULL,
                              user_sigma = NULL,
                              user_m = NULL,
                              compute_higher_order = object$higher_order,
                              type_rational_approximation =
                              object$type_rational_approximation,
                              return_block_list = object$return_block_list,
                              ...) {
  new_object <- object
  d <- object$d

  fem_mesh_matrices <- object$fem_mesh_matrices

  if (is.null(user_nu) && !(object$higher_order) && compute_higher_order) {
    user_nu <- object$nu
  }

  ## get parameters
  if (!is.null(user_nu)) {
    new_object$nu <- rspde_check_user_input(user_nu, "nu")
    nu <- user_nu
    alpha <- nu + d / 2
    m_alpha <- floor(alpha)
    m_order <- m_alpha + 1

    if (compute_higher_order) {
      if (m_order + 1 > length(object$fem_mesh_matrices)) {
        old_m_order <- length(object$fem_mesh_matrices) - 1
        GCi <- object$GCi
        for (i in (old_m_order + 1):m_order) {
          fem_mesh_matrices[[paste0("g", i)]] <- GCi %*%
          fem_mesh_matrices[[paste0("g", i - 1)]]
        }
      }
    }
  }

  new_object[["fem_mesh_matrices"]] <- fem_mesh_matrices

  if (!is.null(user_kappa)) {
    new_object$kappa <- rspde_check_user_input(user_kappa, "kappa")
  }

  if (!is.null(user_sigma)) {
    new_object$sigma <- rspde_check_user_input(user_sigma, "sigma")
  }

  if (!is.null(user_m)) {
    new_object$m <- as.integer(rspde_check_user_input(user_m, "m", 1))
  }

  new_object <- CBrSPDE.matern.operators(
    kappa = new_object$kappa,
    sigma = new_object$sigma,
    nu = new_object$nu,
    G = new_object$G,
    C = new_object$C,
    d = new_object$d,
    m = new_object$m,
    return_block_list = return_block_list,
    type_rational_approximation = type_rational_approximation,
    mesh = NULL,
    fem_mesh_matrices = new_object$fem_mesh_matrices
  )

  return(new_object)
}


#' @name update.rSPDEobj
#' @title Update parameters of rSPDEobj objects
#' @description Function to change the parameters of a rSPDEobj object
#' @param object The operator-based rational SPDE approximation,
#' computed using [matern.operators()] with `type="operator"`
#' @param user_kappa If non-null, update the range parameter
#' of the covariance function.
#' @param user_sigma If non-null, update the standard
#' deviation of the covariance function.
#' @param user_nu If non-null, update the shape parameter
#' of the covariance function.
#' @param user_m If non-null, update the order of the rational
#' approximation, which needs to be a positive integer.
#' @param ... Currently not used.
#' @return It returns an object of class "rSPDEobj. This object contains the
#' same quantities listed in the output of [matern.operators()].
#' @method update rSPDEobj
#' @seealso [simulate.rSPDEobj()], [matern.operators()]
#' @export
#' @examples
#' # Compute the operator-based rational approximation of a
#' # Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' op <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2, type = "operator"
#' )
#' op
#'
#' # Update the range parameter of the model:
#' op <- update(op, user_kappa = 20)
#' op
#'
update.rSPDEobj <- function(object, user_nu = NULL,
                            user_kappa = NULL,
                            user_sigma = NULL,
                            user_m = NULL, ...) {
  new_object <- object

  ## get parameters
  if (!is.null(user_nu)) {
    new_object$nu <- rspde_check_user_input(user_nu, "nu")
  }

  if (!is.null(user_kappa)) {
    new_object$kappa <- rspde_check_user_input(user_kappa, "kappa")
  }

  if (!is.null(user_sigma)) {
    new_object$sigma <- rspde_check_user_input(user_sigma, "sigma")
  }

  if (!is.null(user_m)) {
    new_object$m <- as.integer(rspde_check_user_input(user_m, "m", 1))
  }

  new_object <- matern.operators(
    kappa = new_object$kappa,
    sigma = new_object$sigma,
    nu = new_object$nu,
    G = new_object$G,
    C = new_object$C,
    d = new_object$d,
    m = new_object$m,
    type = "operator"
  )

  return(new_object)
}


#' @name simulate.CBrSPDEobj
#' @title Simulation of a fractional SPDE using the
#' covariance-based rational SPDE approximation
#' @description The function samples a Gaussian random field based using the
#' covariance-based rational SPDE approximation.
#' @param object The covariance-based rational SPDE approximation,
#' computed using [matern.operators()]
#' @param nsim The number of simulations.
#' @param seed An object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param user_kappa If non-null, update the range parameter
#' of the covariance function.
#' @param user_sigma If non-null, update the standard deviation
#' of the covariance function.
#' @param user_nu If non-null, update the shape parameter of the
#' covariance function.
#' @param user_m If non-null, update the order of the rational
#' approximation, which needs to be a positive integer.
#' @param pivot Should pivoting be used for the Cholesky
#' decompositions? Default is TRUE
#' @param ... Currently not used.
#' @return A matrix with the `n` samples as columns.
#' @method simulate CBrSPDEobj
#' @export
#' @examples
#' # Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#' (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model and plot the result
#' Y <- simulate(op_cov)
#' plot(x, Y, type = "l", ylab = "u(x)", xlab = "x")
#'
simulate.CBrSPDEobj <- function(object, nsim = 1,
                                seed = NULL,
                                user_nu = NULL,
                                user_kappa = NULL,
                                user_sigma = NULL,
                                user_m = NULL,
                                pivot = TRUE,
                                ...) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  d <- object$d
  nu_temp <- ifelse(is.null(user_nu), object$nu, user_nu)
  alpha <- nu_temp + d / 2


  ## simulation
  if (alpha %% 1 == 0) { # simulation in integer case
    object <- update.CBrSPDEobj(
      object = object,
      user_nu = user_nu,
      user_kappa = user_kappa,
      user_sigma = user_sigma,
      user_m = user_m,
      compute_higher_order = TRUE
    )


    kappa <- object$kappa
    nu <- object$nu
    tau <- object$tau
    m <- object$m

    alpha <- nu + d / 2

    fem_mesh_matrices <- object$fem_mesh_matrices

    L <- fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]] / kappa^2
    sizeL <- dim(L)[1]


    Q <- rspde.matern.precision.integer(
      kappa = kappa, nu = nu, tau = tau,
      dim = d,
      fem_mesh_matrices = fem_mesh_matrices
    )
    Z <- rnorm(sizeL * nsim)
    dim(Z) <- c(sizeL, nsim)
    if (pivot) {
      LQ <- chol(forceSymmetric(Q), pivot = TRUE)
      reorder <- attr(LQ, "pivot")
      X <- solve(LQ, Z)
      # order back
      orderback <- numeric(length(reorder))
      orderback[reorder] <- seq_len(length(reorder))
      X <- X[orderback, ]
    } else {
      LQ <- chol(forceSymmetric(Q))
      X <- solve(LQ, Z)
    }
  } else {
    object <- update.CBrSPDEobj(
      object = object,
      user_nu = user_nu,
      user_kappa = user_kappa,
      user_sigma = user_sigma,
      user_m = user_m,
      compute_higher_order = FALSE
    )

    m <- object$m
    Q <- object$Q

    Z <- rnorm(dim(Q)[1] * nsim)
    dim(Z) <- c(dim(Q)[1], nsim)
    if (pivot) {
      LQ <- chol(forceSymmetric(Q), pivot = TRUE)
      reorder <- attr(LQ, "pivot")
      X <- solve(LQ, Z)
      # order back
      orderback <- numeric(length(reorder))
      orderback[reorder] <- seq_len(length(reorder))
      X <- X[orderback, ]
    } else {
      LQ <- chol(forceSymmetric(Q))
      X <- solve(LQ, Z)
    }
    A <- Diagonal(dim(Q)[1] / (m + 1))
    Abar <- kronecker(matrix(1, 1, m + 1), A)
    X <- Abar %*% X
  }
  return(X)
}


#' Prediction of a fractional SPDE using a rational SPDE approximation
#'
#' The function is used for computing kriging predictions based on data
#' \eqn{Y_i = u(s_i) + \epsilon_i},
#' where \eqn{\epsilon}{\epsilon} is mean-zero Gaussian measurement noise
#' and \eqn{u(s)}{u(s)} is defined by
#' a fractional SPDE \eqn{L^\beta u(s) = W}{L^\beta u(s) = W}, where
#' \eqn{W}{W} is Gaussian white noise.
#'
#' @param object The rational SPDE approximation, computed using
#' [fractional.operators()],
#' [matern.operators()], or [spde.matern.operators()].
#' @param A A matrix linking the measurement locations to the basis of the
#' FEM approximation of the latent model.
#' @param Aprd A matrix linking the prediction locations to the basis of the
#' FEM approximation of the latent model.
#' @param Y A vector with the observed data, can also be a matrix where the
#' columns are observations
#' of independent replicates of \eqn{u}.
#' @param sigma.e The standard deviation of the Gaussian measurement noise.
#' Put to zero if the model
#' does not have measurement noise.
#' @param compute.variances Set to also TRUE to compute the kriging variances.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list with elements
#' \item{mean }{The kriging predictor (the posterior mean of u|Y).}
#' \item{variance }{The posterior variances (if computed).}
#' @export
#' @method predict rSPDEobj
#'
#' @examples
#' # Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#' sigma.e <- 0.3
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation
#' op <- matern.operators(
#'   kappa = kappa, sigma = sigma,
#'   nu = nu, G = fem$G, C = fem$C, d = 1
#' )
#'
#' # Sample the model
#' u <- simulate(op)
#'
#' # Create some data
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y <- as.vector(A %*% u + sigma.e * rnorm(10))
#'
#' # compute kriging predictions at the FEM grid
#' A.krig <- rSPDE.A1d(x, x)
#' u.krig <- predict(op,
#'   A = A, Aprd = A.krig, Y = Y, sigma.e = sigma.e,
#'   compute.variances = TRUE
#' )
#'
#' plot(obs.loc, Y,
#'   ylab = "u(x)", xlab = "x", main = "Data and prediction",
#'   ylim = c(
#'     min(u.krig$mean - 2 * sqrt(u.krig$variance)),
#'     max(u.krig$mean + 2 * sqrt(u.krig$variance))
#'   )
#' )
#' lines(x, u.krig$mean)
#' lines(x, u.krig$mean + 2 * sqrt(u.krig$variance), col = 2)
#' lines(x, u.krig$mean - 2 * sqrt(u.krig$variance), col = 2)
predict.rSPDEobj <- function(object,
                             A,
                             Aprd,
                             Y,
                             sigma.e,
                             compute.variances = FALSE,
                             ...) {
  Y <- as.matrix(Y)
  if (dim(Y)[1] != dim(A)[1]) {
    stop("the dimensions of A does not match the number of observations")
  }

  out <- list()
  if (length(sigma.e) == 1) {
    if (sigma.e < 0) {
      stop("sigma.e must be non-negative")
    } else if (sigma.e > 0) {
      A <- A %*% object$Pr
      AA <- Aprd %*% object$Pr
      Qhat <- object$Q + (t(A) %*% A) / sigma.e^2
      out$mean <- as.matrix(AA %*% solve(Qhat, t(A) %*% Y / sigma.e^2))
      if (compute.variances) {
        out$variance <- diag(AA %*% solve(Qhat, t(AA)))
      }
    } else { # no nugget
      Ahat <- A %*% object$Pr
      QiAt <- solve(object$Q, t(Ahat))
      AQiA <- Ahat %*% QiAt
      xhat <- solve(object$Q, t(Ahat) %*% solve(AQiA, Y))
      out$mean <- as.vector(Aprd %*% xhat)
      if (compute.variances) {
        AA <- Aprd %*% object$Pr
        M <- object$Q - QiAt %*% solve(AQiA, t(QiAt))
        out$variance <- diag(AA %*% M %*% t(AA))
      }
    }
  } else if (dim(Y)[1] == length(sigma.e)) {
    Q.e <- Diagonal(length(sigma.e), 1 / sigma.e^2)
    A <- A %*% object$Pr
    AA <- Aprd %*% object$Pr
    Qhat <- object$Q + t(A) %*% Q.e %*% A
    out$mean <- as.matrix(AA %*% solve(Qhat, t(A) %*% Q.e %*% Y))
    if (compute.variances) {
      out$variance <- diag(AA %*% solve(Qhat, t(AA)))
    }
  }
  return(out)
}

#' Object-based log-likelihood function for latent Gaussian
#' fractional SPDE model
#'
#' This function evaluates the log-likelihood function for a
#' fractional SPDE model
#' \eqn{L^\beta u(s) = W}{L^\beta u(s) = W} that is observed under
#' Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = x(s_i) + \epsilon_i},
#' where \eqn{\epsilon_i}{\epsilon_i}
#' are iid mean-zero Gaussian variables and \eqn{x(s) =
#' \mu(s) + u(s)}{x(s) = \mu(s) + u(s)}, where
#' \eqn{\mu(s)}{\mu(s)} is the expectation vector of the latent field.
#'
#' @param obj The rational SPDE approximation, computed using
#' [fractional.operators()],
#' [matern.operators()], or [spde.matern.operators()].
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location
#' to the finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @return The log-likelihood value.
#' @export
#' @note This example below shows how the function can be used to evaluate
#' the likelihood of a latent Matern model. Se [matern.loglike()]
#' for an example of how this can be used for maximum
#' likelihood estimation.
#' @seealso [matern.loglike()], [spde.matern.loglike()]
#'
#' @examples
#' # Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#' sigma.e <- 0.3
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation
#' op <- matern.operators(
#'   kappa = kappa, sigma = sigma, nu = nu,
#'   G = fem$G, C = fem$C, d = 1,
#'   type = "operator"
#' )
#'
#' # Sample the model
#' u <- simulate(op)
#'
#' # Create some data
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y <- as.vector(A %*% u + sigma.e * rnorm(10))
#'
#' # compute log-likelihood of the data
#' lik1 <- rSPDE.loglike(op, Y, A, sigma.e)
#' cat(lik1)
rSPDE.loglike <- function(obj,
                          Y,
                          A,
                          sigma.e,
                          mu = 0) {
  Y <- as.matrix(Y)
  if (length(dim(Y)) == 2) {
    n.rep <- dim(Y)[2]
    n <- dim(Y)[1]
  } else {
    n.rep <- 1
    if (length(dim(Y)) == 1) {
      n <- dim(Y)[1]
    } else {
      n <- length(Y)
    }
  }
  if (length(sigma.e) == 1) {
    Q.e <- Diagonal(n) / sigma.e^2
    nugget <- rep(sigma.e^2, n)
  } else {
    if (length(sigma.e) != n) {
      stop("the length of sigma.e does not match the number of observations")
    }
    Q.e <- Diagonal(length(sigma.e), 1 / sigma.e^2)
    nugget <- sigma.e^2
  }



  R <- Matrix::Cholesky(obj$Pl)

  prior.ld <- 4 * c(determinant(R, logarithm = TRUE)$modulus) -
  sum(log(diag(obj$C)))


  A <- A %*% obj$Pr
  Q.post <- obj$Q + t(A) %*% Q.e %*% A


  R.post <- Matrix::Cholesky(Q.post)

  posterior.ld <- 2 * c(determinant(R.post, logarithm = TRUE)$modulus)

  AtY <- t(A) %*% Q.e %*% Y


  mu.post <- mu + solve(R.post, AtY, system = "A")


  lik <- n.rep * (prior.ld - posterior.ld - dim(A)[1] *
  log(2 * pi) - sum(log(nugget))) / 2

  if (n.rep > 1) {
    lik <- lik - 0.5 * sum(colSums((mu.post - mu) * (obj$Q %*% (mu.post - mu))))
    v <- Q.e %*% (Y - A %*% mu.post)
    lik <- lik - 0.5 * sum(colSums((Y - A %*% mu.post) * v))
  } else {
    lik <- lik - 0.5 * (t(mu.post - mu) %*% obj$Q %*% (mu.post - mu) +
    t(Y - A %*% mu.post) %*% Q.e %*% (Y - A %*% mu.post))
  }
  return(as.double(lik))
}


#' @name rSPDE.matern.loglike
#' @title Object-based log-likelihood function for latent Gaussian fractional
#' SPDE model using the rational approximations
#' @description This function evaluates the log-likelihood function for a
#' Gaussian process with a Matern covariance
#' function, that is observed under Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are
#' iid mean-zero Gaussian variables. The latent model is approximated using
#' the a rational approximation
#' of the fractional SPDE model corresponding to the Gaussian process.
#' @param object The rational SPDE approximation,
#' computed using [matern.operators()]
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location to the
#' finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param user_kappa If non-null, update the range parameter of the covariance
#' function.
#' @param user_sigma If non-null, update the standard deviation of the
#' covariance function.
#' @param user_nu If non-null, update the shape parameter of the covariance
#' function.
#' @param user_m If non-null, update the order of the rational approximation,
#' which needs to be a positive integer.
#' @param pivot Should pivoting be used for the Cholesky decompositions? Default
#' is TRUE
#' @return The log-likelihood value.
#' @export
#' @seealso [matern.operators()], [predict.CBrSPDEobj()]
#' @examples
#' # this example illustrates how the function can be used for maximum
#' # likelihood estimation
#'
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' nu <- 0.8
#' kappa <- 5
#' sigma <- 1
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#' (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#'
#' # Compute the covariance-based rational approximation
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model
#' u <- simulate(op_cov, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # Define the negative likelihood function for optimization
#' # using CBrSPDE.matern.loglike
#'
#' # Notice that we are also using sigma instead of tau, so it can be compared
#' # to matern.loglike()
#' mlik_cov <- function(theta, Y, A, op_cov) {
#'   kappa <- exp(theta[1])
#'   sigma <- exp(theta[2])
#'   nu <- exp(theta[3])
#'   return(-rSPDE.matern.loglike(
#'     object = op_cov, Y = Y,
#'     A = A, user_kappa = kappa, user_sigma = sigma,
#'     user_nu = nu, sigma.e = exp(theta[4])
#'   ))
#' }
#'
#' # The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), 1 / sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, mlik_cov,
#'   Y = Y, A = A, op_cov = op_cov,
#'   method = "L-BFGS-B"
#' )
#'
#' print(data.frame(
#'   kappa = c(kappa, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
#'
rSPDE.matern.loglike <- function(object, Y, A, sigma.e, mu = 0,
                                 user_nu = NULL,
                                 user_kappa = NULL,
                                 user_sigma = NULL,
                                 user_m = NULL,
                                 pivot = TRUE) {
  if (inherits(object, "CBrSPDEobj")) {
    return(CBrSPDE.matern.loglike(
      object = object,
      Y = Y, A = A,
      sigma.e = sigma.e,
      mu = mu,
      user_nu = user_nu,
      user_kappa = user_kappa,
      user_sigma = user_sigma,
      user_m = user_m,
      pivot = pivot
    ))
  } else {
    if (inherits(object, "rSPDEobj")) {
      if (object$type == "Matern approximation") {
        object <- update.rSPDEobj(object,
          user_nu = user_nu,
          user_kappa = user_kappa,
          user_sigma = user_sigma,
          user_m = user_m
        )
        return(rSPDE.loglike(obj = object, Y = Y, A = A,
        sigma.e = sigma.e, mu = mu))
      } else {
        stop("The fractional operator should be of type
        'Matern approximation'!")
      }
    } else {
      stop("The object should be of class 'CBrSPDEobj' or 'rSPDEobj'!")
    }
  }
}

#' @name CBrSPDE.matern.loglike
#' @title Object-based log-likelihood function for latent Gaussian fractional
#' SPDE model using the covariance-based rational approximations
#' @description This function evaluates the log-likelihood function for a
#' Gaussian process with a Matern covariance function, that is observed under
#' Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables.
#' The latent model is approximated using
#' the covariance-based rational approximation
#' of the fractional SPDE model corresponding to the Gaussian process.
#' @param object The covariance-based rational SPDE approximation,
#' computed using [matern.operators()]
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location
#' to the finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param user_kappa If non-null, update the range parameter of the
#' covariance function.
#' @param user_sigma If non-null, update the standard deviation of the
#' covariance function.
#' @param user_nu If non-null, update the shape parameter of the
#' covariance function.
#' @param user_m If non-null, update the order of the rational approximation,
#' which needs to be a positive integer.
#' @param pivot Should pivoting be used for the Cholesky decompositions?
#' Default is TRUE
#' @return The log-likelihood value.
#' @noRd
#' @seealso [matern.operators()], [predict.CBrSPDEobj()]
#' @examples
#' # this example illustrates how the function can be used for maximum
#' likelihood estimation
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' nu <- 0.8
#' kappa <- 5
#' sigma <- 1
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi)^(1 / 2) *
#' gamma(nu + 1 / 2)))
#'
#' # Compute the covariance-based rational approximation
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model
#' u <- simulate(op_cov, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # Define the negative likelihood function for optimization using
#' # CBrSPDE.matern.loglike
#' # Notice that we are also using sigma instead of tau, so it can be compared
#' # to matern.loglike()
#' mlik_cov <- function(theta, Y, A, op_cov) {
#'   kappa <- exp(theta[1])
#'   sigma <- exp(theta[2])
#'   nu <- exp(theta[3])
#'   return(-rSPDE.matern.loglike(
#'     object = op_cov, Y = Y,
#'     A = A, user_kappa = kappa, user_sigma = sigma,
#'     user_nu = nu, sigma.e = exp(theta[4])
#'   ))
#' }
#'
#' # The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), 1 / sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, mlik_cov,
#'   Y = Y, A = A, op_cov = op_cov,
#'   method = "L-BFGS-B"
#' )
#'
#' print(data.frame(
#'   kappa = c(kappa, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
CBrSPDE.matern.loglike <- function(object, Y, A, sigma.e, mu = 0,
                                   user_nu = NULL,
                                   user_kappa = NULL,
                                   user_sigma = NULL,
                                   user_m = NULL,
                                   pivot = TRUE) {
  Y <- as.matrix(Y)
  if (length(dim(Y)) == 2) {
    n.rep <- dim(Y)[2]
    n <- dim(Y)[1]
  } else {
    n.rep <- 1
    if (length(dim(Y)) == 1) {
      n <- dim(Y)[1]
    } else {
      n <- length(Y)
    }
  }

  ## get relevant parameters

  object <- update.CBrSPDEobj(
    object = object,
    user_nu = user_nu,
    user_kappa = user_kappa,
    user_sigma = user_sigma,
    user_m = user_m
  )

  m <- object$m

  if (length(sigma.e) == 1) {
    Q.e <- Diagonal(n) / sigma.e^2
    nugget <- rep(sigma.e^2, n)
  } else {
    if (length(sigma.e) != n) {
      stop("the length of sigma.e does not match the number of observations")
    }
    Q.e <- Diagonal(length(sigma.e), 1 / sigma.e^2)
    nugget <- sigma.e^2
  }

  Q.frac <- object$Q.frac

  Q.fracR <- chol(Q.frac, pivot = pivot)

  logdetL <- object$logdetL
  logdetC <- object$logdetC
  Q.int.order <- object$Q.int$order

  if (Q.int.order > 0) {
    logQ <- 2 * sum(log(diag(Q.fracR))) + (Q.int.order) *
    (m + 1) * (logdetL - logdetC)
  } else {
    logQ <- 2 * sum(log(diag(Q.fracR)))
  }

  ## compute Q_x|y
  Q <- object$Q
  Abar <- kronecker(matrix(1, 1, m + 1), A)
  Q_xgiveny <- t(Abar) %*% Q.e %*% Abar + Q
  ## construct mu_x|y

  mu_xgiveny <- t(Abar) %*% Q.e %*% Y
  # upper triangle with reordering


  R <- Matrix::Cholesky(Q_xgiveny)



  mu_xgiveny <- solve(R, mu_xgiveny, system = "A")

  mu_xgiveny <- mu + mu_xgiveny

  ## compute log|Q_xgiveny|
  log_Q_xgiveny <- 2 * determinant(R, logarithm = TRUE)$modulus
  ## compute mu_x|y*Q*mu_x|y
  if (n.rep > 1) {
    mu_part <- sum(colSums((mu_xgiveny - mu) * (Q %*% (mu_xgiveny - mu))))
  } else {
    mu_part <- t(mu_xgiveny - mu) %*% Q %*% (mu_xgiveny - mu)
  }
  ## compute central part
  if (n.rep > 1) {
    central_part <- sum(colSums((Y - Abar %*% mu_xgiveny) * (Q.e %*%
    (Y - Abar %*% mu_xgiveny))))
  } else {
    central_part <- t(Y - Abar %*% mu_xgiveny) %*% Q.e %*% (Y -
    Abar %*% mu_xgiveny)
  }
  ## compute log|Q_epsilon|
  log_Q_epsilon <- -sum(log(nugget))
  ## wrap up
  log_likelihood <- n.rep * (logQ + log_Q_epsilon - log_Q_xgiveny) -
  mu_part - central_part
  if (n.rep > 1) {
    log_likelihood <- log_likelihood - dim(A)[1] * n.rep * log(2 * pi)
  } else {
    log_likelihood <- log_likelihood - length(Y) * log(2 * pi)
  }
  log_likelihood <- log_likelihood / 2

  return(as.double(log_likelihood))
}



#' Parameter-based log-likelihood for a latent Gaussian Matern model
#' using a rational SPDE approximation
#'
#' This function evaluates the log-likelihood function for a Gaussian
#' process with a Matern covariance
#' function, that is observed under Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables.
#' The latent model is approximated using a rational approximation
#' of the fractional SPDE model corresponding to the Gaussian process.
#'
#' @param kappa Range parameter of the latent process.
#' @param sigma Standard deviation of the latent process.
#' @param nu Shape parameter of the latent process.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param G The stiffness matrix of a finite element discretization
#' of the domain.
#' @param C The mass matrix of a finite element discretization of the domain.
#' @param A A matrix linking the measurement locations to the basis of
#' the FEM approximation of the latent model.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a
#' positive integer. The default value is 1.
#' @param type The type of the rational approximation. The options are
#' "covariance" and "operator". The default is "covariance".
#' @param pivot Should pivoting be used for the Cholesky decompositions?
#' Default is TRUE
#'
#' @return The log-likelihood value.
#' @export
#' @seealso [spde.matern.loglike()], [rSPDE.loglike()],
#' [matern.operators()].
#'
#' @examples
#' # this example illustrates how the function can be used for maximum
#' # likelihood estimation
#'
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using the covariance-based
#' # rational approximation
#' nu <- 0.8
#' kappa <- 5
#' sigma <- 1
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi)^(1 / 2) *
#' gamma(nu + 1 / 2)))
#'
#' # Compute the covariance-based rational approximation
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model
#' u <- simulate(op_cov, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # Define the negative likelihood function for optimization
#' # using CBrSPDE.matern.loglike
#' # Notice that we are also using sigma instead of tau, so it can be compared
#' # to matern.loglike()
#' mlik_cov2 <- function(theta, Y, A, C, G) {
#'   kappa <- exp(theta[1])
#'   sigma <- exp(theta[2])
#'   nu <- exp(theta[3])
#'   return(-matern.loglike(
#'     kappa = kappa, sigma = sigma,
#'     nu = nu, sigma.e = exp(theta[4]), Y = Y, A = A,
#'     C = fem$C, G = fem$G, d = 1
#'   ))
#' }
#'
#' # The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, mlik_cov2,
#'   Y = Y, A = A, C = C, G = G,
#'   method = "L-BFGS-B"
#' )
#'
#' print(data.frame(
#'   kappa = c(kappa, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
#'
#' # this example illustrates how the function can be used for
#' # maximum likelihood estimation when using the operator-based
#' # rational approximation
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' nu <- 0.8
#' kappa <- 5
#' sigma <- 1
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation
#' op <- matern.operators(
#'   kappa = kappa, sigma = sigma, nu = nu,
#'   G = fem$G, C = fem$C, d = 1,
#'   type = "operator"
#' )
#'
#' # Sample the model
#' u <- simulate(op, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # define negative likelihood function for optimization using matern.loglike
#' mlik <- function(theta, Y, G, C, A) {
#'   return(-matern.loglike(exp(theta[1]), exp(theta[2]),
#'     exp(theta[3]), exp(theta[4]),
#'     Y = Y, G = G, C = C, A = A, d = 1,
#'     type = "operator"
#'   ))
#' }
#'
#' # The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, mlik,
#'   Y = Y, G = fem$G, C = fem$C, A = A,
#'   method = "L-BFGS-B"
#' )
#'
#' print(data.frame(
#'   kappa = c(kappa, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
matern.loglike <- function(kappa,
                           sigma,
                           nu,
                           sigma.e,
                           Y,
                           G,
                           C,
                           A,
                           mu = 0,
                           d = 2,
                           m = 1,
                           type = c("covariance", "operator"),
                           pivot = TRUE) {
  type <- type[[1]]
  if (!type %in% c("covariance", "operator")) {
    stop("The type should be 'covariance' or 'operator'!")
  }
  if (is.null(d)) {
    stop("the dimension d must be supplied")
  }
  if (type == "covariance") {
    return(CBrSPDE.matern.loglike2(
      kappa = kappa,
      sigma = sigma,
      nu = nu,
      sigma.e = sigma.e,
      mu = mu,
      Y = Y,
      G = G,
      C = C,
      A = A,
      d = d,
      m = m,
      pivot = pivot
    ))
  } else {
    op <- matern.operators(
      kappa = kappa, sigma = sigma, nu = nu,
      G = G, C = C, d = d, m = m,
      type = "operator"
    )
    return(rSPDE.loglike(obj = op, Y = Y, A = A, sigma.e = sigma.e))
  }
}


#' @name CBrSPDE.matern.loglike2
#' @title Parameter-based log-likelihood function for latent Gaussian fractional
#' SPDE model using the covariance-based rational approximations
#' @description This function evaluates the log-likelihood function for a
#' Gaussian process with a Matern covariance function, that is observed under
#' Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables.
#' The latent model is approximated using the covariance-based rational
#' approximation of the fractional SPDE model corresponding to the
#' Gaussian process.
#' @param kappa Range parameter of the latent process.
#' @param sigma Standard deviation of the latent process.
#' @param nu Shape parameter of the latent process.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param G The stiffness matrix of a finite element discretization
#' of the domain.
#' @param C The mass matrix of a finite element discretization of the domain.
#' @param A A matrix linking the measurement locations to the basis of
#' the FEM approximation of the latent model.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be
#' a positive integer.
#' The default value is 2.
#' @param pivot Should pivoting be used for the Cholesky decompositions?
#' Default is TRUE
#' @return The log-likelihood value.
#' @noRd
#' @seealso [matern.operators()], [predict.CBrSPDEobj()]
#' @examples
#' # this example illustrates how the function can be used for maximum
#' # likelihood estimation
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' nu <- 0.8
#' kappa <- 5
#' sigma <- 1
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi)^(1 / 2) *
#' gamma(nu + 1 / 2)))
#'
#' # Compute the covariance-based rational approximation
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model
#' u <- simulate(op_cov, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # Define the negative likelihood function for optimization using
#' # CBrSPDE.matern.loglike2
#' # Notice that we are also using sigma instead of tau, so it can be compared
#' # to matern.loglike()
#' mlik_cov2 <- function(theta, Y, A, C, G) {
#'   kappa <- exp(theta[1])
#'   sigma <- exp(theta[2])
#'   nu <- exp(theta[3])
#'   return(-matern.loglike(
#'     kappa = kappa, sigma = sigma,
#'     nu = nu, sigma.e = exp(theta[4]), Y = Y,
#'     A = A, C = fem$C, G = fem$G, d = 1
#'   ))
#' }
#'
#' # The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, mlik_cov2,
#'   Y = Y, A = A, C = C, G = G,
#'   method = "L-BFGS-B"
#' )
#'
#' print(data.frame(
#'   kappa = c(kappa, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
CBrSPDE.matern.loglike2 <- function(kappa,
                                    sigma,
                                    nu,
                                    sigma.e,
                                    mu = 0,
                                    Y,
                                    G,
                                    C,
                                    A,
                                    d = 2,
                                    m = 2,
                                    pivot = TRUE) {
  obj_cov_rSPDE <- matern.operators(
    C = C,
    G = G,
    nu = nu,
    kappa = kappa,
    sigma = sigma,
    m = m,
    d = d
  )

  return(CBrSPDE.matern.loglike(
    object = obj_cov_rSPDE, Y = Y, A = A, sigma.e = sigma.e, mu = mu,
    user_nu = NULL,
    user_kappa = NULL,
    user_sigma = NULL,
    user_m = NULL,
    pivot = pivot
  ))
}






#' Parameter-based log-likelihood for a latent Gaussian Matern SPDE model
#' using a rational SPDE approximation
#'
#' This function evaluates the log-likelihood function for observations
#' of a Gaussian process defined as the solution to the SPDE
#' \deqn{(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) = W.}
#'
#' The observations are assumed to be generated as
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are
#' iid mean-zero Gaussian variables. The latent model is approximated using a
#' rational approximation of the fractional SPDE model.
#'
#' @param kappa Vector with the, possibly spatially varying, range parameter
#' evaluated at the locations of the mesh used for the finite element
#' discretization of the SPDE.
#' @param tau Vector with the, possibly spatially varying, precision
#' parameter evaluated at the locations of the mesh used for the
#' finite element discretization of the SPDE.
#' @param nu Shape parameter of the covariance function, related
#' to \eqn{\beta} through the equation \eqn{\beta = (\nu + d/2)/2}.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param G The stiffness matrix of a finite element discretization of
#' the domain.
#' @param C The mass matrix of a finite element discretization of the domain.
#' @param A A matrix linking the measurement locations to the basis of the
#' FEM approximation of the latent model.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a
#' positive integer. The default value is 1.
#'
#' @return The log-likelihood value.
#' @export
#' @seealso [matern.loglike()], [rSPDE.loglike()].
#'
#' @examples
#' # this example illustrates how the function can be used for maximum
#' # likelihood estimation
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- rep(0.5, n.x)
#' nu <- 0.8
#' kappa <- rep(1, n.x)
#'
#' # compute rational approximation
#' op <- spde.matern.operators(
#'   kappa = kappa, tau = tau, nu = nu,
#'   G = fem$G, C = fem$C, d = 1
#' )
#'
#' # Sample the model
#' u <- simulate(op, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # define negative likelihood function for optimization using matern.loglike
#' mlik <- function(theta, Y, G, C, A) {
#'   return(-spde.matern.loglike(rep(exp(theta[1]), n.x),
#'     rep(exp(theta[2]), n.x),
#'     exp(theta[3]), exp(theta[4]),
#'     Y = Y, G = G, C = C, A = A, d = 1
#'   ))
#' }
#'
#' #' #The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), 1 / sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, mlik, Y = Y, G = fem$G, C = fem$C, A = A)
#'
#' print(data.frame(
#'   kappa = c(kappa[1], exp(theta$par[1])), tau = c(tau[1], exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
spde.matern.loglike <- function(kappa,
                                tau,
                                nu,
                                sigma.e,
                                Y,
                                G,
                                C,
                                A,
                                d = 2,
                                m = 1) {
  op <- spde.matern.operators(kappa, tau, nu, G, C, d = d, m = m)
  return(rSPDE.loglike(op, Y, A, sigma.e))
}



#' @name predict.CBrSPDEobj
#' @title Prediction of a fractional SPDE using the covariance-based
#' rational SPDE approximation
#' @description The function is used for computing kriging predictions based
#' on data \eqn{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon}{\epsilon}
#' is mean-zero Gaussian measurement noise and \eqn{u(s)}{u(s)} is defined by
#' a fractional SPDE \eqn{(\kappa^2 I - \Delta)^{\alpha/2} (\tau u(s)) = W},
#' where \eqn{W}{W} is Gaussian white noise and \eqn{\alpha = \nu + d/2},
#' where \eqn{d} is the dimension of the domain.
#' @param object The covariance-based rational SPDE approximation,
#' computed using [matern.operators()]
#' @param A A matrix linking the measurement locations to the basis of the FEM
#' approximation of the latent model.
#' @param Aprd A matrix linking the prediction locations to the basis of the
#' FEM approximation of the latent model.
#' @param Y A vector with the observed data, can also be a matrix where the
#' columns are observations
#' of independent replicates of \eqn{u}.
#' @param sigma.e The standard deviation of the Gaussian measurement noise.
#' Put to zero if the model does not have measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param compute.variances Set to also TRUE to compute the kriging variances.
#' @param pivot Should pivoting be used on the Cholesky decompositions?
#' @param ... further arguments passed to or from other methods.
#' @return A list with elements
#' \item{mean }{The kriging predictor (the posterior mean of u|Y).}
#' \item{variance }{The posterior variances (if computed).}
#' @export
#' @method predict CBrSPDEobj
#' @examples
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#' sigma.e <- 0.3
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#'        (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#'
#' # Compute the covariance-based rational approximation
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model
#' u <- simulate(op_cov)
#'
#' # Create some data
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y <- as.vector(A %*% u + sigma.e * rnorm(10))
#'
#' # compute kriging predictions at the FEM grid
#' A.krig <- rSPDE.A1d(x, x)
#' u.krig <- predict(op_cov,
#'   A = A, Aprd = A.krig, Y = Y, sigma.e = sigma.e,
#'   compute.variances = TRUE
#' )
#'
#' plot(obs.loc, Y,
#'   ylab = "u(x)", xlab = "x", main = "Data and prediction",
#'   ylim = c(
#'     min(u.krig$mean - 2 * sqrt(u.krig$variance)),
#'     max(u.krig$mean + 2 * sqrt(u.krig$variance))
#'   )
#' )
#' lines(x, u.krig$mean)
#' lines(x, u.krig$mean + 2 * sqrt(u.krig$variance), col = 2)
#' lines(x, u.krig$mean - 2 * sqrt(u.krig$variance), col = 2)
predict.CBrSPDEobj <- function(object, A, Aprd, Y, sigma.e, mu = 0,
                               compute.variances = FALSE, pivot = TRUE,
                               ...) {
  Y <- as.matrix(Y)
  if (dim(Y)[1] != dim(A)[1]) {
    stop("the dimensions of A does not match the number of observations")
  }

  n <- dim(Y)[1]
  out <- list()

  d <- object$d
  kappa <- object$kappa
  nu <- object$nu
  tau <- object$tau
  m <- object$m

  fem_mesh_matrices <- object$fem_mesh_matrices

  no_nugget <- FALSE

  if (length(sigma.e) == 1) {
    if (sigma.e == 0) {
      no_nugget <- TRUE
    } else {
      Q.e <- Diagonal(n) / sigma.e^2
    }
  } else {
    if (length(sigma.e) != n) {
      stop("the length of sigma.e does not match the number of observations")
    }
    Q.e <- Diagonal(length(sigma.e), 1 / sigma.e^2)
  }

  alpha <- nu + d / 2

  if (!no_nugget) {
    if (alpha %% 1 == 0) { # loglikelihood in integer case
      ## construct Q
      Q <- rspde.matern.precision.integer(
        kappa = kappa, nu = nu,
        tau = tau,
        dim = d,
        fem_mesh_matrices = fem_mesh_matrices
      )

      R <- chol(forceSymmetric(Q), pivot = pivot)
      ## compute Q_x|y
      Q_xgiveny <- (t(A) %*% Q.e %*% A) + Q
      ## construct mu_x|y
      mu_xgiveny <- t(A) %*% Q.e %*% Y
      # upper triangle with reordering
      R <- chol(forceSymmetric(Q_xgiveny), pivot = pivot)
      if (pivot) {
        reorder <- attr(R, "pivot")
        # make it lower triangle
        R <- t(R)
        v <- solve(R, mu_xgiveny[reorder])
        mu_xgiveny <- solve(t(R), v)
        # order back
        orderback <- numeric(length(reorder))
        orderback[reorder] <- seq_len(length(reorder))
        mu_xgiveny <- mu_xgiveny[orderback]
      } else {
        R <- t(R)
        v <- solve(R, mu_xgiveny)
        mu_xgiveny <- solve(t(R), v)
      }

      mu_xgiveny <- mu + mu_xgiveny
      out$mean <- Aprd %*% mu_xgiveny

      if (compute.variances) {
        out$variance <- diag(Aprd %*% solve(Q_xgiveny, t(Aprd)))
      }
    } else { # loglikelihood in non-integer case

      Q <- object$Q

      R <- chol(forceSymmetric(Q), pivot = pivot)
      ## compute Q_x|y
      Q_xgiveny <- kronecker(matrix(1, m + 1, m + 1), t(A) %*% Q.e %*% A) + Q
      ## construct mu_x|y
      Abar <- kronecker(matrix(1, 1, m + 1), A)
      mu_xgiveny <- t(Abar) %*% Q.e %*% Y
      # upper triangle with reordering
      R <- chol(forceSymmetric(Q_xgiveny), pivot = pivot)
      if (pivot) {
        reorder <- attr(R, "pivot")
        # make it lower triangle
        R <- t(R)
        v <- solve(R, mu_xgiveny[reorder])
        mu_xgiveny <- solve(t(R), v)

        orderback <- numeric(length(reorder))
        orderback[reorder] <- seq_len(length(reorder))
        mu_xgiveny <- mu_xgiveny[orderback]
      } else {
        R <- t(R)
        v <- solve(R, mu_xgiveny)
        mu_xgiveny <- solve(t(R), v)
      }
      mu_xgiveny <- mu + mu_xgiveny

      Aprd_bar <- kronecker(matrix(1, 1, m + 1), Aprd)

      out$mean <- Aprd_bar %*% mu_xgiveny

      if (compute.variances) {
        out$variance <- diag(Aprd_bar %*% solve(Q_xgiveny, t(Aprd_bar)))
      }
    }
  } else {
    integer_alpha <- (alpha %% 1 == 0)
    if (integer_alpha) {
      Abar <- A
      Aprd_bar <- Aprd
      Q <- rspde.matern.precision.integer(
        kappa = kappa, nu = nu, tau = tau, dim = d,
        fem_mesh_matrices = fem_mesh_matrices
      )
    } else {
      Abar <- kronecker(matrix(1, 1, m + 1), A)
      Aprd_bar <- kronecker(matrix(1, 1, m + 1), Aprd)
      Q <- object$Q
    }

    QiAt <- solve(Q, t(Abar))
    AQiA <- Abar %*% QiAt
    xhat <- solve(Q, t(Abar) %*% solve(AQiA, Y))

    out$mean <- as.vector(Aprd_bar %*% xhat)
    if (compute.variances) {
      M <- Q - QiAt %*% solve(AQiA, t(QiAt))
      out$variance <- diag(Aprd_bar %*% M %*% t(Aprd_bar))
    }
  }


  return(out)
}



#' @rdname precision.CBrSPDEobj
#' @export
precision <- function(object, ...) {
  UseMethod("precision", object)
}

#' @name precision.CBrSPDEobj
#' @title Get the precision matrix of CBrSPDEobj objects
#' @description Function to get the precision matrix of a CBrSPDEobj object
#' @param object The covariance-based rational SPDE approximation,
#' computed using [matern.operators()]
#' @param user_kappa If non-null, update the range parameter of
#' the covariance function.
#' @param user_sigma If non-null, update the standard deviation of
#' the covariance function.
#' @param user_nu If non-null, update the shape parameter of the
#' covariance function.
#' @param user_m If non-null, update the order of the rational approximation,
#' which needs to be a positive integer.
#' @param ... Currently not used.
#' @return The precision matrix.
#' @method precision CBrSPDEobj
#' @seealso [simulate.CBrSPDEobj()], [matern.operators()]
#' @export
#' @examples
#' # Compute the covariance-based rational approximation of a
#' # Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#' (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Get the precision matrix:
#' prec_matrix <- precision(op_cov)
#'
precision.CBrSPDEobj <- function(object,
                                 user_nu = NULL,
                                 user_kappa = NULL,
                                 user_sigma = NULL,
                                 user_m = NULL,
                                 ...) {
  object <- update.CBrSPDEobj(
    object = object,
    user_nu = user_nu,
    user_kappa = user_kappa,
    user_sigma = user_sigma,
    user_m = user_m
  )

  Q <- object$Q
  return(Q)
}


#' @name rSPDE.construct.matern.loglike
#' @title Constructor of Matern loglikelihood functions.
#' @description This function returns a log-likelihood function for a
#' Gaussian process with a Matern covariance
#' function, that is observed under Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are
#' iid mean-zero Gaussian variables. The latent model is approximated using
#' the a rational approximation
#' of the fractional SPDE model corresponding to the Gaussian process.
#' @param object The rational SPDE approximation,
#' computed using [matern.operators()]
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location to the
#' finite element basis.
#' @param sigma.e IF non-null, the standard deviation of the measurement noise will be kept fixed in 
#' the returned likelihood.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param user_kappa If non-null, the range parameter will be kept fixed in the returned likelihood. 
#' @param user_sigma If non-null, the standard deviation will be kept fixed in the returned likelihood.
#' @param user_nu If non-null, the shape parameter will be kept fixed in the returned likelihood.
#' @param user_m If non-null, update the order of the rational approximation,
#' which needs to be a positive integer.
#' @param log_scale Should the parameters be evaluated in log-scale?
#' @param return_negative_likelihood Return minus the likelihood to turn the maximization into a minimization?
#' @param pivot Should pivoting be used for the Cholesky decompositions? Default
#' is TRUE
#' @return The log-likelihood function. The parameters of the returned function
#' are given in the order sigma, kappa, nu, sigma.e, whenever they are available.
#' @export
#' @seealso [matern.operators()], [predict.CBrSPDEobj()]
#' @examples
#' # this example illustrates how the function can be used for maximum
#' # likelihood estimation
#'
#' set.seed(123)
#' # Sample a Gaussian Matern process on R using a rational approximation
#' nu <- 0.8
#' kappa <- 5
#' sigma <- 1
#' sigma.e <- 0.1
#' n.rep <- 10
#' n.obs <- 100
#' n.x <- 51
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#' (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#'
#' # Compute the covariance-based rational approximation
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' # Sample the model
#' u <- simulate(op_cov, n.rep)
#'
#' # Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs * n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y <- as.matrix(A %*% u + sigma.e * noise)
#'
#' # Define the negative likelihood function for optimization
#' # using CBrSPDE.matern.loglike
#'
#' # Notice that we are also using sigma instead of tau, so it can be compared
#' # to matern.loglike()
#' loglike <- rSPDE.construct.matern.loglike(op_cov, Y, A) 
#' 
#' # The parameters can now be estimated by minimizing mlik with optim
#' \donttest{
#' # Choose some reasonable starting values depending on the size of the domain
#' theta0 <- log(c(sqrt(8), 1 / sqrt(var(c(Y))), 0.9, 0.01))
#'
#' # run estimation and display the results
#' theta <- optim(theta0, loglike,
#'   method = "L-BFGS-B"
#' )
#'
#' print(data.frame(
#'   kappa = c(kappa, exp(theta$par[1])), sigma = c(sigma, exp(theta$par[2])),
#'   nu = c(nu, exp(theta$par[3])), sigma.e = c(sigma.e, exp(theta$par[4])),
#'   row.names = c("Truth", "Estimates")
#' ))
#' }
#'
rSPDE.construct.matern.loglike <- function(object, Y, A, 
                                 sigma.e = NULL, mu = 0,
                                 user_nu = NULL,
                                 user_kappa = NULL,
                                 user_sigma = NULL,
                                 user_m = NULL,
                                 log_scale = TRUE,
                                 return_negative_likelihood = TRUE,
                                 pivot = TRUE){
        param_vector <- likelihood_process_inputs(user_kappa, user_sigma, user_nu, sigma.e)
        
        loglik <- function(theta){
          if(is.null(user_sigma)){
          sigma <- likelihood_process_parameters(theta = theta, 
                  param_vector = param_vector, 
                  which_par = "sigma", 
                  logscale = log_scale)
          } else{
            sigma <- user_sigma
          }
          if(is.null(user_kappa)){
          kappa <- likelihood_process_parameters(theta = theta, 
                  param_vector = param_vector, 
                  which_par = "kappa", 
                  logscale = log_scale)
          } else{
            kappa <- user_kappa
          }
          if(is.null(user_nu)){
          nu <- likelihood_process_parameters(theta = theta, 
                  param_vector = param_vector, 
                  which_par = "nu", 
                  logscale = log_scale)
         
          } else{
            nu <- user_nu
          }
          if(is.null(sigma.e)){
          sigma.e <- likelihood_process_parameters(theta = theta, 
                  param_vector = param_vector, 
                  which_par = "sigma.e", 
                  logscale = log_scale)
          } 
          loglike <- rSPDE.matern.loglike(object = object, Y=Y, A=A,
          sigma.e = sigma.e,
          mu = mu,
          user_kappa = kappa,
          user_nu = nu,
          user_sigma=sigma,
          user_m = user_m,
          pivot = pivot)
          if(return_negative_likelihood){
            return(-loglike)
          } else{
            return(loglike)
          }
        }
        return(loglik)
}


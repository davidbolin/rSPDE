## fractional.operators.R
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

#' Rational approximations of fractional operators
#'
#' `fractional.operators` is used for computing an approximation,
#' which can be used for inference and simulation, of the fractional SPDE
#' \deqn{L^\beta (\tau u(s)) = W.}
#' Here \eqn{L} is a differential operator, \eqn{\beta>0} is
#' the fractional power, \eqn{\tau} is a positive scalar or vector that
#' scales the variance of the solution \eqn{u}, and \eqn{W} is white noise.
#'
#' @param L A finite element discretization of the operator \eqn{L}{L}.
#' @param beta The positive fractional power.
#' @param C The mass matrix of the finite element discretization.
#' @param scale.factor A constant \eqn{c} is a lower bound for the the smallest
#' eigenvalue of the non-discretized operator \eqn{L}{L}.
#' @param m The order of the rational approximation, which needs to be a
#' positive integer. The default value is 1.
#' Higer values gives a more accurate approximation, which are more
#' computationally expensive to use for inference. Currently, the largest value
#' of m that is implemented is 4.
#' @param tau The constant or vector that scales the variance of the solution.
#' The default value is 1.
#'
#' @return `fractional.operators` returns an object of class "rSPDEobj".
#' This object contains the following quantities:
#' \item{Pl}{The operator \eqn{P_l}.}
#' \item{Pr}{The operator \eqn{P_r}.}
#' \item{C}{The mass lumped mass matrix.}
#' \item{Ci}{The inverse of `C`.}
#' \item{m}{The order of the rational approximation.}
#' \item{beta}{The fractional power.}
#' \item{type}{String indicating the type of approximation.}
#' \item{Q}{The matrix `t(Pl) %*% solve(C,Pl)`.}
#' \item{type}{String indicating the type of approximation.}
#' \item{Pl.factors}{List with elements that can be used to assemble \eqn{P_l}.}
#' \item{Pr.factors}{List with elements that can be used to assemble \eqn{P_r}.}
#' @export
#' @seealso [matern.operators()], [spde.matern.operators()],
#' [matern.operators()]
#' @details The approximation is based on a rational approximation of
#' the fractional operator, resulting in an
#' approximate model on the form \deqn{P_l u(s) = P_r W,}
#' where \eqn{P_j = p_j(L)} are non-fractional operators defined in terms of
#' polynomials \eqn{p_j} for \eqn{j=l,r}. The order of \eqn{p_r} is given by
#' `m` and the order of \eqn{p_l} is \eqn{m + m_\beta}
#' where \eqn{m_\beta} is the integer part of \eqn{\beta} if \eqn{\beta>1} and
#' \eqn{m_\beta = 1} otherwise.
#'
#' The discrete approximation can be written as \eqn{u = P_r x} where
#' \eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})} and \eqn{Q = P_l^T C^{-1} P_l}.
#' Note that the matrices \eqn{P_r} and \eqn{Q} may be be ill-conditioned
#' for \eqn{m>1}. In this case, the methods in [operator.operations()]
#' should be used for operations involving the matrices, since these methods
#' are more numerically stable.
#'
#' @examples
#' # Compute rational approximation of a Gaussian process with a
#' # Matern covariance function on R
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
#'   (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#' op <- fractional.operators(
#'   L = fem$G + kappa^2 * fem$C, beta = (nu + 1 / 2) / 2,
#'   C = fem$C, scale.factor = kappa^2, tau = tau
#' )
#'
#' v <- t(rSPDE.A1d(x, 0.5))
#' c.approx <- Sigma.mult(op, v)
#'
#' # plot the result and compare with the true Matern covariance
#' plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
#'   type = "l", ylab = "C(h)",
#'   xlab = "h", main = "Matern covariance and rational approximations"
#' )
#' lines(x, c.approx, col = 2)
#'
fractional.operators <- function(L,
                                 beta,
                                 C,
                                 scale.factor,
                                 m = 1,
                                 tau = 1) {
  if (min(tau) < 0) {
    stop("tau should be positive")
  }
  if ((m %% 1) != 0 || m < 0) {
    stop("m must be a positive integer")
  }
  if (scale.factor <= 0) {
    stop("the scaling factor must be positive")
  }

  C <- Matrix::Diagonal(dim(C)[1], rowSums(C))
  Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
  I <- Matrix::Diagonal(dim(C)[1])
  L <- L / scale.factor
  CiL <- Ci %*% L
  if (beta %% 1 == 0) { # not fractional
    Pr <- I
    Pl <- L
    if (beta > 1) {
      for (i in 1:beta) {
        Pl <- Pl %*% CiL
      }
    }
    Pl.roots <- Pl.factors <- NULL
    Pl.k <- beta
    Pl.scaling <- scale.factor^beta
    Pr.roots <- Pr.factors <- NULL
  } else {
    roots <- get.roots(m, beta)
    Pl.roots <- roots$rb
    Pr.roots <- roots$rc
    m_beta <- max(1, floor(beta))
    Pr.factors <- Pl.factors <- list()
    # construct Pl
    Pl <- I - CiL * roots$rb[1]
    Pl.factors[[1]] <- I - CiL * roots$rb[1]
    if (length(roots$rb) > 1) {
      for (i in 2:length(roots$rb)) {
        Pl <- Pl %*% (I - CiL * roots$rb[i])
        Pl.factors[[i]] <- I - CiL * roots$rb[i]
      }
    }
    Lp <- C
    if (m_beta > 1) {
      for (i in 1:(m_beta - 1)) {
        Lp <- Lp %*% CiL
      }
    }
    Pl.k <- m_beta - 1
    Pl.scaling <- scale.factor^beta / roots$factor
    Pl <- Lp %*% Pl

    # construct Pr
    Pr <- I - CiL * roots$rc[1]
    Pr.factors[[1]] <- I - CiL * roots$rc[1]
    if (length(roots$rc) > 1) {
      for (i in 2:length(roots$rc)) {
        Pr <- Pr %*% (I - CiL * roots$rc[i])
        Pr.factors[[i]] <- I - CiL * roots$rc[i]
      }
    }
  }
  Pl <- Pl * Pl.scaling
  Phi <- Matrix::Diagonal(dim(C)[1], 1 / tau)
  output <- list(
    Q = t(Pl) %*% Ci %*% Pl,
    Pl = Pl,
    Pr = Phi %*% Pr,
    Ci = Ci,
    C = C,
    CiL = CiL,
    L = L,
    Pl.factors = list(
      scaling = Pl.scaling,
      roots = Pl.roots,
      factor = Pl.factors,
      k = Pl.k
    ),
    Pr.factors = list(
      roots = Pr.roots,
      factor = Pr.factors,
      Phi = Phi
    ),
    m = m,
    beta = beta,
    type = "fractional approximation"
  )
  class(output) <- "rSPDEobj"
  return(output)
}

#' Rational approximations of stationary Gaussian Matern random fields
#'
#' `matern.operators` is used for computing a rational SPDE approximation
#' of a stationary Gaussian random fields on \eqn{R^d} with a Matern covariance
#' function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}
#' (\kappa h)^\nu K_\nu(\kappa h)}{C(h) = (\sigma^2/(2^{\nu-1}\Gamma(\nu))
#' (\kappa h)^\nu K_\nu(\kappa h).}
#'
#' @param kappa Parameter kappa of the covariance function.
#' @param range Range parameter of the covariance function (will be used if kappa is not provided).
#' @param sigma Standard deviation of the covariance function.
#' @param tau Parameter tau of the covariance function (will be used if sigma is not provided).
#' @param nu Shape parameter of the covariance function.
#' @param G The stiffness matrix of a finite element discretization of the
#' domain of interest. Does not need to be given if `mesh` is used.
#' @param C The mass matrix of a finite element discretization of the domain
#' of interest. Does not need to be given if `mesh` is used.
#' @param mesh An optional inla mesh. `d`, `C` and `G`
#' must be given if `mesh` is not given.
#' @param d The dimension of the domain. Does not need to be given if
#' `mesh` is used.
#' @param m The order of the rational approximation, which needs to be a
#' positive integer. The default value is 1.
#' @param type The type of the rational approximation. The options are
#' "covariance" and "operator". The default is "covariance".
#' @param compute_higher_order Logical. Should the higher order finite
#' element matrices be computed?
#' @param return_block_list Logical. For `type = "covariance"`,
#' should the block parts of the precision matrix be returned
#' separately as a list?
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are
#' "chebfun", "brasil" or "chebfunLB".
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @return If `type` is "covariance", then `matern.operators`
#' returns an object of class "CBrSPDEobj".
#' This object is a list containing the
#' following quantities:
#' \item{C}{The mass lumped mass matrix.}
#' \item{Ci}{The inverse of `C`.}
#' \item{GCi}{The stiffness matrix G times `Ci`}
#' \item{Gk}{The stiffness matrix G along with the higher-order
#' FEM-related matrices G2, G3, etc.}
#' \item{fem_mesh_matrices}{A list containing the mass
#' lumped mass matrix, the stiffness matrix and
#' the higher-order FEM-related matrices.}
#' \item{m}{The order of the rational approximation.}
#' \item{alpha}{The fractional power of the precision operator.}
#' \item{type}{String indicating the type of approximation.}
#' \item{d}{The dimension of the domain.}
#' \item{nu}{Shape parameter of the covariance function.}
#' \item{kappa}{Range parameter of the covariance function}
#' \item{tau}{Scale parameter of the covariance function.}
#' \item{sigma}{Standard deviation of the covariance function.}
#' \item{type}{String indicating the type of approximation.}
#'
#' If `type` is "operator", then `matern.operators`
#'  returns an object of class "rSPDEobj". This object contains the
#' quantities listed in the output of [fractional.operators()],
#' the `G` matrix, the dimension of the domain, as well as the
#' parameters of the covariance function.
#'
#' @details If `type` is "covariance", we use the
#' covariance-based rational approximation of the fractional operator.
#' In the SPDE approach, we model \eqn{u} as the solution of the following SPDE:
#' \deqn{L^{\alpha/2}(\tau u) = \mathcal{W},}
#' where
#' \eqn{L  = -\Delta +\kappa^2 I} and \eqn{\mathcal{W}} is the standard
#' Gaussian white noise. The covariance operator of \eqn{u} is given
#' by \eqn{L^{-\alpha}}. Now, let \eqn{L_h} be a finite-element
#' approximation of \eqn{L}. We can use
#' a rational approximation of order \eqn{m} on \eqn{L_h^{-\alpha}} to
#' obtain the following approximation:
#' \deqn{L_{h,m}^{-\alpha} = L_h^{-m_\alpha} p(L_h^{-1})q(L_h^{-1})^{-1},}
#' where \eqn{m_\alpha = \lfloor \alpha\rfloor}, \eqn{p} and \eqn{q} are
#' polynomials arising from such rational approximation.
#' From this approximation we construct an approximate precision
#' matrix for \eqn{u}.
#'
#' If `type` is "operator", the approximation is based on a
#' rational approximation of the fractional operator
#' \eqn{(\kappa^2 -\Delta)^\beta}, where \eqn{\beta = (\nu + d/2)/2}.
#' This results in an approximate model of the form \deqn{P_l u(s) = P_r W,}
#' where \eqn{P_j = p_j(L)} are non-fractional operators defined in terms
#' of polynomials \eqn{p_j} for \eqn{j=l,r}. The order of \eqn{p_r} is given
#' by `m` and the order of \eqn{p_l} is \eqn{m + m_\beta}
#' where \eqn{m_\beta} is the integer part of \eqn{\beta} if \eqn{\beta>1} and
#' \eqn{m_\beta = 1} otherwise.
#'
#' The discrete approximation can be written as \eqn{u = P_r x} where
#' \eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})} and \eqn{Q = P_l^T C^{-1} P_l}.
#' Note that the matrices \eqn{P_r} and \eqn{Q} may be be
#' ill-conditioned for \eqn{m>1}. In this case, the methods in
#' [operator.operations()] should be used for operations involving
#' the matrices, since these methods are more numerically stable.
#' @export
#' @seealso [fractional.operators()],
#' [spde.matern.operators()],
#' [matern.operators()]
#'
#' @examples
#' # Compute the covariance-based rational approximation of a
#' # Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' nobs <- 101
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' op_cov <- matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, sigma = sigma, d = 1, m = 2
#' )
#'
#' v <- t(rSPDE.A1d(x, 0.5))
#' # Compute the precision matrix
#' Q <- op_cov$Q
#' # A matrix here is the identity matrix
#' A <- Diagonal(nobs)
#' # We need to concatenate 3 A's since we are doing a covariance-based rational
#' # approximation of order 2
#' Abar <- cbind(A, A, A)
#' w <- rbind(v, v, v)
#' # The approximate covariance function:
#' c_cov.approx <- (Abar) %*% solve(Q, w)
#' c.true <- folded.matern.covariance.1d(rep(0.5, length(x)),
#'    abs(x), kappa, nu, sigma)
#'
#' # plot the result and compare with the true Matern covariance
#' plot(x, c.true,
#'   type = "l", ylab = "C(h)",
#'   xlab = "h", main = "Matern covariance and rational approximations"
#' )
#' lines(x, c_cov.approx, col = 2)
#'
#'
#' # Compute the operator-based rational approximation of a Gaussian
#' # process with a Matern covariance function on R
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
#'   kappa = kappa, sigma = sigma, nu = nu,
#'   G = fem$G, C = fem$C, d = 1,
#'   type = "operator"
#' )
#'
#' v <- t(rSPDE.A1d(x, 0.5))
#' c.approx <- Sigma.mult(op, v)
#' c.true <- folded.matern.covariance.1d(rep(0.5, length(x)),
#'   abs(x), kappa, nu, sigma)
#'
#' # plot the result and compare with the true Matern covariance
#' plot(x, c.true,
#'   type = "l", ylab = "C(h)",
#'   xlab = "h", main = "Matern covariance and rational approximation"
#' )
#' lines(x, c.approx, col = 2)
matern.operators <- function(kappa = NULL,
                             sigma = NULL,
                             tau = NULL,
                             range = NULL,
                             nu,
                             G = NULL,
                             C = NULL,
                             d = NULL,
                             mesh = NULL,
                             m = 1,
                             type = c("covariance", "operator"),
                             compute_higher_order = FALSE,
                             return_block_list = FALSE,
                             type_rational_approximation = c("chebfun",
                             "brasil", "chebfunLB"),
                             fem_mesh_matrices = NULL) {
  type <- type[[1]]
  nu <- min(nu, 10)
  if (!type %in% c("covariance", "operator")) {
    stop("The type should be 'covariance' or 'operator'!")
  }
  if (is.null(d) && is.null(mesh)) {
    stop("You should give either the dimension d or the mesh!")
  }

  if ((is.null(C) || is.null(G)) && is.null(mesh)) {
    stop("You should either provide mesh, or provide C *and* G!")
  }

  if(is.null(sigma)&&is.null(tau)){
    stop("You must either provide sigma or tau.")
  }

  if(is.null(range)&&is.null(kappa)){
    stop("You must either provide range or kappa.")
  }

    if (!is.null(mesh)) {
      d <- get_inla_mesh_dimension(inla_mesh = mesh)
      fem <- INLA::inla.mesh.fem(mesh)
      C <- fem$c0
      G <- fem$g1
    }

    if(!is.null(range)){
      kappa <- sqrt(8 * nu) / range
    }

    if(is.null(tau)){
      tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    }

    if(!is.null(tau)){
    sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2 * nu) *
    (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    }




  if (type == "operator") {
  
    beta <- (nu + d / 2) / 2
    operators <- fractional.operators(
      L = G + C * kappa^2,
      beta = beta,
      C = C,
      scale.factor = kappa^2,
      m = m,
      tau = tau
    )
    output <- operators
    output$kappa <- kappa
    output$range <- range
    output$tau <- tau
    output$sigma <- sigma
    output$nu <- nu
    output$d <- d
    output$G <- G
    output$type <- "Matern approximation"
    output$stationary <- TRUE
    return(output)
  } else {
    type_rational_approximation <- type_rational_approximation[[1]]
    out <- CBrSPDE.matern.operators(
      C = C, G = G, mesh = mesh, nu = nu, kappa = kappa, sigma = sigma,
      m = m, d = d, compute_higher_order = compute_higher_order,
      return_block_list = return_block_list,
      type_rational_approximation = type_rational_approximation
    )
    out$range <- range
    out$tau <- tau
    return(out)
  }
}



#' @name CBrSPDE.matern.operators
#' @title Covariance-based rational approximations of stationary Gaussian
#' Matern random fields
#' @description `CBrSPDE.matern.operators` is used for computing a
#' covariance-based rational SPDE approximation of stationary Gaussian random
#' fields on \eqn{R^d} with a Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h)}{C(h) =
#' (\sigma^2/(2^{\nu-1}\Gamma(\nu))(\kappa h)^\nu K_\nu(\kappa h)}
#' @param kappa Range parameter of the covariance function.
#' @param sigma Standard deviation of the covariance function.
#' @param nu Shape parameter of the covariance function.
#' @param G The stiffness matrix of a finite element discretization
#' of the domain of interest.
#' @param C The mass matrix of a finite element discretization of
#' the domain of interest.
#' @param mesh An inla mesh.
#' @param d The dimension of the domain.
#' @param m The order of the rational approximation, which needs
#' to be a positive integer.
#' The default value is 2.
#' @return `CBrSPDE.matern.operators` returns an object of
#' class "CBrSPDEobj". This object is a list containing the
#' following quantities:
#' \item{C}{The mass lumped mass matrix.}
#' \item{Ci}{The inverse of `C`.}
#' \item{GCi}{The stiffness matrix G times `Ci`}
#' \item{Gk}{The stiffness matrix G along with the higher-order
#' FEM-related matrices G2, G3, etc.}
#' \item{fem_mesh_matrices}{A list containing the mass lumped mass
#' matrix, the stiffness matrix and
#' the higher-order FEM-related matrices.}
#' \item{m}{The order of the rational approximation.}
#' \item{alpha}{The fractional power of the precision operator.}
#' \item{type}{String indicating the type of approximation.}
#' \item{d}{The dimension of the domain.}
#' \item{nu}{Shape parameter of the covariance function.}
#' \item{kappa}{Range parameter of the covariance function}
#' \item{tau}{Scale parameter of the covariance function.}
#' \item{sigma}{Standard deviation of the covariance function.}
#' \item{type}{String indicating the type of approximation.}
#' @noRd
#' @seealso [matern.operators()], [spde.matern.operators()]
#' @details We use the covariance-based rational approximation of the
#' fractional operator. In the SPDE approach, we model \eqn{u} as the
#' solution of the following SPDE:
#' \deqn{L^{\alpha/2}(\tau u) = \mathcal{W},}
#' where
#' \eqn{L  = -\Delta +\kappa^2 I} and \eqn{\mathcal{W}} is
#' the standard Gaussian white noise. The covariance operator of \eqn{u}
#' is given by \eqn{L^{-\alpha}}. Now, let \eqn{L_h} be a finite-element
#' approximation of \eqn{L}. We can use a rational approximation of order
#' \eqn{m} on \eqn{L_h^{-\alpha}} to obtain the following approximation:
#' \deqn{L_{h,m}^{-\alpha} = L_h^{-m_\alpha}
#' p(L_h^{-1})q(L_h^{-1})^{-1},}
#' where \eqn{m_\alpha = \max\{1,\lfloor \alpha\rfloor\}},
#' \eqn{p} and \eqn{q} are polynomials arising from such rational approximation.
#' From this approximation we construct an approximate precision
#' matrix for \eqn{u}.
#'
#'
#' @examples
#' # Compute the covariance-based rational approximation of a
#' # Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' nobs <- 101
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#' (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#' op_cov <- CBrSPDE.matern.operators(
#'   C = fem$C, G = fem$G, nu = nu,
#'   kappa = kappa, tau = tau, d = 1, m = 2
#' )
#'
#' v <- t(rSPDE.A1d(x, 0.5))
#' # Compute the precision matrix
#' Q <- op_cov$Q
#' # A matrix here is the identity matrix
#' A <- Diagonal(nobs)
#' # We need to concatenate 3 A's since we are doing a covariance-based rational
#' # approximation of order 2
#' Abar <- cbind(A, A, A)
#' w <- rbind(v, v, v)
#' # The approximate covariance function:
#' c_cov.approx <- (Abar) %*% solve(Q, w)
#'
#' # plot the result and compare with the true Matern covariance
#' plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
#'   type = "l", ylab = "C(h)",
#'   xlab = "h", main = "Matern covariance and rational approximations"
#' )
#' lines(x, c_cov.approx, col = 2)
CBrSPDE.matern.operators <- function(C,
                                     G,
                                     mesh,
                                     nu,
                                     kappa,
                                     sigma,
                                     m = 2,
                                     d,
                                     compute_higher_order = FALSE,
                                     return_block_list = FALSE,
                                     type_rational_approximation = c("chebfun",
                                     "brasil", "chebfunLB"),
                                     fem_mesh_matrices = NULL) {
  type_rational_approximation <- type_rational_approximation[[1]]

  if (is.null(fem_mesh_matrices)) {
    if (!is.null(mesh)) {
      ## get alpha, m_alpha
      d <- get_inla_mesh_dimension(inla_mesh = mesh)
      alpha <- nu + d / 2
      m_alpha <- floor(alpha)
      m_order <- m_alpha + 1
      tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
      (4 * pi)^(d / 2) * gamma(nu + d / 2)))

      if (d > 1) {
        if (compute_higher_order) {
          fem <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
        } else {
          fem <- INLA::inla.mesh.fem(mesh)
        }

        C <- fem$c0
        G <- fem$g1
        Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
        GCi <- G %*% Ci

        Gk <- list()
        Gk[[1]] <- G
        for (i in 2:m_order) {
          Gk[[i]] <- fem[[paste0("g", i)]]
        }
      } else if (d == 1) {
        fem <- INLA::inla.mesh.fem(mesh, order = 2)
        C <- fem$c0
        G <- fem$g1
        Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
        GCi <- G %*% Ci

        Gk <- list()
        Gk[[1]] <- G
        if (compute_higher_order) {
          Gk[[2]] <- fem$g2
          if (m_order > 2) {
            for (i in 3:m_order) {
              Gk[[i]] <- GCi %*% Gk[[i - 1]]
            }
          }
        }
      }
    } else {
      ## get alpha, m_alpha
      alpha <- nu + d / 2
      m_alpha <- floor(alpha)
      m_order <- m_alpha + 1
      tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
      (4 * pi)^(d / 2) * gamma(nu + d / 2)))

      ## get lumped mass matrix
      C <- Matrix::Diagonal(dim(C)[1], rowSums(C))

      ## get G_k matrix: k is up to m_alpha if alpha is integer,
      # k is up to m_alpha + 1 otherwise.
      # inverse lumped mass matrix
      Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))

      GCi <- G %*% Ci
      # create a list to store all the G_k matrix

      Gk <- list()

      Gk[[1]] <- G
      # determine how many G_k matrices we want to create
      if (compute_higher_order) {
        for (i in 2:m_order) {
          Gk[[i]] <- GCi %*% Gk[[i - 1]]
        }
      }
    }

    # create a list contains all the finite element related matrices
    fem_mesh_matrices <- list()
    fem_mesh_matrices[["c0"]] <- C
    fem_mesh_matrices[["g1"]] <- G

    if (compute_higher_order) {
      for (i in 1:m_order) {
        fem_mesh_matrices[[paste0("g", i)]] <- Gk[[i]]
      }
    }
  } else {
    C <- fem_mesh_matrices$c0
    G <- fem_mesh_matrices$g1
    Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
    GCi <- G %*% Ci
    alpha <- nu + d / 2
    m_alpha <- floor(alpha)
    m_order <- m_alpha + 1
    if (!is.null(mesh)) {
      d <- get_inla_mesh_dimension(inla_mesh = mesh)
    }
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    Gk <- list()
    Gk[[1]] <- G
    for (i in 2:m_order) {
      Gk[[i]] <- fem_mesh_matrices[[paste0("g", i)]]
    }
  }


  L <- (G + kappa^2 * C) / kappa^2

  Lchol <- chol(L)
  logdetL <- 2 * sum(log(diag(Lchol)))

  logdetC <- sum(log(diag(C)))

  CiL <- GCi / kappa^2 + Diagonal(dim(GCi)[1])

  if (m_alpha == 0) {
    aux_mat <- Diagonal(dim(L)[1])
  } else {
    aux_mat <- CiL
  }


  if (return_block_list) {
    Q.int <- aux_mat

    if(alpha %% 1 == 0){
      Q.frac <- Matrix::Diagonal(dim(L)[1])
      Q <- L * kappa^2

      if(alpha > 1){
        CinvL <- Matrix::Diagonal(dim(C)[1], 1/diag(C)) %*% L * kappa^2

        for(k in 1:(alpha-1)){
          Q <- Q %*% CinvL
        }
      } 

      Q <- tau^2 * Q
      Q.int <- Q
    } else{
      if(m == 0){
        stop("Return block list does not work with m = 0, either increase m or set return_block_list to FALSE.")
      }
      Q.frac <- rspde.matern.precision(
        kappa = kappa, nu = nu, tau = tau,
        rspde.order = m, dim = d,
        fem_mesh_matrices = fem_mesh_matrices, only_fractional = TRUE,
        return_block_list = TRUE,
        type_rational_approx = type_rational_approximation
      )

     Q <- Q.frac

      if (m_alpha > 0) {
        for (j in seq_len(length(Q))) {
          for (i in 1:m_alpha) {
            Q[[j]] <- Q.int %*% Q[[j]]
          }
        }
      }
      Q.int <- list(Q.int = Q.int, order = m_alpha)
    }

  } else {
    Q.int <- list(Q.int = kronecker(Diagonal(m + 1), aux_mat), order = m_alpha)

    if(alpha %% 1 == 0){
      Q.frac <- Matrix::Diagonal(dim(L)[1])
      Q <- L * kappa^2

      if(alpha > 1){
        CinvL <- Matrix::Diagonal(dim(C)[1], 1/diag(C)) %*% L * kappa^2

        for(k in 1:(alpha-1)){
          Q <- Q %*% CinvL
        }
      } 

      Q <- tau^2 * Q
      Q.int <- list(Q.int = Q, order = m_alpha)
    } else if (m > 0){
        Q.frac <- rspde.matern.precision(
          kappa = kappa, nu = nu, tau = tau,
          rspde.order = m, dim = d,
          fem_mesh_matrices = fem_mesh_matrices, only_fractional = TRUE,
          type_rational_approx = type_rational_approximation
        )

        Q <- Q.frac

        if (m_alpha > 0) {
          for (i in 1:m_alpha) {
            Q <- Q.int$Q.int %*% Q
          }
        }
    } else{
      Q <- markov.Q(alpha/2, kappa, d, list(C=C, G=G))
      Q.frac <- Matrix::Diagonal(dim(L)[1])
    }
  }

  ## output
  output <- list(
    C = C, G = G, L = L, Ci = Ci, GCi = GCi, Gk = Gk,
    fem_mesh_matrices = fem_mesh_matrices,
    alpha = alpha, nu = nu, kappa = kappa,
    tau = tau, m = m, d = d,
    sigma = sigma,
    logdetL = logdetL, logdetC = logdetC,
    Q.frac = Q.frac, Q.int = Q.int,
    Q = Q, sizeC = dim(C)[1],
    higher_order = compute_higher_order,
    type_rational_approximation = type_rational_approximation,
    return_block_list = return_block_list,
    stationary = TRUE
  )
  output$type <- "Covariance-Based Matern SPDE Approximation"
  class(output) <- "CBrSPDEobj"
  return(output)
}


#' Rational approximations of non-stationary Gaussian SPDE Matern random fields
#'
#' `spde.matern.operators` is used for computing a rational SPDE
#' approximation of a Gaussian random
#' fields on \eqn{R^d} defined as a solution to the SPDE
#' \deqn{(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) = W.}
#'
#' @param kappa Vector with the, possibly spatially varying, range parameter
#' evaluated at the locations of the mesh used for the finite element
#' discretization of the SPDE.
#' @param tau Vector with the, possibly spatially varying, precision
#' parameter evaluated at the locations
#' of the mesh used for the finite element discretization of the SPDE.
#' @param theta Theta parameter that connects B.tau and B.kappa to tau and kappa through a log-linear regression, in case the parameterization is `spde`,
#' and that connects B.sigma and B.range to tau and kappa in case the parameterization is `matern`.
#' @param B.sigma Matrix with specification of log-linear model for \eqn{\sigma}. Will be used if `parameterization = 'matern'`.
#' @param B.range Matrix with specification of log-linear model for \eqn{\rho}, which is a range-like parameter (it is exactly the range parameter in the stationary case). Will be used if `parameterization = 'matern'`.
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The default is `matern`.
#' @param B.tau Matrix with specification of log-linear model for \eqn{\tau}. Will be used if `parameterization = 'spde'`.
#' @param B.kappa Matrix with specification of log-linear model for \eqn{\kappa}. Will be used if `parameterization = 'spde'`.
#' @param nu Shape parameter of the covariance function, related to
#' \eqn{\beta} through the equation
#' \eqn{\beta = (\nu + d/2)/2}.
#' @param G The stiffness matrix of a finite element discretization of
#' the domain of interest.
#' @param C The mass matrix of a finite element discretization of the
#' domain of interest.
#' @param mesh An optional inla mesh. `d`, `C` and `G`
#' must be given if `mesh` is not given.
#' @param d The dimension of the domain. Does not need to be given if
#' `mesh` is used.
#' @param m The order of the rational approximation, which needs to be a
#' positive integer. The default value is 1.
#' @param type The type of the rational approximation. The options are
#' "covariance" and "operator". The default is "covariance".
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are
#' "chebfun", "brasil" or "chebfunLB".
#'
#'
#' @details The approximation is based on a rational approximation of the
#' fractional operator \eqn{(\kappa(s)^2 -\Delta)^\beta}, where
#' \eqn{\beta = (\nu + d/2)/2}. This results in an approximate model
#' on the form \deqn{P_l u(s) = P_r W,} where \eqn{P_j = p_j(L)} are
#' non-fractional operators defined in terms of polynomials \eqn{p_j} for
#' \eqn{j=l,r}. The order of \eqn{p_r} is given by `m` and the order
#' of \eqn{p_l} is \eqn{m + m_\beta} where \eqn{m_\beta} is the integer
#' part of \eqn{\beta} if \eqn{\beta>1} and \eqn{m_\beta = 1} otherwise.
#'
#' The discrete approximation can be written as \eqn{u = P_r x} where
#' \eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})}
#' and \eqn{Q = P_l^T C^{-1} P_l}. Note that the matrices \eqn{P_r} and
#' \eqn{Q} may be be ill-conditioned for \eqn{m>1}.
#' In this case, the metehods in [operator.operations()]
#' should be used for operations involving the matrices, since
#' these methods are more numerically stable.
#'
#' @return `spde.matern.operators` returns an object of
#' class "rSPDEobj. This object contains the
#' quantities listed in the output of [fractional.operators()]
#' as well as the smoothness parameter \eqn{\nu}.
#' @export
#' @seealso [fractional.operators()],
#' [spde.matern.operators()],
#' [matern.operators()]
#'
#' @examples
#' # Sample non-stationary Matern field on R
#' tau <- 1
#' nu <- 0.8
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # define a non-stationary range parameter
#' kappa <- seq(from = 2, to = 20, length.out = length(x))
#'
#' # compute rational approximation
#' op <- spde.matern.operators(
#'   kappa = kappa, tau = tau, nu = nu,
#'   G = fem$G, C = fem$C, d = 1
#' )
#'
#' # sample the field
#' u <- simulate(op)
#'
#' # plot the sample
#' plot(x, u, type = "l", ylab = "u(s)", xlab = "s")
#'
spde.matern.operators <- function(kappa = NULL,
                                  tau = NULL,
                                  theta = NULL,
                                  B.tau = matrix(c(0, 1, 0), 1, 3), 
                                  B.kappa = matrix(c(0, 0, 1), 1, 3), 
                                  B.sigma = matrix(c(0, 1, 0), 1, 3), 
                                  B.range = matrix(c(0, 0, 1), 1, 3),
                                  nu,
                                  parameterization = c("matern", "spde"),
                                  G = NULL,
                                  C = NULL,
                                  d = NULL,
                                  mesh = NULL,
                                  m = 1,
                                  type = c("covariance", "operator"),
                                  type_rational_approximation = c("chebfun",
                             "brasil", "chebfunLB")) {
  type <- type[[1]]
  if (!type %in% c("covariance", "operator")) {
    stop("The type should be 'covariance' or 'operator'!")
  }

  parameterization <- parameterization[[1]]

  if (!parameterization %in% c("matern", "spde")) {
    stop("parameterization should be either 'matern' or 'spde'!")
  }

  if (is.null(d) && is.null(mesh)) {
    stop("You should give either the dimension d or the mesh!")
  }

  if ((is.null(C) || is.null(G)) && is.null(mesh)) {
    stop("You should either provide mesh, or provide C *and* G!")
  }

  if (!is.null(mesh)) {
      d <- get_inla_mesh_dimension(inla_mesh = mesh)
      fem <- INLA::inla.mesh.fem(mesh)
      C <- fem$c0
      G <- fem$g1
  }

  if(!is.null(theta)){

    if(any(dim(B.tau) != dim(B.kappa))){
      stop("B.tau and B.kappa must have the same dimensions!")
    }
    if(any(dim(B.sigma) != dim(B.range))){
      stop("B.sigma and B.range must have the same dimensions!")
    }

    if(parameterization == "matern"){
      B_matrices <- convert_B_matrices(B.sigma, B.range, ncol(C), nu, d)
      B.tau <- B_matrices[["B.tau"]]
      B.kappa <- B_matrices[["B.kappa"]]
    } else{
      B.tau <- prepare_B_matrices(B.tau, ncol(C), 
          ncol(B.tau)-1)
      B.kappa <- prepare_B_matrices(B.kappa, ncol(C), 
          ncol(B.kappa)-1)
    }
    
      new_theta <- c(1, theta)

      tau <- exp(B.tau %*% new_theta)
      kappa <- exp(B.kappa %*% new_theta)

  } else if(is.null(kappa) || is.null(tau)){
    stop("If theta is NULL, then kappa and tau must not be NULL!")
  } 

  alpha <- nu + d / 2
  
  if (nu < 0) {
    stop("nu must be positive")
  }
  
  if(type == "operator"){
      beta <- (nu + d / 2) / 2
      kp <- Matrix::Diagonal(dim(C)[1], kappa^2)

      operators <- fractional.operators(
        L = G + C %*% kp,
        beta = beta,
        C = C,
        scale.factor = min(kappa)^2,
        m = m,
        tau = tau
      )
      output <- operators
      output$beta <- beta
      output$B.tau <- B.tau
      output$B.kappa <- B.kappa
      output$kappa <- kappa
      output$tau <- tau
      output$theta <- theta
      output$stationary <- FALSE
      output$type <- "Matern SPDE approximation"
  } else{
    type_rational_approximation <- type_rational_approximation[[1]]
    m_alpha <- floor(alpha)
    m_order <- m_alpha + 1

    L <- Matrix::Diagonal(dim(C)[1], kappa^2 * diag(C))
    L <- L + G

    tau_matrix <- Matrix::Diagonal(dim(C)[1], tau)
    
    if(alpha %% 1 == 0){
      Q <- L

      if(alpha > 1){
        CinvL <- Matrix::Diagonal(dim(C)[1], 1/diag(C)) %*% L

        for(k in 1:(alpha-1)){
          Q <- Q %*% CinvL
        }
      } 

      Q <- tau_matrix %*% Q %*% tau_matrix
    } else{
      factor <- min(kappa^2)
      L <- L/factor

      if(m_alpha > 0){
        CinvL <- Matrix::Diagonal(dim(C)[1], 1/diag(C)) %*% L
      }

      mt <- get_rational_coefficients(m, type_rational_approximation)
      rspde.order <- m
      row_nu <- round(1000*cut_decimals(alpha))
      r <- unlist(mt[row_nu, 2:(1+rspde.order)])
      p <- unlist(mt[row_nu, (2+rspde.order):(1+2*rspde.order)])
      k_rat <- unlist(mt[row_nu, 2+2*rspde.order])

      Q <- (L - p[1] * C)/r[1]
        if(m_alpha > 0){
          for(m in 1:m_alpha){
            Q <- Q %*% CinvL
          }
        }

      Q <- tau_matrix %*% Q %*% tau_matrix

      if(rspde.order > 1){
        for(k_ind in 2:rspde.order){
          Q_tmp <- (L - p[k_ind] * C)/r[k_ind]
          if(m_alpha > 0){
            for(count_m in 1:m_alpha){
              Q_tmp <- Q_tmp %*% CinvL
            }
          }
        }
        Q_tmp <- tau_matrix %*% Q_tmp %*% tau_matrix
        Q <- bdiag(Q, Q_tmp)
      }

      # K part

      if(m_alpha == 0){
        Q_tmp <- C
      } else if(m_alpha == 1){
        Q_tmp <- L
      } else{
        Q_tmp <- L
        for(count_m in 1:(m_alpha-1)){
          Q_tmp <- Q_tmp %*% CinvL
        }
      }

      Q_tmp <- tau_matrix %*% Q_tmp %*% tau_matrix
      Q_tmp <- Q_tmp/k_rat

      Q <- bdiag(Q, Q_tmp)

      Q <- factor^(alpha) * Q
    }
  

  ## output
  output <- list(
    C = C, G = G,
    alpha = alpha, nu = nu, kappa = kappa,
    B.tau = B.tau, B.kappa = B.kappa,
    tau = tau, theta = theta, m = m, d = d,
    Q = Q, sizeC = dim(C)[1],
    type_rational_approximation = type_rational_approximation,
    stationary = FALSE
  )
  output$type <- "Covariance-Based Matern SPDE Approximation"
  class(output) <- "CBrSPDEobj"
  }

  return(output)
}

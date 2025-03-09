## util.R
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


#' Observation matrix for finite element discretization on R
#'
#' A finite element discretization on R can be written as
#' \eqn{u(s) = \sum_i^n u_i \varphi_i(s)}{u(s) = \sum_i^n u_i \varphi_i(s)}
#' where \eqn{\varphi_i(s)} is a piecewise linear
#' "hat function" centered at location
#' \eqn{x_i}{x_i}. This function computes an
#' \eqn{m\times n}{m x n} matrix \eqn{A}{A}
#' that links the basis function in the expansion to specified locations
#' \eqn{s = (s_1,\ldots, s_m)} in the domain through
#' \eqn{A_ij = \varphi_j(s_i)}{A_ij = \varphi_j(s_i)}.
#'
#' @param x The locations of the nodes in the FEM discretization.
#' @param loc The locations \eqn{(s_1,\ldots, s_m)}
#'
#' @return The sparse matrix `A`.
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso [rSPDE.fem1d()]
#'
#' @examples
#' # create mass and stiffness matrices for a FEM discretization on [0,1]
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # create the observation matrix for some locations in the domain
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
rSPDE.A1d <- function(x, loc) {
  if (min(loc) < min(x) || max(loc) > max(x)) {
    stop("locations outside support of basis")
  }

  n.x <- length(x)
  n.loc <- length(loc)
  i <- as.vector(cbind(1:n.loc, 1:n.loc))
  j <- matrix(0, n.loc, 2)
  vals <- matrix(1, n.loc, 2)
  for (ii in seq_len(n.loc)) {
    j[ii, 1] <- sum(sum((loc[ii] - x) >= 0))
    vals[ii, 1] <- loc[ii] - x[j[ii, 1]]
    j[ii, 2] <- j[ii, 1] + 1
    if (j[ii, 2] <= n.x) {
      vals[ii, 2] <- x[j[ii, 2]] - loc[ii]
    } else {
      j[ii, 2] <- j[ii, 2] - 2
    }
  }
  j <- as.vector(j)
  vals <- as.vector(matrix(1 - vals / rowSums(vals)))

  A <- sparseMatrix(i = i, j = j, x = vals, dims = c(n.loc, n.x))
  return(A)
}


#' Finite element calculations for problems on R
#'
#' This function computes mass and stiffness matrices
#' for a FEM approximation on R, assuming
#' Neumann boundary conditions.
#' These matrices are needed when discretizing the
#' operators in rational approximations.
#'
#' @param x Locations of the nodes in the FEM approximation.
#'
#' @return The function returns a list with the following elements
#' \item{G }{The stiffness matrix with elements \eqn{(\nabla \phi_i, \nabla \phi_j)}.}
#' \item{C }{The mass matrix with elements \eqn{(\phi_i, \phi_j)}.}
#' \item{Cd }{Mass lumped mass matrix.}
#' \item{B }{Matrix with elements \eqn{(\nabla \phi_i, \phi_j)}.}
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso [rSPDE.A1d()]
#' @examples
#' # create mass and stiffness matrices for a FEM discretization on [0,1]
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
rSPDE.fem1d <- function(x) {
  n <- length(x)
  d <- c(Inf, diff(x))
  dm1 <- c(d[2:n], Inf)
  G <- -bandSparse(
    n = n, m = n, k = c(-1, 0, 1),
    diagonals = cbind(1 / dm1, -(1 / dm1 + 1 / d), 1 / dm1)
  )
  C <- bandSparse(
    n = n, m = n, k = c(-1, 0, 1),
    diagonals = cbind(dm1 / 6, (dm1 + d) / 3, c(d[2:n], Inf) / 6)
  )
  C[1, 1:2] <- c(d[2], d[2] / 2) / 3
  C[n, (n - 1):n] <- c(d[n] / 2, d[n]) / 3

  Cd <- Diagonal(rowSums(C),n=n)

  B <- bandSparse(n = n, m = n, k = c(-1, 0, 1),
                  diagonals = cbind(rep(0.5,n), rep(0,n), rep(-0.5,n)))

  return(list(G = G, C = C, Cd = Cd, B = B))
}

#' Finite element calculations for problems in 2D
#'
#' This function computes mass and stiffness matrices for a mesh in 2D, assuming
#' Neumann boundary conditions.
#'
#' @param FV Matrix where each row defines a triangle
#' @param P Locations of the nodes in the mesh.
#'
#' @return The function returns a list with the following elements
#' \item{G }{The stiffness matrix with elements \eqn{(\nabla \phi_i, \nabla \phi_j)}.}
#' \item{C }{The mass matrix with elements \eqn{(\phi_i, \phi_j)}.}
#' \item{Cd }{The mass lumped matrix with diagonal elements \eqn{(\phi_i, 1)}.}
#' \item{Hxx }{Matrix with elements \eqn{(\partial_x \phi_i, \partial_x \phi_j)}.}
#' \item{Hyy }{Matrix with elements \eqn{(\partial_y \phi_i, \partial_y \phi_j)}.}
#' \item{Hxy }{Matrix with elements \eqn{(\partial_x \phi_i, \partial_y \phi_j)}.}
#' \item{Hyx }{Matrix with elements \eqn{(\partial_y \phi_i, \partial_x \phi_j)}.}
#' \item{Bx }{Matrix with elements \eqn{(\partial_x \phi_i, \phi_j)}.}
#' \item{By }{Matrix with elements \eqn{(\partial_y \phi_i, \phi_j)}.}
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso [rSPDE.fem1d()]
#' @examples
#' P <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' FV <- rbind(c(1, 2, 3), c(2, 3, 4))
#' fem <- rSPDE.fem2d(FV, P)
rSPDE.fem2d <- function(FV, P) {
  d <- ncol(FV) - 1
  if (d != 2) {
    stop("Only 2d supported")
  }
  if (ncol(P) != d) {
    P <- t(P)
  }
  if (ncol(P) != d) {
    stop("Wrong dimension of P")
  }

  nV <- nrow(P)
  nF <- nrow(FV)
  Gi <- matrix(0, nrow = nF * 3, ncol = 3)
  Gj <- Gz <- Ci <- Cj <- Cz <- Gxx <- Gxy <- Gyx <- Gyy <- Bxz <- Byz <- Gi

  Mxx <- matrix(c(1, -1, 0, -1, 1, 0, 0, 0, 0), 3, 3)
  Myy <- matrix(c(1, 0, -1, 0, 0, 0, -1, 0, 1), 3, 3)
  Mxy <- matrix(c(1, -1, 0, 0, 0, 0, -1, 1, 0), 3, 3)
  Myx <- matrix(c(1, 0, -1, -1, 0, 1, 0, 0, 0), 3, 3)
  for (f in 1:nF) {
    dd <- 3 * (f - 1) + (1:3)
    Gi[dd, ] <- Ci[dd, ] <- FV[f, ] %*% t(rep(1, 3))
    Gj[dd, ] <- Cj[dd, ] <- t(Gi[dd, ])

    xy <- t(P[FV[f, ], ])
    m1 <- rbind(rep(1, 3), xy)
    m2 <- rbind(rep(0, 2), diag(1, 2))
    m <- solve(m1, m2)
    ddet <- abs(det(m1))
    Gz[dd, ] <- ddet * (m %*% t(m)) / 2
    Cz[dd, ] <- ddet * (rep(1, 3) + diag(3)) / 24

    Bk <- matrix(c(
      xy[1, 2] - xy[1, 1],
      xy[2, 2] - xy[2, 1],
      xy[1, 3] - xy[1, 1],
      xy[2, 3] - xy[2, 1]
    ), 2, 2)

    Bki <- solve(Bk)
    Cxx <- Bki %*% matrix(c(1, 0, 0, 0), 2, 2) %*% t(Bki)
    Cyy <- Bki %*% matrix(c(0, 0, 0, 1), 2, 2) %*% t(Bki)
    Cxy <- Bki %*% matrix(c(0, 0, 1, 0), 2, 2) %*% t(Bki)
    Cyx <- Bki %*% matrix(c(0, 1, 0, 0), 2, 2) %*% t(Bki)

    Gxx[dd, ] <- ddet * (Cxx[1, 1] * Mxx + Cxx[1, 2] * Mxy + Cxx[2, 1] * Myx + Cxx[2, 2] * Myy) / 2
    Gyy[dd, ] <- ddet * (Cyy[1, 1] * Mxx + Cyy[1, 2] * Mxy + Cyy[2, 1] * Myx + Cyy[2, 2] * Myy) / 2
    Gxy[dd, ] <- ddet * (Cxy[1, 1] * Mxx + Cxy[1, 2] * Mxy + Cxy[2, 1] * Myx + Cxy[2, 2] * Myy) / 2
    Gyx[dd, ] <- ddet * (Cyx[1, 1] * Mxx + Cyx[1, 2] * Mxy + Cyx[2, 1] * Myx + Cyx[2, 2] * Myy) / 2

    ab1 <- solve(matrix(c(xy[1,2]-xy[1,1], xy[1,3]-xy[1,1],
                          xy[2,2]-xy[2,1], xy[2,3]-xy[2,1]),2,2),rep(1,2))
    ab2 <- solve(matrix(c(xy[1,1]-xy[1,2], xy[1,3]-xy[1,2],
                          xy[2,1]-xy[2,2], xy[2,3]-xy[2,2]),2,2),rep(1,2))
    ab3 <- solve(matrix(c(xy[1,1]-xy[1,3], xy[1,2]-xy[1,3],
                          xy[2,1]-xy[2,3], xy[2,2]-xy[2,3]),2,2),rep(1,2))

    Bxz[dd, ] <-  -c(ab1[1], ab2[1], ab3[1]) * ddet/6
    Byz[dd, ] <-  -c(ab1[2], ab2[2], ab3[2]) * ddet/6

  }

  G <- Matrix::sparseMatrix(
    i = as.vector(Gi), j = as.vector(Gj),
    x = as.vector(Gz), dims = c(nV, nV)
  )
  Hxx <- Matrix::sparseMatrix(
    i = as.vector(Gi), j = as.vector(Gj),
    x = as.vector(Gxx), dims = c(nV, nV)
  )
  Hyy <- Matrix::sparseMatrix(
    i = as.vector(Gi), j = as.vector(Gj),
    x = as.vector(Gyy), dims = c(nV, nV)
  )
  Hxy <- Matrix::sparseMatrix(
    i = as.vector(Gi), j = as.vector(Gj),
    x = as.vector(Gxy), dims = c(nV, nV)
  )
  Hyx <- Matrix::sparseMatrix(
    i = as.vector(Gi), j = as.vector(Gj),
    x = as.vector(Gyx), dims = c(nV, nV)
  )

  Bx <- Matrix::sparseMatrix(
      i = as.vector(Gi), j = as.vector(Gj),
      x = as.vector(Bxz), dims = c(nV, nV)
  )

  By <- Matrix::sparseMatrix(
      i = as.vector(Gi), j = as.vector(Gj),
      x = as.vector(Byz), dims = c(nV, nV)
  )

  Ce <- Matrix::sparseMatrix(
    i = as.vector(Ci), j = as.vector(Cj),
    x = as.vector(Cz), dims = c(nV, nV)
  )
  C <- Matrix::Diagonal(n = nV, x = Matrix::colSums(Ce))
  return(list(G = G, C = Ce, Cd = C,
              Hxx = Hxx, Hyy = Hyy, Hxy = Hxy, Hyx = Hyx,
              Bx = Bx, By = By))
}
#' Warnings free loading of add-on packages
#'
#' Turn off all warnings for require(), to allow clean completion
#' of examples that require unavailable Suggested packages.
#'
#' @param package The name of a package, given as a character string.
#' @param lib.loc a character vector describing the location of R library trees
#' to search through, or `NULL`.  The default value of `NULL`
#' corresponds to all libraries currently known to `.libPaths()`.
#' Non-existent library trees are silently ignored.
#' @param character.only a logical indicating whether `package` can be
#' assumed to be a character string.
#'
#' @return `require.nowarnings` returns (invisibly)
#' `TRUE` if it succeeds, otherwise `FALSE`
#' @details `require(package)` acts the same as
#' `require(package, quietly = TRUE)` but with warnings turned off.
#' In particular, no warning or error is given if the package is unavailable.
#' Most cases should use `requireNamespace(package,
#' quietly = TRUE)` instead,
#' which doesn't produce warnings.
#' @seealso [require()]
#' @export
#' @examples
#' ## This should produce no output:
#' if (require.nowarnings(nonexistent)) {
#'   message("Package loaded successfully")
#' }
#'
require.nowarnings <- function(
    package, lib.loc = NULL,
    character.only = FALSE) {
  if (!character.only) {
    package <- as.character(substitute(package))
  }
  suppressWarnings(
    require(package,
      lib.loc = lib.loc,
      quietly = TRUE,
      character.only = TRUE
    )
  )
}

#' @name get.initial.values.rSPDE
#' @title Initial values for log-likelihood optimization in rSPDE models
#' with a latent stationary Gaussian Matern model
#' @description Auxiliar function to obtain domain-based initial values for
#' log-likelihood optimization in rSPDE models
#' with a latent stationary Gaussian Matern model
#' @param mesh An in INLA mesh
#' @param mesh.range The range of the mesh.
#' @param graph.obj A `metric_graph` object. To be used in case both `mesh` and `mesh.range` are `NULL`.
#' @param dim The dimension of the domain.
#' @param B.sigma Matrix with specification of log-linear model for \eqn{\sigma}. Will be used if `parameterization = 'matern'`.
#' @param B.range Matrix with specification of log-linear model for \eqn{\rho}, which is a range-like parameter (it is exactly the range parameter in the stationary case). Will be used if `parameterization = 'matern'`.
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu (smoothness). `spde` uses kappa, tau and nu (smoothness). The default is `matern`.
#' @param B.tau Matrix with specification of log-linear model for \eqn{\tau}. Will be used if `parameterization = 'spde'`.
#' @param B.kappa Matrix with specification of log-linear model for \eqn{\kappa}. Will be used if `parameterization = 'spde'`.
#' @param nu The smoothness parameter.
#' @param include.nu Should we also provide an initial guess for nu?
#' @param n.spde The number of basis functions in the mesh model.
#' @param log.scale Should the results be provided in log scale?
#' @param nu.upper.bound Should an upper bound for nu be considered?
#' @return A vector of the form (theta_1,theta_2,theta_3) or where
#' theta_1 is the initial guess for tau, theta_2 is the initial guess for kappa
#' and theta_3 is the initial guess for nu.
#' @export
#'

get.initial.values.rSPDE <- function(mesh = NULL, mesh.range = NULL,
                                     graph.obj = NULL,
                                     n.spde = 1,
                                     dim = NULL, B.tau = NULL, B.kappa = NULL,
                                     B.sigma = NULL, B.range = NULL, nu = NULL,
                                     parameterization = c("matern", "spde"),
                                     include.nu = TRUE, log.scale = TRUE,
                                     nu.upper.bound = NULL) {
  if (is.null(mesh) && is.null(mesh.range) && is.null(graph.obj)) {
    stop("You should either provide mesh, mesh.range or graph_obj!")
  }

  parameterization <- parameterization[[1]]

  if (!parameterization %in% c("matern", "spde")) {
    stop("parameterization should be either 'matern' or 'spde'!")
  }

  if (is.null(mesh) && is.null(graph.obj) && is.null(dim)) {
    stop("If you don't provide mesh, you have to provide dim!")
  }

  if (!is.null(mesh)) {
    if (!inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d"))) {
      stop("The mesh should be created using fmesher!")
    }

    dim <- fmesher::fm_manifold_dim(mesh)
  }

  if (!is.null(graph.obj)) {
    if (!inherits(graph.obj, "metric_graph")) {
      stop("graph_obj should be a metric_graph object.")
    }
    dim <- 1
  }

  if (include.nu) {
    if (!is.null(nu.upper.bound)) {
      nu <- min(1, nu.upper.bound / 2)
    } else {
      nu <- 1
    }
  } else {
    if (is.null(nu)) {
      stop("If include.nu is FALSE, then nu must be provided!")
    }
  }

  if (parameterization == "matern") {
    if (is.null(B.sigma)) {
      B.sigma <- matrix(c(0, 1, 0), 1, 3)
    }
    if (is.null(B.range)) {
      B.range <- matrix(c(0, 0, 1), 1, 3)
    }

    if (is.null(graph.obj)) {
      param <- get_parameters_rSPDE(
        mesh = mesh,
        alpha = nu + dim / 2,
        B.tau = B.tau,
        B.kappa = B.kappa,
        B.sigma = B.sigma,
        B.range = B.range,
        nu.nominal = nu,
        alpha.nominal = nu + dim / 2,
        parameterization = parameterization,
        prior.std.dev.nominal = 1,
        prior.range.nominal = NULL,
        prior.tau = NULL,
        prior.kappa = NULL,
        theta.prior.mean = NULL,
        theta.prior.prec = 0.1,
        mesh.range = mesh.range,
        d = dim,
        n.spde = n.spde
      )
    } else {
      param <- get_parameters_rSPDE_graph(
        graph_obj = graph.obj,
        alpha = nu + 1 / 2,
        B.tau = B.tau,
        B.kappa = B.kappa,
        B.sigma = B.sigma,
        B.range = B.range,
        nu.nominal = nu,
        alpha.nominal = nu + 1 / 2,
        parameterization = parameterization,
        prior.std.dev.nominal = 1,
        prior.range.nominal = NULL,
        prior.tau = NULL,
        prior.kappa = NULL,
        theta.prior.mean = NULL,
        theta.prior.prec = 0.1
      )
    }

    initial <- param$theta.prior.mean
  } else {
    if (is.null(B.tau)) {
      B.tau <- matrix(c(0, 1, 0), 1, 3)
    }
    if (is.null(B.kappa)) {
      B.kappa <- matrix(c(0, 0, 1), 1, 3)
    }
    if (is.null(graph.obj)) {
      param <- get_parameters_rSPDE(
        mesh = mesh,
        alpha = nu + dim / 2,
        B.tau = B.tau,
        B.kappa = B.kappa,
        B.sigma = B.sigma,
        B.range = B.range,
        nu.nominal = nu,
        alpha.nominal = nu + dim / 2,
        parameterization = parameterization,
        prior.std.dev.nominal = 1,
        prior.range.nominal = NULL,
        prior.tau = NULL,
        prior.kappa = NULL,
        theta.prior.mean = NULL,
        theta.prior.prec = 0.1,
        mesh.range = mesh.range,
        d = dim,
        n.spde = n.spde
      )
    } else {
      param <- get_parameters_rSPDE_graph(
        graph_obj = graph.obj,
        alpha = nu + 1 / 2,
        B.tau = B.tau,
        B.kappa = B.kappa,
        B.sigma = B.sigma,
        B.range = B.range,
        nu.nominal = nu,
        alpha.nominal = nu + 1 / 2,
        parameterization = parameterization,
        prior.std.dev.nominal = 1,
        prior.range.nominal = NULL,
        prior.tau = NULL,
        prior.kappa = NULL,
        theta.prior.mean = NULL,
        theta.prior.prec = 0.1
      )
    }

    initial <- param$theta.prior.mean
  }

  if (include.nu) {
    initial <- c(initial, log(nu))
  }

  if (log.scale) {
    return(initial)
  } else {
    return(exp(initial))
  }
}




#' @name cut_decimals
#' @title Approximation function for covariance-based rSPDE models
#' @description Approximation function to be used to compute the
#' precision matrix for covariance-based rSPDE models
#' @param nu A real number
#' @return An approximation
#' @noRd

cut_decimals <- function(nu) {
  temp <- nu - floor(nu)
  if (temp < 10^(-3)) {
    temp <- 10^(-3)
  }
  if (temp > 0.999) {
    temp <- 0.999
  }
  return(temp)
}

#' @name check_class_inla_rspde
#' @title Check if the object inherits from inla_rspde class
#' @description Check if the object inherits from inla_rspde class
#' @param model A model to test if it inherits from inla_rspde
#' @return Gives an error if the object does not inherit from inla_rspde
#' @noRd

check_class_inla_rspde <- function(model) {
  if (!inherits(model, c("inla_rspde", "inla_rspde_matern1d", "inla_rspde_fintrinsic"))) {
    stop("You should provide a rSPDE model!")
  }
}

#' @name fem_mesh_order_1d
#' @title Get fem_mesh_matrices for 1d inla.mesh objects
#' @description Get fem_mesh_matrices for 1d inla.mesh objects
#' @param inla_mesh An INLA mesh
#' @param m_order the order of the FEM matrices
#' @return A list with fem_mesh_matrices
#' @noRd


fem_mesh_order_1d <- function(inla_mesh, m_order) {
  # fem_mesh <- rSPDE.fem1d(inla_mesh[["loc"]])
  # mesh_1d <- fmesher::fm_mesh_1d(inla_mesh[["loc"]])
  # fem_mesh <- fmesher::fm_fem(mesh_1d)
  mesh_1d <- fm_mesh_1d(inla_mesh[["loc"]])
  fem_mesh <- fm_fem(mesh_1d)
  C <- fem_mesh$c0
  C <- Matrix::Diagonal(dim(C)[1], rowSums(C))
  C <- as(C, "TsparseMatrix")
  G <- fem_mesh$g1
  Gk <- list()
  Ci <- C
  Ci@x <- 1 / (C@x)

  GCi <- G %*% Ci
  Gk[[1]] <- G
  # determine how many G_k matrices we want to create
  if (m_order > 1) {
    for (i in 2:m_order) {
      Gk[[i]] <- GCi %*% Gk[[i - 1]]
    }
  }

  # create a list contains all the finite element related matrices
  fem_mesh_matrices <- list()
  fem_mesh_matrices[["c0"]] <- C

  for (i in 1:m_order) {
    fem_mesh_matrices[[paste0("g", i)]] <- Gk[[i]]
  }
  return(fem_mesh_matrices)
}

#' @name generic_fem_mesh_order
#' @title Get fem_mesh_matrices from C and G matrices
#' @description Get fem_mesh_matrices from C and G matrices
#' @param fem_matrices A list with objects C and G
#' @param m_order the order of the FEM matrices
#' @return A list with fem_mesh_matrices
#' @noRd


generic_fem_mesh_order <- function(fem_matrices, m_order) {
  C <- fem_matrices$C
  C <- Matrix::Diagonal(dim(C)[1], rowSums(C))
  C <- INLA::inla.as.sparse(C)
  # C <- as(C,"TsparseMatrix")
  G <- fem_matrices$G
  Gk <- list()
  Ci <- C
  Ci@x <- 1 / (C@x)

  GCi <- G %*% Ci
  Gk[[1]] <- G
  # determine how many G_k matrices we want to create
  if (m_order > 1) {
    for (i in 2:m_order) {
      Gk[[i]] <- GCi %*% Gk[[i - 1]]
    }
  }

  # create a list contains all the finite element related matrices
  fem_mesh_matrices <- list()
  fem_mesh_matrices[["c0"]] <- C

  for (i in 1:m_order) {
    fem_mesh_matrices[[paste0("g", i)]] <- Gk[[i]]
  }
  return(fem_mesh_matrices)
}


#' @name get.sparsity.graph.rspde
#' @title Sparsity graph for rSPDE models
#' @description Creates the sparsity graph for rSPDE models
#' @param mesh An INLA mesh, optional
#' @param fem_mesh_matrices A list containing the FEM-related matrices.
#' The list should contain elements C, G, G_2, G_3, etc. Optional,
#' should be provided if mesh is not provided.
#' @param dim The dimension, optional. Should be provided if mesh
#' is not provided.
#' @param nu The smoothness parameter
#' @param force_non_integer Should nu be treated as non_integer?
#' @param rspde.order The order of the covariance-based rational SPDE approach.
#' @return The sparsity graph for rSPDE models to be used in R-INLA interface.
#' @noRd

get.sparsity.graph.rspde <- function(mesh = NULL,
                                     fem_mesh_matrices = NULL,
                                     nu,
                                     force_non_integer = FALSE,
                                     rspde.order = 2,
                                     dim = NULL) {
  if (!is.null(mesh)) {
    dim <- fmesher::fm_manifold_dim(mesh)
    if (!fmesher::fm_manifold(mesh, c("R1", "R2"))) {
      # FL: Is this actually required? Is fm_fem() etc support not sufficient?
      stop("The mesh should be from a flat manifold.")
    }
  } else if (is.null(dim)) {
    stop("If an INLA mesh is not provided, you should provide the dimension!")
  }
  sharp <- TRUE
  alpha <- nu + dim / 2

  m_alpha <- floor(alpha)

  integer_alpha <- (alpha %% 1 == 0)

  if (force_non_integer) {
    integer_alpha <- FALSE
  }

  if (!is.null(fem_mesh_matrices)) {
    if (integer_alpha) {
      return(fem_mesh_matrices[[paste0("g", m_alpha)]])
    } else {
      if (sharp) {
        if (m_alpha > 0) {
          return(bdiag(
            kronecker(
              diag(rep(1, rspde.order)),
              fem_mesh_matrices[[paste0("g", m_alpha + 1)]]
            ),
            fem_mesh_matrices[[paste0("g", m_alpha)]]
          ))
        } else {
          return(bdiag(
            kronecker(
              diag(rep(1, rspde.order)),
              fem_mesh_matrices[["g1"]]
            ),
            fem_mesh_matrices[["c0"]]
          ))
        }
      } else {
        return(kronecker(
          diag(rep(1, rspde.order + 1)),
          fem_mesh_matrices[[paste0("g", m_alpha + 1)]]
        ))
      }
    }
  } else if (!is.null(mesh)) {
    if (integer_alpha) {
      # fem_mesh_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha)
      # fem_mesh_matrices <- fmesher::fm_fem(mesh, order = m_alpha)
      fem_mesh_matrices <- fm_fem(mesh, order = m_alpha)
      return(fem_mesh_matrices[[paste0("g", m_alpha)]])
    } else {
      if (dim == 2) {
        # fem_mesh_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
        # fem_mesh_matrices <- fmesher::fm_fem(mesh, order = m_alpha + 1)
        fem_mesh_matrices <- fm_fem(mesh, order = m_alpha + 1)
      } else {
        fem_mesh_matrices <- fem_mesh_order_1d(mesh, m_order = m_alpha + 1)
      }


      if (sharp) {
        if (m_alpha > 0) {
          return(bdiag(
            kronecker(
              diag(rep(1, rspde.order)),
              fem_mesh_matrices[[paste0("g", m_alpha + 1)]]
            ),
            fem_mesh_matrices[[paste0("g", m_alpha)]]
          ))
        } else {
          return(bdiag(
            kronecker(
              diag(rep(1, rspde.order)),
              fem_mesh_matrices[["g1"]]
            ),
            fem_mesh_matrices[["c0"]]
          ))
        }
      } else {
        return(kronecker(
          diag(rep(1, rspde.order + 1)),
          fem_mesh_matrices[[paste0("g", m_alpha + 1)]]
        ))
      }
    }
  } else {
    stop("You should provide either mesh or fem_mesh_matrices!")
  }
}


#' @name build_sparse_matrix_rspde
#' @title Create sparse matrix from entries and graph
#' @description Create sparse matrix from entries and graph
#' @param entries The entries of the precision matrix
#' @param graph The sparsity graph of the precision matrix
#' @return index for rSPDE models.
#' @noRd

build_sparse_matrix_rspde <- function(entries, graph) {
  if (!is.null(graph)) {
    # graph <- as(graph, "dgTMatrix")
    graph <- as(graph, "TsparseMatrix")
    idx <- which(graph@i <= graph@j)
    Q <- Matrix::sparseMatrix(
      i = graph@i[idx], j = graph@j[idx], x = entries,
      symmetric = TRUE, index1 = FALSE
    )
  }
  return(Q)
}


#' @name analyze_sparsity_rspde
#' @title Analyze sparsity of matrices in the rSPDE approach
#' @description Auxiliar function to analyze sparsity of matrices
#' in the rSPDE approach
#' @param nu.upper.bound Upper bound for the smoothness parameter
#' @param dim The dimension of the domain
#' @param rspde.order The order of the rational approximation
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @param include_lower_order Logical. Should the lower-order terms
#' be included? They are needed for the cases
#' when alpha = nu + d/2 is integer or for when sharp is set to TRUE.
#' @param include_higher_order Logical. Should be included for when nu
#' is estimated or for when alpha = nu + d/2 is not an integer.
#' @return A list containing informations on sparsity of the precision matrices
#' @noRd

analyze_sparsity_rspde <- function(nu.upper.bound, dim, rspde.order,
                                   fem_mesh_matrices,
                                   include_lower_order = TRUE,
                                   include_higher_order = TRUE) {
  beta <- nu.upper.bound / 2 + dim / 4

  m_alpha <- floor(2 * beta)

  positions_matrices <- list()

  C_list <- symmetric_part_matrix(fem_mesh_matrices$c0)
  G_1_list <- symmetric_part_matrix(fem_mesh_matrices$g1)
  if (m_alpha < 2) {
    G_2_list <- symmetric_part_matrix(fem_mesh_matrices[["g2"]])
  }
  if (m_alpha > 1) {
    for (j in 2:(m_alpha)) {
      assign(
        paste0("G_", j, "_list"),
        symmetric_part_matrix(fem_mesh_matrices[[paste0("g", j)]])
      )
    }
  }

  if (include_higher_order) {
    assign(
      paste0("G_", m_alpha + 1, "_list"),
      symmetric_part_matrix(fem_mesh_matrices[[paste0(
        "g",
        m_alpha + 1
      )]])
    )

    positions_matrices[[1]] <- match(
      C_list$M,
      get(paste0("G_", m_alpha + 1, "_list"))[["M"]]
    )
  }

  idx_matrices <- list()

  idx_matrices[[1]] <- C_list$idx

  if (m_alpha > 0) {
    for (i in 1:m_alpha) {
      if (include_higher_order) {
        positions_matrices[[i + 1]] <- match(get(paste0(
          "G_", i,
          "_list"
        ))[["M"]], get(paste0(
          "G_", m_alpha + 1,
          "_list"
        ))[["M"]])
      }
      idx_matrices[[i + 1]] <- get(paste0("G_", i, "_list"))[["idx"]]
    }
  }

  if (include_higher_order) {
    idx_matrices[[m_alpha + 2]] <- get(paste0(
      "G_", m_alpha + 1,
      "_list"
    ))[["idx"]]
  }

  if (include_lower_order) {
    positions_matrices_less <- list()
    if (m_alpha > 0) {
      positions_matrices_less[[1]] <- match(C_list$M, get(paste0(
        "G_",
        m_alpha, "_list"
      ))[["M"]])
    } else {
      positions_matrices_less[[1]] <- match(C_list$M, get(paste0(
        "G_",
        1, "_list"
      ))[["M"]])
    }

    if (m_alpha > 1) {
      for (i in 1:(m_alpha - 1)) {
        positions_matrices_less[[i + 1]] <- match(get(paste0(
          "G_", i,
          "_list"
        ))[["M"]], get(paste0("G_", m_alpha, "_list"))[["M"]])
      }
    } else if (m_alpha == 1) {
      positions_matrices_less[[2]] <- seq_len(length(get(paste0(
        "G_",
        m_alpha, "_list"
      ))[["M"]]))
    }
  } else {
    positions_matrices_less <- NULL
  }

  return(list(
    positions_matrices = positions_matrices,
    idx_matrices = idx_matrices,
    positions_matrices_less = positions_matrices_less
  ))
}

#' @name symmetric_part_matrix
#' @title Gets the upper triangular part of a matrix
#' @description Gets the upper triangular part of a matrix
#' @param M A matrix or a sparse matrix
#' @return A sparse matrix formed by the upper triangular part of `M`.
#' @noRd

symmetric_part_matrix <- function(M) {
  # M <- as(M, "dgTMatrix")
  M <- as(M, "TsparseMatrix")
  idx <- which(M@i <= M@j)
  sM <- cbind(M@i[idx], M@j[idx])
  colnames(sM) <- NULL
  return(list(M = split(sM, seq(nrow(sM))), idx = idx))
}


#' @name get.roots
#' @title Get roots of the polynomials used in the operator based rational
#' approximation.
#' @description Get list with rational coefficients
#' @param order order of the rational approximation
#' @param beta value of beta to get the coefficients for.
#' @param type_interp Type of interpolation. Options are "linear" or "spline".
#' @return A list with coefficients.
#' @noRd
get.roots <- function(order, beta, type_interp = "linear") {
  if(!(order %in% c(1,2,3,4))) {
    stop("order must be one of the values 1,2,3,4.")
  }
  if (beta > 2) {
    beta <- beta - floor(beta - 1)
  }
  mt <- get(paste0("m", order, "table"))
  rb <- rep(0, order + 1)
  rc <- rep(0, order)
  if(type_interp == "linear"){
      if(order == 1) {
          rc = approx(mt$beta, mt[[paste0("rc")]], beta)$y
      } else {
          rc = sapply(1:order, function(i) {
              approx(mt$beta, mt[[paste0("rc.", i)]], beta)$y
          })
      }
      rb = sapply(1:(order+1), function(i) {
           approx(mt$beta, mt[[paste0("rb.", i)]], xout = beta)$y
      })
      factor = approx(mt$beta, mt$factor, xout = beta)$y
  } else if(type_interp == "spline") {
      if(order == 1) {
          rc = spline(mt$beta, mt[[paste0("rc")]], xout = beta)$y
      } else {
          rc = sapply(1:order, function(i) {
              spline(mt$beta, mt[[paste0("rc.", i)]], xout = beta)$y
          })
      }
      rb = sapply(1:(order+1), function(i) {
          spline(mt$beta, mt[[paste0("rb.", i)]], xout = beta)$y
      })
      factor = spline(mt$beta, mt$factor, xout = beta)$y
  } else {
      stop("invalid type. The options are 'linear' and 'spline'.")
  }
  return(list(rb = rb, rc = rc, factor = factor))
}

#' @name get_rational_coefficients
#' @title Get matrix with rational coefficients
#' @description Get matrix with rational coefficients
#' @param order order of the rational approximation
#' @param type_rational_approx Type of the rational
#' approximation. Options are "mix", "chebfun", "brasil", "chebfunLB" and "operator"
#' @return A matrix with rational approximations.
#' @noRd

get_rational_coefficients <- function(order, type_rational_approx) {
  if (type_rational_approx == "chebfun") {
    mt <- get(paste0("m", order, "t"))
  } else if (type_rational_approx == "brasil") {
    mt <- get(paste0("m_brasil", order, "t"))
  } else if (type_rational_approx == "chebfunLB") {
    mt <- get(paste0("m_chebfun", order, "t"))
  } else if(type_rational_approx == "mix"){
    mt_brasil <- get(paste0("m_brasil", order, "t"))
    mt_chebfun <- get(paste0("m", order, "t"))
    mt <- matrix(nrow = nrow(mt_brasil), ncol = ncol(mt_brasil))
    mt[1:500,] <- mt_brasil[1:500,]
    mt[501:999] <- mt_chebfun[501:999,]
  } else{
    stop("The options are 'mix', 'chebfun', 'brasil' and 'chebfunLB'!")
  }
  return(mt)
}


#' @name interp_rational_coefficients
#' @title Get list with interpolated rational coefficients
#' @description Get list with interpolated rational coefficients for specific
#' value of alpha.
#' @param order order of the rational approximation
#' @param type_rational_approx Type of the rational
#' approximation. Options are "chebfun", "brasil"
#' and "chebfunLB"
#' @param type_interp Type of interpolation. Options are "linear"
#' (linear interpolation), "log" (log-linear interpolation), "spline" (spline
#' interpolation) and "logspline" (log-spline interpolation).
#' @param alpha Value of alpha for the coefficients.
#' @return A list with rational approximations.
#' @noRd
interp_rational_coefficients <- function(order,
                                         type_rational_approx,
                                         type_interp = "spline",
                                         alpha){
    mt <- get_rational_coefficients(order = order,
                                    type_rational_approx=type_rational_approx)
    alpha <- cut_decimals(alpha)
    if(type_interp == "linear"){
        r = sapply(1:order, function(i) {
            approx(mt$alpha, mt[[paste0("r", i)]], alpha)$y
        })
        p = sapply(1:order, function(i) {
            approx(mt$alpha, mt[[paste0("p", i)]], alpha)$y
        })
        k = approx(mt$alpha, mt$k, cut_decimals(alpha))$y
    } else if (type_interp == "log"){
        r = sapply(1:order, function(i) {
            exp(approx(mt$alpha, log(mt[[paste0("r", i)]]), alpha)$y)
        })
        p = sapply(1:order, function(i) {
            -exp(approx(mt$alpha, log(-mt[[paste0("p", i)]]), alpha)$y)
        })
        k = exp(approx(mt$alpha, log(mt$k), alpha)$y)
    } else if(type_interp == "spline") {
        r = sapply(1:order, function(i) {
            spline(mt$alpha, mt[[paste0("r", i)]], xout = alpha)$y
        })
        p = sapply(1:order, function(i) {
            spline(mt$alpha, mt[[paste0("p", i)]], xout = alpha)$y
        })
        k = spline(mt$alpha, mt$k, xout = alpha)$y
    } else if(type_interp == "logspline") {
        r = sapply(1:order, function(i) {
            exp(spline(mt$alpha, log(mt[[paste0("r", i)]]), xout = alpha)$y)
        })
        p = sapply(1:order, function(i) {
            -exp(spline(mt$alpha, log(-mt[[paste0("p", i)]]), xout = alpha)$y)
        })
        k = exp(spline(mt$alpha, log(mt$k), xout = alpha)$y)
    } else {
        stop("invalid type. The options are 'linear', 'log', 'spline' and 'logspline'.")
    }
    return(list(k=k, r=r, p=p))
}

#' Changing the type of the rational approximation
#'
#' @param x A `CBrSPDE` or an `rpsde.inla` object
#' @param value The type of rational approximation.
#' The current options are "chebfun", "brasil" and "chebfunLB"
#'
#' @return An object of the same class with the new rational approximation.
#' @export
#'
`rational.type<-` <- function(x, value) {
  object <- x

  type_rational_approximation <- value
  type_rational_approximation <- type_rational_approximation[[1]]
  if (!(type_rational_approximation %in% c("chebfun", "brasil", "chebfunLB"))) {
    stop('The possible types are "chebfun", "brasil" and "chebfunLB"!')
  }
  if (inherits(object, "CBrSPDEobj")) {
    model <- update(x, type_rational_approximation = value)
  } else if (inherits(object, "inla_rspde")) {
    nu.upper.bound <- object$nu.upper.bound
    prior.nu.dist <- object$prior.nu.dist
    mesh <- object$mesh
    nu <- object[["nu"]]
    rspde.order <- object$rspde.order
    parameterization <- object$parameterization
    theta.prior.prec <- object$theta.prior.prec
    theta.prior.mean <- object$theta.prior.mean
    start.theta <- object$start.theta
    prior.nu <- object$prior.nu
    start.nu <- object$start.nu
    debug <- object$debug


    model <- rspde.matern(mesh,
      nu.upper.bound = nu.upper.bound,
      rspde.order = rspde.order,
      nu = nu,
      debug = debug,
      parameterization = parameterization,
      theta.prior.mean = theta.prior.mean,
      theta.prior.prec = theta.prior.prec,
      start.theta = start.theta,
      prior.nu = prior.nu,
      start.nu = start.nu,
      prior.nu.dist = prior.nu.dist,
      type.rational.approx = type_rational_approximation
    )
  } else {
    stop("The object must be of class 'CBrSPDE' or 'inla_rspde'!")
  }
  return(model)
}



#' Get type of rational approximation.
#'
#' @param object A `CBrSPDEobj` object or an `inla_rspde` object.
#'
#' @return The type of rational approximation.
#' @export
#'
rational.type <- function(object) {
  if (inherits(object, "CBrSPDEobj")) {
    return(object$type_rational_approximation)
  } else if (inherits(object, "inla_rspde")) {
    return(object$type.rational.approx)
  } else if (inherits(object, "rSPDEobj")) {
    return("chebfun")
  } else {
    stop("Not a valid rSPDE object!")
  }
}


#' Changing the order of the rational approximation
#'
#' @param x A `CBrSPDE` or an `rpsde.inla` object
#' @param value The order of rational approximation.
#'
#' @return An object of the same class with the new order
#' of rational approximation.
#' @export
#'
`rational.order<-` <- function(x, value) {
  object <- x

  rspde.order <- value
  rspde.order <- rspde.order[[1]]

  if (inherits(object, "CBrSPDEobj") || inherits(object, "rSPDEobj")) {
    model <- update(object, m = rspde.order)
  } else if (inherits(object, "inla_rspde")) {
    if (rspde.order > 0 && object$integer.nu) {
      warning("The order was not changed since there is no
      rational approximation (an integer model was
      considered).")
      return(object)
    }
    nu.upper.bound <- object$nu.upper.bound
    prior.nu.dist <- object$prior.nu.dist
    mesh <- object$mesh
    nu <- object[["nu"]]
    parameterization <- object$parameterization
    theta.prior.prec <- object$theta.prior.prec
    theta.prior.mean <- object$theta.prior.mean
    start.theta <- object$start.theta
    prior.nu <- object$prior.nu
    start.nu <- object$start.nu
    type_rational_approximation <- object$type.rational.approx
    debug <- object$debug


    model <- rspde.matern(mesh,
      nu.upper.bound = nu.upper.bound,
      rspde.order = rspde.order,
      nu = nu,
      debug = debug,
      parameterization = parameterization,
      theta.prior.mean = theta.prior.mean,
      theta.prior.prec = theta.prior.prec,
      start.theta = start.theta,
      prior.nu = prior.nu,
      start.nu = start.nu,
      prior.nu.dist = prior.nu.dist,
      type.rational.approx = type_rational_approximation
    )
  } else if (!is.null(attr(object, "inla_rspde_Amatrix"))) {
    n_temp <- ncol(object)
    old_rspde.order <- attr(object, "rspde.order")
    orig_dim <- n_temp / (old_rspde.order + 1)
    A <- object[, 1:orig_dim]
    Abar <- kronecker(matrix(1, 1, rspde.order + 1), A)
    attr(Abar, "inla_rspde_Amatrix") <- TRUE
    attr(Abar, "rspde.order") <- rspde.order
    integer_nu <- attr(object, "integer_nu")
    if (integer_nu && rspde.order > 0) {
      warning("The order was not changed since there is
      no rational approximation (an integer model was
      considered).")
      return(object)
    }
    attr(Abar, "integer_nu") <- integer_nu
    return(Abar)
  } else if (inherits(object, "inla_rspde_index")) {
    integer_nu <- attr(object, "integer_nu")

    if (integer_nu && rspde.order > 0) {
      warning("The order was not changed since there is
      no rational approximation (an integer model was
      considered).")
      return(object)
    }

    n_mesh <- attr(object, "n.mesh")
    name <- attr(object, "name")
    n.group <- attr(object, "n.group")
    n.repl <- attr(object, "n.repl")

    factor_rspde <- rspde.order + 1

    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")

    out <- list()
    out[[name]] <- as.vector(sapply(1:factor_rspde, function(i) {
      rep(rep(((i - 1) * n_mesh + 1):(i * n_mesh),
        times = n.group
      ), times = n.repl)
    }))
    out[[name.group]] <- rep(rep(rep(1:n.group, each = n_mesh),
      times = n.repl
    ), times = factor_rspde)
    out[[name.repl]] <- rep(rep(1:n.repl, each = n_mesh * n.group),
      times = factor_rspde
    )
    class(out) <- c("inla_rspde_index", class(out))
    attr(out, "rspde.order") <- rspde.order
    attr(out, "integer_nu") <- integer_nu
    attr(out, "n.mesh") <- n_mesh
    attr(out, "name") <- name
    attr(out, "n.group") <- n.group
    attr(out, "n.repl") <- n.repl
    return(out)
  } else {
    stop("The object must be of class 'CBrSPDE' or 'inla_rspde'!")
  }
  return(model)
}


#' Get the order of rational approximation.
#'
#' @param object A `CBrSPDEobj` object or an `inla_rspde` object.
#'
#' @return The order of rational approximation.
#' @export
#'
rational.order <- function(object) {
  if (inherits(object, "CBrSPDEobj") || inherits(object, "rSPDEobj")) {
    return(object$m)
  } else if (inherits(object, "inla_rspde")) {
    return(object$rspde.order)
  } else if (!is.null(attr(object, "inla_rspde_Amatrix"))) {
    return(attr(object, "rspde.order"))
  } else if (inherits(object, "inla_rspde_index")) {
    return(attr(object, "rspde.order"))
  } else {
    stop("Not a valid rSPDE object!")
  }
}


#' Check user input.
#'
#' @param param A parameter to validate.
#' @param label Label for the parameter (used in error messages).
#' @param lower_bound Optional lower bound for the parameter.
#' @param dim Expected dimension of the parameter (default is 1 for scalar).
#' @param upper_bound Optional upper bound for the parameter.
#'
#' @return The validated parameter.
#' @noRd
#'
rspde_check_user_input <- function(param, label, lower_bound = NULL, dim = 1, upper_bound = NULL) {
  if (!is.numeric(param)) {
    stop(paste(label, "should be a numeric value!"))
  }
  
  if (length(param) != dim) {
    if (dim == 1) {
      stop(paste(label, "should be a single numeric value!"))
    } else {
      stop(paste(label, "should have a length of", dim, "!"))
    }
  }
  
  if (!is.null(lower_bound) && any(param < lower_bound)) {
    stop(paste(label, "should be greater than or equal to", lower_bound, "!"))
  }
  
  if (!is.null(upper_bound) && any(param > upper_bound)) {
    stop(paste(label, "should be less than or equal to", upper_bound, "!"))
  }
  
  return(param)
}


#' Process inputs likelihood
#'
#' @param kappa kappa
#' @param tau tau
#' @param nu nu
#' @param sigma.e sigma.e
#'
#' @return List with the positions
#' @noRd

likelihood_process_inputs_spde <- function(kappa, tau, nu, sigma.e) {
  param_vector <- c("tau", "kappa", "nu", "sigma.e")
  if (!is.null(tau)) {
    param_vector <- setdiff(param_vector, "tau")
  }
  if (!is.null(kappa)) {
    param_vector <- setdiff(param_vector, "kappa")
  }
  if (!is.null(nu)) {
    param_vector <- setdiff(param_vector, "nu")
  }
  if (!is.null(sigma.e)) {
    param_vector <- setdiff(param_vector, "sigma.e")
  }
  if (length(param_vector) == 0) {
    stop("You should leave at least one parameter free.")
  }
  return(param_vector)
}

#' Process inputs likelihood
#'
#' @param kappa kappa
#' @param tau tau
#' @param nu nu
#' @param sigma.e sigma.e
#'
#' @return List with the positions
#' @noRd

likelihood_process_inputs_matern <- function(range, sigma, nu, sigma.e) {
  param_vector <- c("sigma", "range", "nu", "sigma.e")
  if (!is.null(sigma)) {
    param_vector <- setdiff(param_vector, "sigma")
  }
  if (!is.null(range)) {
    param_vector <- setdiff(param_vector, "range")
  }
  if (!is.null(nu)) {
    param_vector <- setdiff(param_vector, "nu")
  }
  if (!is.null(sigma.e)) {
    param_vector <- setdiff(param_vector, "sigma.e")
  }
  if (length(param_vector) == 0) {
    stop("You should leave at least one parameter free.")
  }
  return(param_vector)
}

#' Process parameters likelihood
#'
#' @param theta vector of parameters
#' @param param_vector vector of parameters to be used
#' @param which_par which parameter to consider
#' @param logscale log scale?
#'
#' @return The value in the correct scale
#' @noRd

likelihood_process_parameters <- function(theta, param_vector, which_par, logscale) {
  coord_par <- which(which_par == param_vector)
  if (logscale) {
    param_value <- exp(theta[[coord_par]])
  } else {
    param_value <- theta[[coord_par]]
  }
  return(param_value)
}


#' @noRd
# Get priors and starting values
# Based on INLA::param2.matern.orig()

get_parameters_rSPDE <- function(
    mesh, alpha,
    B.tau,
    B.kappa,
    B.sigma,
    B.range,
    nu.nominal,
    alpha.nominal,
    parameterization,
    prior.std.dev.nominal,
    prior.range.nominal,
    prior.tau,
    prior.kappa,
    theta.prior.mean,
    theta.prior.prec,
    mesh.range = NULL,
    d = NULL,
    n.spde = NULL) {
  if (!is.null(mesh)) {
    if (!inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d"))) {
      stop("The mesh should be created using fmesher!")
    }

    d <- fmesher::fm_manifold_dim(mesh)
    n.spde <- fmesher::fm_dof(mesh)
  } else {
    if (is.null(d)) {
      stop("If you do not provide the mesh, you must provide the dimension!")
    }
    if (is.null(n.spde)) {
      stop("If you do not provide the mesh, you must provide n.spde!")
    }
  }

  if (is.null(B.tau) && is.null(B.sigma)) {
    stop("One of B.tau or B.sigma must not be NULL.")
  }
  if (is.null(B.kappa) && is.null(B.range)) {
    stop("One of B.kappa or B.range must not be NULL.")
  }



  if (parameterization == "spde") {
    n.theta <- ncol(B.kappa) - 1L

    B.kappa <- prepare_B_matrices(
      B.kappa, n.spde,
      n.theta
    )
    B.tau <- prepare_B_matrices(B.tau, n.spde, n.theta)
  } else if (parameterization == "matern") {
    n.theta <- ncol(B.sigma) - 1L

    B.sigma <- prepare_B_matrices(
      B.sigma, n.spde,
      n.theta
    )
    B.range <- prepare_B_matrices(
      B.range, n.spde,
      n.theta
    )

    B.kappa <- cbind(
      0.5 * log(8 * nu.nominal) - B.range[, 1],
      -B.range[, -1, drop = FALSE]
    )

    B.tau <- cbind(
      0.5 * (lgamma(nu.nominal) - lgamma(alpha.nominal) -
        d / 2 * log(4 * pi)) - nu.nominal * B.kappa[, 1] -
        B.sigma[, 1],
      -nu.nominal * B.kappa[, -1, drop = FALSE] -
        B.sigma[, -1, drop = FALSE]
    )
  } else if (parameterization == "matern2") {
    n.theta <- ncol(B.sigma) - 1L

    B.sigma <- prepare_B_matrices(
      B.sigma, n.spde,
      n.theta
    )
    B.range <- prepare_B_matrices(
      B.range, n.spde,
      n.theta
    )

    B.kappa <- -B.range

    B.tau <- cbind(
      0.5 * (lgamma(nu.nominal) - lgamma(alpha.nominal) -
        d / 2 * log(4 * pi)) - nu.nominal * B.kappa[, 1] -
        0.5 * B.sigma[, 1],
      -nu.nominal * B.kappa[, -1, drop = FALSE] -
        0.5 * B.sigma[, -1, drop = FALSE]
    )
  }


  if (is.null(theta.prior.prec)) {
    theta.prior.prec <- diag(0.1, n.theta, n.theta)
  } else {
    theta.prior.prec <- as.matrix(theta.prior.prec)
    if (ncol(theta.prior.prec) == 1) {
      theta.prior.prec <- diag(
        as.vector(theta.prior.prec),
        n.theta, n.theta
      )
    }
    if ((nrow(theta.prior.prec) != n.theta) || (ncol(theta.prior.prec) !=
      n.theta)) {
      stop(paste(
        "Size of theta.prior.prec is (", paste(dim(theta.prior.prec),
          collapse = ",", sep = ""
        ), ") but should be (",
        paste(c(n.theta, n.theta), collapse = ",", sep = ""),
        ")."
      ))
    }
  }


  if (is.null(theta.prior.mean)) {
    if (is.null(prior.range.nominal)) {
      if (is.null(mesh.range)) {
        mesh.range <- ifelse(d == 2, (max(c(diff(range(mesh$loc[
          ,
          1
        ])), diff(range(mesh$loc[, 2])), diff(range(mesh$loc[
          ,
          3
        ]))))), diff(mesh$interval))
      }
      prior.range.nominal <- mesh.range * 0.2
    }
    if (is.null(prior.kappa)) {
      prior.kappa <- sqrt(8 * nu.nominal) / prior.range.nominal
    }
    if (is.null(prior.tau)) {
      prior.tau <- sqrt(gamma(nu.nominal) / gamma(alpha.nominal) / ((4 *
        pi)^(d / 2) * prior.kappa^(2 * nu.nominal) * prior.std.dev.nominal^2))
    }
    if (n.theta > 0) {
      if (parameterization == "spde") {
        theta.prior.mean <- qr.solve(rbind(
          B.tau[, -1, drop = FALSE],
          B.kappa[, -1, drop = FALSE]
        ), c(log(prior.tau) -
          B.tau[, 1], log(prior.kappa) - B.kappa[, 1]))
      } else if (parameterization == "matern") {
        theta.prior.mean <- qr.solve(rbind(
          B.sigma[, -1, drop = FALSE],
          B.range[, -1, drop = FALSE]
        ), c(log(prior.std.dev.nominal) -
          B.sigma[, 1], log(prior.range.nominal) - B.range[, 1]))
      } else if (parameterization == "matern2") {
        theta.prior.mean <- qr.solve(rbind(
          B.sigma[, -1, drop = FALSE],
          B.range[, -1, drop = FALSE]
        ), c(2 * log(prior.std.dev.nominal) -
          B.sigma[, 1], -log(prior.kappa) - B.range[, 1]))
      }
    } else {
      theta.prior.mean <- rep(0, n.theta)
    }
  }
  param <- list(
    B.tau = B.tau,
    B.kappa = B.kappa, theta.prior.mean = theta.prior.mean,
    theta.prior.prec = theta.prior.prec
  )
  return(param)
}

#' @noRd
# Check B matrices and adjust the number of lines
# Based on INLA:::inla.spde.homogenise_B_matrix()

prepare_B_matrices <- function(B, n.spde, n.theta) {
  if (!is.numeric(B)) {
    stop("B matrix must be numeric.")
  }
  if (is.matrix(B)) {
    if ((nrow(B) != 1) && (nrow(B) != n.spde)) {
      stop(paste("B matrix must have either 1 or", as.character(n.spde), "rows."))
    }
    if ((ncol(B) != 1) && (ncol(B) != 1 + n.theta)) {
      stop(paste("B matrix must have 1 or", as.character(1 +
        n.theta), "columns."))
    }
    if (ncol(B) == 1) {
      return(cbind(as.vector(B), matrix(0, n.spde, n.theta)))
    } else if (ncol(B) == 1 + n.theta) {
      if (nrow(B) == 1) {
        return(matrix(as.vector(B), n.spde, 1 + n.theta,
          byrow = TRUE
        ))
      } else if (nrow(B) == n.spde) {
        return(B)
      }
    }
  } else {
    if ((length(B) == 1) || (length(B) == n.spde)) {
      return(cbind(B, matrix(0, n.spde, n.theta)))
    } else if (length(B) == 1 + n.theta) {
      return(matrix(B, n.spde, 1 + n.theta, byrow = TRUE))
    } else {
      stop(paste(
        "Length of B must be 1,", as.character(1 + n.theta),
        "or", as.character(n.spde)
      ))
    }
  }
  stop("Unrecognised structure for B matrix")
}



#' @noRd
# Get priors and starting values
# Based on INLA::param2.matern.orig()

get_parameters_rSPDE_graph <- function(
    graph_obj, alpha,
    B.tau,
    B.kappa,
    B.sigma,
    B.range,
    nu.nominal,
    alpha.nominal,
    parameterization,
    prior.std.dev.nominal,
    prior.range.nominal,
    prior.tau,
    prior.kappa,
    theta.prior.mean,
    theta.prior.prec) {
  if (!inherits(graph_obj, "metric_graph")) {
    stop("The graph object should be of class metric_graph!")
  }
  if (is.null(B.tau) && is.null(B.sigma)) {
    stop("One of B.tau or B.sigma must not be NULL.")
  }
  if (is.null(B.kappa) && is.null(B.range)) {
    stop("One of B.kappa or B.range must not be NULL.")
  }

  d <- 1
  n.spde <- nrow(graph_obj$mesh$C)

  if (parameterization == "spde") {
    n.theta <- ncol(B.kappa) - 1L
    B.kappa <- prepare_B_matrices(
      B.kappa, n.spde,
      n.theta
    )
    B.tau <- prepare_B_matrices(B.tau, n.spde, n.theta)
  } else if (parameterization == "matern") {
    n.theta <- ncol(B.sigma) - 1L
    B.sigma <- prepare_B_matrices(
      B.sigma, n.spde,
      n.theta
    )
    B.range <- prepare_B_matrices(
      B.range, n.spde,
      n.theta
    )

    B.kappa <- cbind(
      0.5 * log(8 * nu.nominal) - B.range[, 1],
      -B.range[, -1, drop = FALSE]
    )

    B.tau <- cbind(0.5 * (lgamma(nu.nominal) - lgamma(alpha.nominal) -
      d / 2 * log(4 * pi)) - nu.nominal * B.kappa[, 1] -
      B.sigma[, 1], -nu.nominal * B.kappa[, -1, drop = FALSE] -
      B.sigma[, -1, drop = FALSE])
  }


  if (is.null(theta.prior.prec)) {
    theta.prior.prec <- diag(0.1, n.theta, n.theta)
  } else {
    theta.prior.prec <- as.matrix(theta.prior.prec)
    if (ncol(theta.prior.prec) == 1) {
      theta.prior.prec <- diag(
        as.vector(theta.prior.prec),
        n.theta, n.theta
      )
    }
    if ((nrow(theta.prior.prec) != n.theta) || (ncol(theta.prior.prec) !=
      n.theta)) {
      stop(paste(
        "Size of theta.prior.prec is (", paste(dim(theta.prior.prec),
          collapse = ",", sep = ""
        ), ") but should be (",
        paste(c(n.theta, n.theta), collapse = ",", sep = ""),
        ")."
      ))
    }
  }


  if (is.null(theta.prior.mean)) {
    if (is.null(prior.range.nominal)) {
      if (is.null(graph_obj$geo_dist)) {
        graph_obj$compute_geodist(obs = FALSE)
      } else if (is.null(graph_obj$geo_dist[[".vertices"]])) {
        graph_obj$compute_geodist(obs = FALSE)
      }
      finite_geodist <- is.finite(graph_obj$geo_dist[[".vertices"]])
      finite_geodist <- graph_obj$geo_dist[[".vertices"]][finite_geodist]
      prior.range.nominal <- max(finite_geodist) * 0.2
    }
    if (is.null(prior.kappa)) {
      prior.kappa <- sqrt(8 * nu.nominal) / prior.range.nominal
    }
    if (is.null(prior.tau)) {
      prior.tau <- sqrt(gamma(nu.nominal) / gamma(alpha.nominal) / (4 *
        pi * prior.kappa^(2 * nu.nominal) * prior.std.dev.nominal^2))
    }
    if (n.theta > 0) {
      if (parameterization == "spde") {
        theta.prior.mean <- qr.solve(rbind(
          B.tau[, -1, drop = FALSE],
          B.kappa[, -1, drop = FALSE]
        ), c(log(prior.tau) -
          B.tau[, 1], log(prior.kappa) - B.kappa[, 1]))
      } else if (parameterization == "matern") {
        theta.prior.mean <- qr.solve(rbind(
          B.sigma[, -1, drop = FALSE],
          B.range[, -1, drop = FALSE]
        ), c(log(prior.std.dev.nominal) -
          B.sigma[, 1], log(prior.range.nominal) - B.range[, 1]))
      }
    } else {
      theta.prior.mean <- rep(0, n.theta)
    }
  }
  param <- list(
    B.tau = B.tau,
    B.kappa = B.kappa, theta.prior.mean = theta.prior.mean,
    theta.prior.prec = theta.prior.prec
  )
  return(param)
}



#' @noRd

# Function to convert B.sigma and B.range to B.tau and B.kappa

convert_B_matrices <- function(B.sigma, B.range, n.spde, nu.nominal, d) {
  n.theta <- ncol(B.sigma) - 1L

  alpha.nominal <- nu.nominal + d / 2

  B.sigma <- prepare_B_matrices(
    B.sigma, n.spde,
    n.theta
  )
  B.range <- prepare_B_matrices(
    B.range, n.spde,
    n.theta
  )

  B.kappa <- cbind(
    0.5 * log(8 * nu.nominal) - B.range[, 1],
    -B.range[, -1, drop = FALSE]
  )

  B.tau <- cbind(
    0.5 * (lgamma(nu.nominal) - lgamma(alpha.nominal) -
      d / 2 * log(4 * pi)) - nu.nominal * B.kappa[, 1] -
      B.sigma[, 1],
    -nu.nominal * B.kappa[, -1, drop = FALSE] -
      B.sigma[, -1, drop = FALSE]
  )

  return(list(B.tau = B.tau, B.kappa = B.kappa))
}

#' @noRd
# Change parameterization in rspde_lme to matern

change_parameterization_lme <- function(likelihood, d, nu, par, hessian, 
                                        fixed_params = c(tau = FALSE, kappa = FALSE)) {
  tau <- par[1]
  kappa <- par[2]

  C1 <- sqrt(8 * nu)
  C2 <- sqrt(gamma(nu) / ((4 * pi)^(d / 2) * gamma(nu + d / 2)))

  sigma <- C2 / (tau * kappa^nu)
  range <- C1 / kappa

  # Initialize result vectors
  coeff <- c(sigma, range)
  std_random <- rep(NA, 2)  # Default to NA for fixed parameters

  # If both parameters are fixed, we are done - return NAs for std errors
  if (all(fixed_params)) {
    return(list(coeff = coeff, std_random = std_random))
  }

  grad_par <- matrix(c(
    -C2 / (kappa^nu * sigma^2), 0,
    nu * range^(nu - 1) * C2 / (sigma * C1^nu),
    -C1 / range^2
  ), nrow = 2, ncol = 2)

  grad_par <- grad_par[!fixed_params, , drop=FALSE]

  new_observed_fisher <- t(grad_par) %*% hessian %*% (grad_par)

  # No need to include the additional term as the gradient is approximately zero.
  # from some numerical experiments, the approximation without the additional term
  # seems to be better in general.

  # hess_par <- matrix(c(2*C2/(kappa^nu * sigma^3), 0,
  #                     -nu * C2/((sigma^2) * (C1^nu)) * range^(nu-1),
  #                     2*C1/range^3) , ncol=2, nrow=2)

  # if(!improve_gradient){
  #   grad_lik <- numDeriv::grad(likelihood, log(par), method = "simple", method.args = gradient_args)
  # } else{
  #   grad_lik <- numDeriv::grad(likelihood, log(par), method = "Richardson", method.args = gradient_args)
  # }

  # grad_lik <- c(1/tau, 1/kappa) * grad_lik

  # add_mat <- diag(grad_lik) %*% hess_par

  # add_mat <- 0.5 * (add_mat + t(add_mat))

  # new_observed_fisher <- new_observed_fisher + add_mat

  inv_fisher <- tryCatch(solve(new_observed_fisher), error = function(e) matrix(NA, nrow(new_observed_fisher), ncol(new_observed_fisher)))

  if (!any(is.na(inv_fisher))) {
    std_random <- sqrt(diag(inv_fisher))
    
    # Set standard errors to NA for fixed parameters
    if (fixed_params[1]) std_random[1] <- NA  # sigma
    if (fixed_params[2]) std_random[2] <- NA  # range
  }

  # new_lik <- function(theta){
  #       sigma <- exp(theta[1])
  #       range <- exp(theta[2])

  #       kappa <- C1/range
  #       tau <- C2/(sigma * kappa^nu)
  #       return(likelihood(log(c(tau,kappa))))
  # }

  # hess_tmp <- numDeriv::hessian(new_lik, log(c(sigma,range)))

  # hess_tmp <- diag(c(1/sigma, 1/range)) %*% hess_tmp %*% diag(c(1/sigma, 1/range))

  return(list(coeff = coeff, std_random = std_random))
}



#' @noRd
#'

return_same_input_type_matrix_vector <- function(v, orig_v) {
  if (isS4(orig_v)) {
    return(v)
  } else {
    v_out <- as.matrix(v)
    dim(v_out) <- dim(orig_v)
    return(v_out)
  }
}



#' find indices of the rows with all NA's in lists
#' @noRd
#'
idx_not_all_NA <- function(data_list) {
  data_list[[".edge_number"]] <- NULL
  data_list[[".distance_on_edge"]] <- NULL
  data_list[[".coord_x"]] <- NULL
  data_list[[".coord_y"]] <- NULL
  data_list[[".group"]] <- NULL
  data_names <- names(data_list)
  n_data <- length(data_list[[data_names[1]]])
  idx_non_na <- logical(n_data)
  for (i in 1:n_data) {
    na_idx <- lapply(data_list, function(dat) {
      return(is.na(dat[i]))
    })
    idx_non_na[i] <- !all(unlist(na_idx))
  }
  return(idx_non_na)
}

#' find indices of the rows with at least one NA's in lists
#' @noRd
#'
idx_not_any_NA <- function(data_list) {
  data_list[[".edge_number"]] <- NULL
  data_list[[".distance_on_edge"]] <- NULL
  data_list[[".coord_x"]] <- NULL
  data_list[[".coord_y"]] <- NULL
  data_list[[".group"]] <- NULL
  data_names <- names(data_list)
  n_data <- length(data_list[[data_names[1]]])
  idx_non_na <- logical(n_data)
  for (i in 1:n_data) {
    na_idx <- lapply(data_list, function(dat) {
      return(is.na(dat[i]))
    })
    idx_non_na[i] <- !any(unlist(na_idx))
  }
  return(idx_non_na)
}


#' @noRd
#'

select_indexes <- function(data, idx) {
  if (inherits(data, "SpatialPointsDataFrame")) {
    data <- data[idx, , drop = FALSE]
  } else {
    data <- lapply(data, function(dat) {
      if (is.null(dim(dat))) {
        return(dat[idx])
      } else {
        return(dat[idx, , drop = FALSE])
      }
    })
  }
  return(data)
}


#' Create train and test splits for cross-validation
#'
#' @description
#' Creates train and test splits for cross-validation by handling multiple data types
#' and supporting k-fold, leave-one-out (LOO), and leave-percentage-out (LPO) methods.
#' Handles missing values and maintains data structure across multiple datasets.
#'
#' @param data_list A list of datasets, one per likelihood. Each dataset can be a data.frame, 
#'        SpatialPointsDataFrame, or metric_graph_data object
#' @param cv_type Type of cross-validation: "k-fold", "loo", or "lpo". Default is "k-fold"
#' @param k Number of folds for k-fold CV. Default is 5
#' @param percentage Training data percentage for LPO CV (1-99). Default is 20
#' @param number_folds Number of folds for LPO CV. Default is 10
#'
#' @return A list where each element contains:
#'   \item{train}{Indices for training data mapped to original datasets}
#'   \item{test}{Indices for test data mapped to original datasets}
#'
#' @details
#' The function handles NA values by removing rows with any missing values before
#' creating splits. For multiple datasets, indices are mapped back to their original
#' positions in each dataset.
#' @export 

create_train_test_indices <- function(data_list, cv_type = c("k-fold", "loo", "lpo"),
                                    k = 5, percentage = 20, number_folds = 10) {
  # First concatenate all data
  if (inherits(data_list[[1]], "metric_graph_data")) {
    data_list <- lapply(data_list, as.data.frame)
  } 
  
  data <- do.call(rbind, data_list)
  
  # Get indices for concatenated data as before
  idx <- seq_len(nrow(data))
    
  # Get cumulative sizes to map back to individual datasets
  n_samples <- sapply(data_list, nrow)
  cum_sizes <- cumsum(c(0, n_samples))
  
  # Function to map concatenated indices to individual dataset indices
  map_to_likelihood_indices <- function(indices) {
    lapply(seq_along(data_list), function(i) {
      likelihood_indices <- indices[indices > cum_sizes[i] & indices <= cum_sizes[i + 1]]
      likelihood_indices - cum_sizes[i]
    })
  }
  
  if (cv_type == "k-fold") {
    folds <- cut(sample(idx), breaks = k, label = FALSE)
    fold_list <- lapply(1:k, function(i) {
      test_idx <- which(folds == i, arr.ind = TRUE)
      train_idx <- idx[-test_idx]
      test_idx <- idx[test_idx]
      
      list(
        train = map_to_likelihood_indices(train_idx),
        test = map_to_likelihood_indices(test_idx)
      )
    })
  } else if (cv_type == "loo") {
    fold_list <- lapply(seq_along(idx), function(i) {
      list(
        train = map_to_likelihood_indices(idx[-i]),
        test = map_to_likelihood_indices(idx[i])
      )
    })
  } else if (cv_type == "lpo") {
    fold_list <- lapply(1:number_folds, function(i) {
      test_idx <- sample(idx, size = (1 - percentage / 100) * length(idx))
      train_idx <- idx[-match(test_idx, idx)]
      
      list(
        train = map_to_likelihood_indices(train_idx),
        test = map_to_likelihood_indices(test_idx)
      )
    })
  }
  
  return(fold_list)
}

# Check for required packages
#' @noRd
check_packages <- function(packages, func) {
    are_installed <-vapply(packages,
                           function(x) {
                               requireNamespace(x, quietly = TRUE)
                               },
                           TRUE
        )
    if (any(!are_installed)) {
        stop(paste0("Needed package(s) ",
                    paste0("'", packages[!are_installed], "'", collapse = ", "),
                    " not installed, but are needed by ", func)
             )
    }
}

#' @noRd
# Get appropriate shared library
get_shared_library <- function(shared_lib) {
  if (shared_lib == "INLA") {
    return(INLA::inla.external.lib("rSPDE"))
  }
  if (shared_lib == "rSPDE") {
    rspde_lib <- system.file("shared", package = "rSPDE")
    return(ifelse(Sys.info()["sysname"] == "Windows",
                 paste0(rspde_lib, "/rspde_cgeneric_models.dll"),
                 paste0(rspde_lib, "/rspde_cgeneric_models.so")))
  }
  if (shared_lib == "detect") {
    rspde_lib_local <- system.file("shared", package = "rSPDE")
    lib_path <- ifelse(Sys.info()["sysname"] == "Windows",
                      paste0(rspde_lib_local, "/rspde_cgeneric_models.dll"),
                      paste0(rspde_lib_local, "/rspde_cgeneric_models.so"))
    return(if (file.exists(lib_path)) lib_path else INLA::inla.external.lib("rSPDE"))
  }
  stop("'shared_lib' must be 'INLA', 'rSPDE', or 'detect'")
}

#' @noRd
set_prior <- function(prior, default_mean, default_precision, p = 1) {
  # Validate default parameters
  if (!is.numeric(default_mean) || length(default_mean) != p) {
    stop(paste("default_mean must be a numeric vector of length equal to",p,"."))
  }
  if (!is.numeric(default_precision) || length(default_precision) != p || any(default_precision <= 0)) {
    stop(paste("default_precision must be a positive numeric vector of length equal to",p,"."))
  }

  # Return default prior if none is provided
  if (is.null(prior)) {
    return(list(mean = default_mean, precision = default_precision))
  }

  # Ensure prior only contains allowed elements
  allowed_elements <- c("mean", "precision")
  invalid_elements <- setdiff(names(prior), allowed_elements)
  if (length(invalid_elements) > 0) {
    warning(sprintf("Invalid elements in prior: %s. Only 'mean' and 'precision' are allowed.",
                    paste(invalid_elements, collapse = ", ")))
  }

  # Validate and set 'mean'
  if (!is.null(prior$mean)) {
    if (!is.numeric(prior$mean) || length(prior$mean) != p) {
      stop(sprintf("'mean' must be a numeric vector of length %d.", p))
    }
  } else {
    prior$mean <- default_mean  # Use default mean if not provided
  }

  # Validate and set 'precision'
  if (!is.null(prior$precision)) {
    if (!is.numeric(prior$precision) || length(prior$precision) != p || any(prior$precision <= 0)) {
      stop(sprintf("'precision' must be a positive numeric vector of length %d.", p))
    }
  } else {
    prior$precision <- default_precision  # Use default precision if not provided
  }

  return(prior)
}

handle_prior_nu <- function(prior.nu, nu.upper.bound, nu.prec.inc = 0.01, prior.nu.dist = "lognormal") {
  if (is.null(prior.nu)) {
    prior.nu <- list()
  }
  
  # Check and set loglocation
  if (is.null(prior.nu$loglocation)) {
    prior.nu$loglocation <- log(min(1, nu.upper.bound / 2))
  } else if (length(prior.nu$loglocation) != 1) {
    warning("'prior.nu$loglocation' has length > 1. Only the first element will be used.")
    prior.nu$loglocation <- prior.nu$loglocation[1]
  }
  
  # Check and set mean
  if (is.null(prior.nu[["mean"]])) {
    prior.nu[["mean"]] <- min(1, nu.upper.bound / 2)
  } else if (length(prior.nu[["mean"]]) != 1) {
    warning("'prior.nu$mean' has length > 1. Only the first element will be used.")
    prior.nu[["mean"]] <- prior.nu[["mean"]][1]
  }
  
  # Check and set prec
  if (is.null(prior.nu$prec)) {
    mu_temp <- prior.nu[["mean"]] / nu.upper.bound
    prior.nu$prec <- max(1 / mu_temp, 1 / (1 - mu_temp)) + nu.prec.inc
  } else if (length(prior.nu$prec) != 1) {
    warning("'prior.nu$prec' has length > 1. Only the first element will be used.")
    prior.nu$prec <- prior.nu$prec[1]
  }
  
  # Check and set logscale
  if (is.null(prior.nu[["logscale"]])) {
    prior.nu[["logscale"]] <- 1
  } else if (length(prior.nu[["logscale"]]) != 1) {
    warning("'prior.nu$logscale' has length > 1. Only the first element will be used.")
    prior.nu[["logscale"]] <- prior.nu[["logscale"]][1]
  }
  
  # Determine starting value for nu
  if (prior.nu.dist == "beta") {
    start.nu <- prior.nu[["mean"]]
  } else if (prior.nu.dist == "lognormal") {
    start.nu <- exp(prior.nu[["loglocation"]])
  } else {
    stop("prior.nu.dist should be either 'beta' or 'lognormal'.")
  }
  
  # Validate start.nu range
  if (start.nu > nu.upper.bound || start.nu < 0) {
    if (prior.nu.dist == "beta") {
      stop("The 'mean' element of 'prior.nu' should be a number between 0 and nu.upper.bound!")
    } else {
      stop("The 'loglocation' element of 'prior.nu' should be a number less than log(nu.upper.bound)!")
    }
  }
  
  return(list(prior.nu = prior.nu, start.nu = start.nu))
}

#' Transform Anisotropic SPDE Model Parameters to Original Scale
#'
#' @description
#' This function takes a vector of transformed parameters and applies the appropriate
#' transformations to return them in the original scale for use in anisotropic SPDE models.
#'
#' @param theta A numeric vector of length 4 or 5, containing the transformed parameters in this order:
#' \describe{
#'   \item{lhx}{The logarithmic representation of hx.}
#'   \item{lhy}{The logarithmic representation of hy.}
#'   \item{logit_hxy}{The logit-transformed representation of hxy.}
#'   \item{lsigma}{The logarithmic representation of sigma.}
#'   \item{lnu (optional)}{The logarithmic representation of nu. If not provided, nu is not returned.}
#' }
#' @param nu_upper_bound (optional) A numeric value representing the upper bound for the smoothness parameter nu.
#' This is only used, and must be provided, if `lnu` is provided.
#'
#' @return A named list with the parameters in the original scale:
#' \describe{
#'   \item{hx}{The original scale for hx (exponential of lhx).}
#'   \item{hy}{The original scale for hy (exponential of lhy).}
#'   \item{hxy}{The original scale for hxy (inverse logit transformation of logit_hxy).}
#'   \item{sigma}{The original scale for sigma (exponential of lsigma).}
#'   \item{nu (optional)}{The original scale for nu (using the forward_nu transformation). Only included if `lnu` is provided.}
#' }
#' @export
#'
#' @examples
#' # With lnu
#' theta <- c(log(0.1), log(0.2), log((0.3 + 1) / (1 - 0.3)), log(0.5), log(1))
#' nu_upper_bound <- 2
#' transform_parameters_anisotropic(theta, nu_upper_bound)
#'
#' # Without lnu
#' theta <- c(log(0.1), log(0.2), log((0.3 + 1) / (1 - 0.3)), log(0.5))
#' transform_parameters_anisotropic(theta)
transform_parameters_anisotropic <- function(theta, nu_upper_bound = NULL) {
  if (!(length(theta) %in% c(4, 5))) {
    stop("Theta must be a numeric vector of length 4 or 5.")
  }
  
  # Functions for transformations
  adjusted_inv_logit <- function(z) {
    (2 / (1 + exp(-z))) - 1
  }
  
  forward_nu <- function(lnu, nu_upper_bound) {
    exp(lnu) / (1 + exp(lnu)) * nu_upper_bound
  }
  
  # Extract parameters
  lhx <- theta[1]
  lhy <- theta[2]
  logit_hxy <- theta[3]
  lsigma <- theta[4]
  
  # Transform parameters to original scale
  hx <- exp(lhx)
  hy <- exp(lhy)
  hxy <- adjusted_inv_logit(logit_hxy)
  sigma <- exp(lsigma)
  
  # Prepare the output
  result <- list(hx = hx, hy = hy, hxy = hxy, sigma = sigma)
  
  # If lnu is provided, compute nu
  if (length(theta) == 5) {
    if (is.null(nu_upper_bound)) {
      stop("nu_upper_bound must be provided if lnu is included in theta.")
    }
    lnu <- theta[5]
    nu <- forward_nu(lnu, nu_upper_bound)
    result$nu <- nu
  }
  
  return(result)
}


#' @noRd

find_inla_lib_path <- function() {
    # First check if INLA is installed
    if (!requireNamespace("INLA", quietly = TRUE)) {
        warning("INLA package is not installed")
        return(NULL)
    }
    
    # Get the base INLA bin directory
    inla_bin_path <- system.file("bin", package = "INLA")
    
    if (inla_bin_path == "") {
        warning("INLA bin directory not found")
        return(NULL)
    }
    
    # Determine OS and architecture
    os <- .Platform$OS.type
    arch <- R.Version()$arch
    
    # Initialize path
    lib_path <- NULL
    
    if (os == "windows") {
        # For Windows - always use windows/64bit or windows/32bit
        base_path <- file.path(inla_bin_path, "windows")
        if (dir.exists(base_path)) {
            lib_path <- if (grepl("64", arch)) {
                file.path(base_path, "64bit")
            } else {
                file.path(base_path, "32bit")
            }
        }
    } else if (os == "unix") {
        if (Sys.info()["sysname"] == "Darwin") {
            # For macOS - special case for ARM64
            if (grepl("arm64|aarch64", arch)) {
                lib_path <- file.path(inla_bin_path, "mac.arm64")
            } else {
                # For Intel Mac
                base_path <- file.path(inla_bin_path, "mac")
                if (dir.exists(base_path)) {
                    lib_path <- if (grepl("64", arch)) {
                        file.path(base_path, "64bit")
                    } else {
                        file.path(base_path, "32bit")
                    }
                }
            }
        } else {
            # For Linux - always use linux/64bit or linux/32bit
            base_path <- file.path(inla_bin_path, "linux")
            if (dir.exists(base_path)) {
                lib_path <- if (grepl("64", arch)) {
                    file.path(base_path, "64bit")
                } else {
                    file.path(base_path, "32bit")
                }
            }
        }
    }
    
    if (is.null(lib_path)) {
        warning("Could not determine appropriate library path")
        return(NULL)
    }
    
    if (!dir.exists(lib_path)) {
        warning(sprintf("Directory does not exist: %s", lib_path))
        return(NULL)
    }
    
    return(lib_path)
}

#' @noRd 
rspde_check_cgeneric_symbol <- function(model) {
    # Ensure the required fields exist in the model object
    if (!"f" %in% names(model) || !"cgeneric" %in% names(model$f) || 
        !"shlib" %in% names(model$f$cgeneric) || !"model" %in% names(model$f$cgeneric)) {
        stop("There was a problem with the model creation.")
    }
    
    # Extract the shared library path and the symbol name
    shlib <- model$f$cgeneric$shlib
    symbol <- model$f$cgeneric$model
    
    # Check if the shared library exists
    if (!file.exists(shlib)) {
        stop(paste("The shared library", shlib, "does not exist."))
    }
    
    # Get R_HOME library path
    r_lib_path <- file.path(R.home("lib"))
    
    # Get INLA library path
    inla_lib_path <- find_inla_lib_path()
    
    # Set up library path environment variable based on OS
    if (.Platform$OS.type == "windows") {
        current_path <- Sys.getenv("PATH")
        new_path <- if (current_path == "") {
            paste(r_lib_path, inla_lib_path, sep = ";")
        } else {
            paste(current_path, r_lib_path, inla_lib_path, sep = ";")
        }
        Sys.setenv(PATH = new_path)
    } else {
        current_lib_path <- Sys.getenv("LD_LIBRARY_PATH")
        new_lib_path <- if (current_lib_path == "") {
            paste(r_lib_path, inla_lib_path, sep = ":")
        } else {
            paste(current_lib_path, r_lib_path, inla_lib_path, sep = ":")
        }
        Sys.setenv(LD_LIBRARY_PATH = new_lib_path)
    }
    
    # Use the `dyn.load` and `is.loaded` functions to check for the symbol
    tryCatch({
        dyn.load(shlib) # Load the shared library
        if (is.loaded(symbol)) {
            dyn.unload(shlib) # Unload if the symbol is available
            return(invisible(TRUE)) # Return silently
        } else {
            warning(paste0("The symbol '", symbol, "' is not available in the shared library. Please install the latest testing version of INLA. 
      If the problem persists after installing the latest testing version of INLA, please open an issue at https://github.com/davidbolin/rSPDE/issues, 
      requesting that this model be added to INLA."))
        }
        dyn.unload(shlib) # Ensure the library is unloaded
    }, error = function(e) {
        warning(paste0("Error while loading the shared library or checking the symbol: ", e$message, 
                       ". Please install the latest testing version of INLA. If the problem persists after installing the 
                   latest testing version of INLA, please open an issue at https://github.com/davidbolin/rSPDE/issues, 
                   requesting that this model be added to INLA."))
    })
    
    # Restore original environment variables
    if (.Platform$OS.type == "windows") {
        Sys.setenv(PATH = current_path)
    } else {
        Sys.setenv(LD_LIBRARY_PATH = current_lib_path)
    }
}

#' @noRd
match_with_tolerance <- function(input, loc, tolerance = 1e-6) {
  # Initialize a vector to store matched indices
  matched_indices <- integer(length(input))
  
  for (i in seq_along(input)) {
    # Find the indices in loc that match the current input element within the tolerance
    match_idx <- which(abs(loc - input[i]) <= tolerance)
    
    if (length(match_idx) == 0) {
      # If no match is found, throw an error
      stop(sprintf("Error: The input location %.10f is not present in the original locations used to create the model object.", input[i]))
    } else if (length(match_idx) > 1) {
      # Handle the case where multiple matches are found
      warning(sprintf("Warning: Multiple matches found for input location %.10f. Using the first match.", input[i]))
      match_idx <- match_idx[1]
    }
    
    # Store the matched index
    matched_indices[i] <- match_idx
  }
  
  return(matched_indices)
}


#' @noRd 
merge_with_tolerance <- function(original_data, new_data, by, tolerance = 1e-5) {
  # Ensure column names match by adding missing columns
  all_columns <- union(names(original_data), names(new_data))
  original_data[setdiff(all_columns, names(original_data))] <- NA
  new_data[setdiff(all_columns, names(new_data))] <- NA
  
  # Extract reference columns
  original_loc <- original_data[[by]]
  new_loc <- new_data[[by]]
  
  # Initialize the merged dataset
  merged_data <- original_data
  
  # Match rows from new_data to original_data within the tolerance
  for (i in seq_along(new_loc)) {
    diffs <- abs(original_loc - new_loc[i])
    if (any(diffs <= tolerance)) {
      # Find the closest match in original_data
      matched_index <- which.min(diffs)
      merged_row <- merged_data[matched_index, ]
      new_row <- new_data[i, ]
      
      # Exclude the `by` column from the merge
      columns_to_merge <- setdiff(names(new_data), by)
      
      # Check for conflicts and replace missing values in merged_row with new_row
      for (col in columns_to_merge) {
        if (!is.na(new_row[[col]])) {
          if (!is.na(merged_row[[col]]) && merged_row[[col]] != new_row[[col]]) {
            warning(sprintf(
              "Conflicting values in column '%s' for location '%s': original='%s', new='%s'. Using new value.",
              col, new_loc[i], merged_row[[col]], new_row[[col]]
            ))
          }
          merged_row[[col]] <- new_row[[col]]
        }
      }
      
      # Replace the row in merged_data
      merged_data[matched_index, ] <- merged_row
    } else {
      # Add unmatched rows from new_data directly
      merged_data <- rbind(merged_data, new_data[i, ])
    }
  }
  
  # Remove duplicates based on the `by` column
  merged_data <- merged_data[!duplicated(merged_data[[by]]), ]
  
  return(merged_data)
}


#' Transform Spacetime SPDE Model Parameters to Original Scale
#'
#' @description
#' This function takes a vector of transformed parameters and applies the appropriate
#' transformations to return them in the original scale for use in spacetime SPDE models.
#'
#' @param theta A numeric vector containing the transformed parameters in this order:
#' \describe{
#'   \item{lkappa}{The logarithmic representation of kappa.}
#'   \item{lsigma}{The logarithmic representation of sigma.}
#'   \item{lgamma}{The logarithmic representation of gamma.}
#'   \item{logit_rho (optional)}{The logit-transformed representation of rho, if drift = 1.}
#'   \item{logit_rho2 (optional)}{The logit-transformed representation of rho2, if drift = 1 and d = 2.}
#' }
#' @param st_model A list containing the spacetime model parameters:
#' \describe{
#'   \item{d}{The dimension (e.g., 1 or 2).}
#'   \item{bound}{The bound for rho and rho2.}
#'   \item{is_bounded}{A logical value indicating if rho and rho2 are bounded.}
#'   \item{drift}{A logical value indicating if drift is included in the model.}
#' }
#'
#' @return A named list with the parameters in the original scale:
#' \describe{
#'   \item{kappa}{The original scale for kappa (exponential of lkappa).}
#'   \item{sigma}{The original scale for sigma (exponential of lsigma).}
#'   \item{gamma}{The original scale for gamma (exponential of lgamma).}
#'   \item{rho (optional)}{The original scale for rho.}
#'   \item{rho2 (optional)}{The original scale for rho2, if d = 2.}
#' }
#' @export
transform_parameters_spacetime <- function(theta, st_model) {
  if (!is.list(st_model) || !all(c("d", "bound_rho", "is_bounded", "drift") %in% names(st_model))) {
    stop("st_model must be a list containing 'd', 'bound_rho', 'is_bounded', and 'drift'.")
  }
  
  # Extract model parameters
  d <- st_model$d
  bound <- as.double(st_model$bound_rho)
  is_bounded <- st_model$is_bounded
  drift <- st_model$drift
  
  # Functions for transformations
  adjusted_inv_logit <- function(z, L) {
    if (L <= 0) stop("Bound L must be positive.")
    L * (2 / (1 + exp(-z)) - 1)
  }
  
  # Transform required parameters
  lkappa <- theta[1]
  lsigma <- theta[2]
  lgamma <- theta[3]
  kappa <- exp(lkappa)
  sigma <- exp(lsigma)
  gamma <- exp(lgamma)
  
  result <- list(kappa = kappa, sigma = sigma, gamma = gamma)
  
  # Include rho and rho2 if drift is included
  if (drift) {
    if (is_bounded) {
      logit_rho <- theta[4]
      rho <- adjusted_inv_logit(logit_rho, bound)
    } else {
      rho <- theta[4]
    }
    result$rho <- rho
    
    # Include rho2 if d = 2
    if (d == 2) {
      if (is_bounded) {
        logit_rho2 <- theta[5]
        rho2 <- adjusted_inv_logit(logit_rho2, bound)
      } else {
        rho2 <- theta[5]
      }
      result$rho2 <- rho2
    } else {
      result$rho2 <- 0.0
    }
  } else {
    result$rho <- 0.0
    result$rho2 <- 0.0
  }
  
  return(result)
}

#' Remove the "(fixed)" suffix from parameter names
#'
#' @param param_list A list containing parameter names that might have "(fixed)" suffix
#' @return The same list with "(fixed)" removed from names
#' @noRd

clean_fixed_param_names <- function(param_list) {
  # Function to clean a single name
  clean_name <- function(name) {
    gsub(" \\(fixed\\)$", "", name)
  }
  
  # Process the list recursively
  if (is.list(param_list)) {
    # For each element in the list
    for (i in seq_along(param_list)) {
      if (is.list(param_list[[i]])) {
        # If it's a nested list, process recursively
        param_list[[i]] <- clean_fixed_param_names(param_list[[i]])
      } else if (!is.null(names(param_list[[i]]))) {
        # If it's a named vector, clean its names
        names(param_list[[i]]) <- sapply(names(param_list[[i]]), clean_name)
      }
    }
    
    # Also clean the names of the list itself
    if (!is.null(names(param_list))) {
      names(param_list) <- sapply(names(param_list), clean_name)
    }
  } else if (!is.null(names(param_list))) {
    # If it's a named vector, clean its names
    names(param_list) <- sapply(names(param_list), clean_name)
  }
  
  return(param_list)
}

#' Extracts starting values from previous_fit
#'
#' @noRd 
extract_starting_values <- function(previous_fit, fix_coeff = FALSE, model_options = NULL) {
  # Validate previous_fit
  if (is.null(previous_fit) || !inherits(previous_fit, "rspde_lme")) {
    return(model_options)
  }
  
  # Determine prefix based on fix_coeff
  prefix <- if (fix_coeff) "fix_" else "start_"
  
  # Initialize model_options_tmp
  model_options_tmp <- list()
  
  # Check if it's a non-stationary model (previous_fit$stationary is FALSE and inherits from 'CBrSPDEobj' or 'rSPDEobj')
  is_nonstationary <- !isTRUE(previous_fit$stationary) && 
                     (inherits(previous_fit$latent_model, "CBrSPDEobj") || 
                      inherits(previous_fit$latent_model, "rSPDEobj"))
  
  if (is_nonstationary) {
    # Handle the non-stationary case with Theta parameters
    random_effects <- previous_fit$coeff$random_effects
    param_names <- names(random_effects)
    
    # Extract all Theta parameters and put them in a vector
    theta_indices <- grep("^Theta", param_names)
    if (length(theta_indices) > 0) {
      theta_values <- random_effects[theta_indices]
      model_options_tmp[[paste0(prefix, "theta")]] <- unname(theta_values)
    }
    
    # Add other parameters (not Theta)
    non_theta_indices <- setdiff(seq_along(random_effects), theta_indices)
    for (i in non_theta_indices) {
      param_name <- param_names[i]
      model_options_tmp[[paste0(prefix, tolower(param_name))]] <- random_effects[[i]]
    }
  } else {
    # Regular case - use original approach
    model_options_tmp <- setNames(
      as.list(previous_fit$coeff$random_effects),
      paste0(prefix, names(previous_fit$coeff$random_effects))
    )
  }
  
  # Add sigma_e with appropriate prefix
  model_options_tmp[[paste0(prefix, "sigma_e")]] <- previous_fit$coeff$measurement_error[[1]]
  
  # If fix_coeff is TRUE, check for additional parameters
  if (fix_coeff) {
    # Check and add alpha if non-NULL
    if (!is.null(previous_fit$alpha)) {
      model_options_tmp$fix_alpha <- previous_fit$alpha
    }
    
    # Check and add beta if non-NULL
    if (!is.null(previous_fit$beta)) {
      model_options_tmp$fix_beta <- previous_fit$beta
    }
    
    # Check and add nu if non-NULL
    if (!is.null(previous_fit$nu)) {
      model_options_tmp$fix_nu <- previous_fit$nu
    }
  }
  
  # If user provided model_options, combine them
  if (!is.null(model_options)) {
    # Overwrite extracted options with user-provided options
    for (name in names(model_options)) {
      model_options_tmp[[name]] <- model_options[[name]]
    }
  }
  
  return(clean_fixed_param_names(model_options_tmp))
}

#' @noRd 

process_model_options <- function(model, model_options) {
  # Only process for standard non-spacetime/non-special models
  if (!inherits(model, "spacetimeobj") && 
      !inherits(model, "CBrSPDEobj2d") && 
      !inherits(model, "intrinsicCBrSPDEobj")) {
    
    d <- model$d
    
    # Process nu parameter - priority order:
    # fix_nu > fix_alpha > start_nu > start_alpha > model$nu
    if (!is.null(model_options$fix_nu)) {
      nu <- model_options$fix_nu
      model_options$fix_alpha <- nu + d/2
    } else if (!is.null(model_options$fix_alpha)) {
      nu <- model_options$fix_alpha - d/2
      model_options$fix_nu <- nu
    } else if (!is.null(model_options$start_nu)) {
      nu <- model_options$start_nu
      model_options$start_alpha <- nu + d/2
    } else if (!is.null(model_options$start_alpha)) {
      nu <- model_options$start_alpha - d/2
      model_options$start_nu <- nu
    } else if (!is.null(model$nu)) {
      nu <- model$nu
    } else {
      stop("There was some error processing model_options. Cannot determine 'nu' parameter.")
    }
    
    # Process range/kappa parameters - priority order for kappa:
    # fix_range > start_range > fix_kappa > start_kappa > model$kappa
    if (!is.null(model_options$fix_range)) {
      model_options$fix_kappa <- sqrt(8 * nu) / model_options$fix_range
    } else if (!is.null(model_options$start_range)) {
      model_options$start_kappa <- sqrt(8 * nu) / model_options$start_range
    } else if (is.null(model_options$fix_kappa) && is.null(model_options$start_kappa) && is.null(model$kappa)) {
      stop("There was some error processing model_options. Cannot determine 'kappa' parameter.")
    }
    
    # Get final kappa for calculations
    if (!is.null(model_options$fix_kappa)) {
      kappa <- model_options$fix_kappa
    } else if (!is.null(model_options$start_kappa)) {
      kappa <- model_options$start_kappa
    } else {
      kappa <- model$kappa
    }
    
    # Process sigma/tau parameters - priority order for tau:
    # fix_sigma > start_sigma > fix_tau > start_tau > model$tau
    if (!is.null(model_options[["fix_sigma"]])) {
      model_options[["fix_tau"]] <- sqrt(gamma(nu) / 
                                   (model_options[["fix_sigma"]]^2 * kappa^(2 * nu) * 
                                    (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    } else if (!is.null(model_options[["start_sigma"]])) {
      model_options[["start_tau"]] <- sqrt(gamma(nu) / 
                                     (model_options[["start_sigma"]]^2 * kappa^(2 * nu) * 
                                      (4 * pi)^(d / 2) * gamma(nu + d / 2)))
    } else if (is.null(model_options[["fix_tau"]]) && is.null(model_options[["start_tau"]]) && is.null(model$tau)) {
      stop("There was some error processing model_options. Cannot determine 'tau' parameter.")
    }
  }
  
  return(model_options)
}

#' @noRd 

general_checks_lme <- function(model, model_options) {
  # Check alpha and beta combination for intrinsic models
  if (inherits(model, "intrinsicCBrSPDEobj") && !is.null(model_options$fix_alpha) && !is.null(model_options$fix_beta)) {
    if (model_options$fix_alpha + model_options$fix_beta <= model$d/2) {
      stop("One must have alpha + beta > d/2.")
    }
  }
  
  # Check fix_alpha for non-intrinsic models and non-space-time models
  if (!inherits(model, "intrinsicCBrSPDEobj") && 
      !inherits(model, "spacetimeobj") && 
      !is.null(model_options$fix_alpha) && 
      model_options$fix_alpha <= model$d / 2) {
    stop(paste("model_options$fix_alpha must be greater than dim/2 =", model$d / 2))
  }
  
  # Check start_alpha for non-intrinsic models
  if (!inherits(model, "intrinsicCBrSPDEobj") && 
      !inherits(model, "spacetimeobj") &&   
      !is.null(model_options$start_alpha) && 
      model_options$start_alpha <= model$d / 2) {
    stop(paste("model_options$start_alpha must be greater than dim/2 =", model$d / 2))
  }
  
  invisible(NULL)
}

#' @noRd 

get_model_starting_values <- function(model) {
  # Check model inheritance types
  spacetime <- inherits(model, "spacetimeobj")
  anisotropic <- inherits(model, "CBrSPDEobj2d")
  intrinsic <- inherits(model, "intrinsicCBrSPDEobj")
  
  # Non-stationary, non-spacetime, non-intrinsic case
  if (!model$stationary && !spacetime && !intrinsic) {
    if (is.null(model$theta)) {
      return(NULL)  # Will be handled by model_options in the second function
    }
    
    starting_values <- model$theta
    
  # Stationary, non-spacetime case
  } else if (model$stationary && !spacetime) {
    if (anisotropic) {
      # Initialize anisotropic parameters
      starting_values <- c(
        log(model$sigma), 
        log(model$hx), 
        log(model$hy), 
        -log(2/(model$hxy+1) - 1)
      )
      names(starting_values) <- c("sigma", "hx", "hy", "hxy")
      
    } else if (intrinsic) {
      if (model$alpha == 0) {
        starting_values <- log(c(model$tau))
        names(starting_values) <- c("tau")
      } else {
        starting_values <- log(c(model$tau, max(c(model$kappa, 1e-5))))
        names(starting_values) <- c("tau", "kappa")
      }
    } else {
      starting_values <- log(c(model$tau, max(c(model$kappa, 1e-5))))
      names(starting_values) <- c("tau", "kappa")
    }
  
  # Spacetime case
  } else {
    if (model$alpha == 0) {
      starting_values <- c(log(model$kappa), log(model$sigma), log(model$gamma))
      names(starting_values) <- c("kappa", "sigma", "gamma")
    } else {
      starting_values <- c(log(model$kappa), log(model$sigma), log(model$gamma), model$rho)
      names(starting_values) <- c("kappa", "sigma", "gamma", "rho")
    }
  }
  
  return(starting_values)
}

#' @noRd 

update_starting_values <- function(starting_values, model, model_options) {
  # Check model inheritance types
  spacetime <- inherits(model, "spacetimeobj")
  anisotropic <- inherits(model, "CBrSPDEobj2d")
  intrinsic <- inherits(model, "intrinsicCBrSPDEobj")
  
  # For non-stationary models, handle theta specially
  if (!model$stationary && !spacetime && !intrinsic) {
    if (is.null(starting_values)) {
      # No model$theta, must use model_options
      if (is.null(model_options$start_theta) && is.null(model_options$fix_theta)) {
        stop("For models given by spde.matern.operators(), either model$theta must be non-null or model_options$start_theta or model_options$fix_theta must be non-NULL!")
      }
      
      if (!is.null(model_options$fix_theta)) {
        # Check length if parameterization provided
        if (!is.null(parameterization) && parameterization == "spde") {
          if (length(model_options$fix_theta) != ncol(model$B.tau) - 1) {
            stop("Length of fix_theta must be equal to ncol(model$B.tau) - 1 when parameterization is 'spde'")
          }
        } else if (!is.null(parameterization)) {
          if (length(model_options$fix_theta) != ncol(model$B.sigma) - 1) {
            stop("Length of fix_theta must be equal to ncol(model$B.sigma) - 1 when parameterization is not 'spde'")
          }
        }
        starting_values <- model_options$fix_theta
      } else {
        # Use start_theta
        if (!is.null(parameterization) && parameterization == "spde") {
          if (length(model_options$start_theta) != ncol(model$B.tau) - 1) {
            stop("Length of start_theta must be equal to ncol(model$B.tau) - 1 when parameterization is 'spde'")
          }
        } else if (!is.null(parameterization)) {
          if (length(model_options$start_theta) != ncol(model$B.sigma) - 1) {
            stop("Length of start_theta must be equal to ncol(model$B.sigma) - 1 when parameterization is not 'spde'")
          }
        }
        starting_values <- model_options$start_theta
      }
      return(starting_values)
    } else if (!is.null(model_options$fix_theta)) {
      return(model_options$fix_theta)
    } else if (!is.null(model_options$start_theta)) {
      return(model_options$start_theta)
    } else {
      return(starting_values)
    }
  }
  
  # For all other models, loop through parameter names and update values
  if (!is.null(names(starting_values))) {
    for (param_name in names(starting_values)) {
      fix_param <- paste0("fix_", param_name)
      start_param <- paste0("start_", param_name)
      
      if (!is.null(model_options[[fix_param]])) {
        # Special handling for hxy which has a transformation
        if (param_name == "hxy") {
          starting_values[param_name] <- -log(2/(model_options[[fix_param]]+1) - 1)
        } else {
          starting_values[param_name] <- log(model_options[[fix_param]])
        }
      } else if (!is.null(model_options[[start_param]])) {
        if (param_name == "hxy") {
          starting_values[param_name] <- -log(2/(model_options[[start_param]]+1) - 1)
        } else {
          starting_values[param_name] <- log(model_options[[start_param]])
        }
      }
      
      # Special case: rho doesn't need log transformation
      if (param_name == "rho" && !is.null(model_options[[fix_param]])) {
        starting_values[param_name] <- model_options[[fix_param]]
      } else if (param_name == "rho" && !is.null(model_options[[start_param]])) {
        starting_values[param_name] <- model_options[[start_param]]
      }
    }
  }
  
  return(starting_values)
}

#' Get starting values for LME model
#'
#' @param model The model object
#' @param model_options Options including fixed or starting values
#' @param y_resp Response variable
#' @param starting_values_aux Auxiliary starting values for latent field parameters
#' @return Named vector of starting values for optimization
#' @noRd
#' 
get_starting_values_lme <- function(model, model_options, y_resp, starting_values_aux) {
  # Check model inheritance types
  intrinsic <- inherits(model, "intrinsicCBrSPDEobj")
  spacetime <- inherits(model, "spacetimeobj")
  
  # Set up estimate_nu and nu based on fix_nu or fix_alpha
  estimate_nu <- TRUE
  nu <- NULL
  if (!is.null(model_options$fix_nu)) {
    estimate_nu <- FALSE
    nu <- model_options$fix_nu
  } else if (!is.null(model_options$fix_alpha)) {
    estimate_nu <- FALSE
    nu <- model_options$fix_alpha - model$d/2
  }

  if(spacetime){
    estimate_nu <- FALSE
  } else if(!is.null(model_options$fix_nu) || !is.null(model_options$fix_alpha)){
    estimate_nu <- FALSE
  }
    
  # Handle starting values for alpha and beta - prioritize fix_ values over start_ values
  start_alpha <- if (!is.null(model_options$fix_alpha)) {
    model_options$fix_alpha
  } else {
    model_options$start_alpha
  }
  
  start_beta <- if (!is.null(model_options$fix_beta)) {
    model_options$fix_beta
  } else {
    model_options$start_beta
  }
  
  start_nu <- if (!is.null(model_options$fix_nu)) {
    model_options$fix_nu
  } else {
    model_options$start_nu
  }
  
  start_sigma_e <- if (!is.null(model_options$fix_sigma_e)) {
    model_options$fix_sigma_e
  } else {
    model_options$start_sigma_e
  }
  
  # Determine starting values based on model type and parameters
  if (intrinsic) {
    if (is.null(start_alpha)) {
      start_alpha = 1
    }
    if (is.null(start_beta)) {
      start_beta = 0.9
    }
    
    estimate_alpha <- is.null(model_options$fix_alpha)
    estimate_beta <- is.null(model_options$fix_beta)
    
    if (estimate_alpha && estimate_beta) {
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp)), 
        nu = log(start_alpha - model$d/2), # alpha > d/2, alpha = d/2 + exp(theta)
        beta = beta2theta(start_beta, model$d)
      )
      start_values <- c(start_values, starting_values_aux)
    } else if (estimate_alpha) {
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp)), 
        nu = log(start_alpha - model$d/2)
      )
      start_values <- c(start_values, starting_values_aux)
    } else if (estimate_beta) {
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp)), 
        beta = beta2theta(start_beta, model$d)
      )
      start_values <- c(start_values, starting_values_aux)
    } else {
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp))
      )
      start_values <- c(start_values, starting_values_aux)
    }
  } else if (estimate_nu) {
    if (is.null(start_nu) && is.null(start_alpha)) {
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp)),
        nu = log(model$nu)
      )
      start_values <- c(start_values, starting_values_aux)
    } else if (!is.null(start_nu)) {
      if (!is.numeric(start_nu)) {
        stop("start_nu must be numeric.")
      }
      if (length(start_nu) > 1) {
        stop("start_nu must have length 1.")
      }
      if (start_nu <= 0) {
        stop("start_nu must be positive")
      }
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp)),
        nu = log(start_nu)
      )
      start_values <- c(start_values, starting_values_aux)
    } else {
      if (!is.numeric(start_alpha)) {
        stop("start_alpha must be numeric.")
      }
      if (length(start_alpha) > 1) {
        stop("start_alpha must have length 1.")
      }
      if (start_alpha <= model$d / 2) {
        stop(paste("start_alpha must be greater than dim/2 =", model$d / 2))
      }
      start_values <- c(
        sigma_e = log(0.1 * sd(y_resp)),
        nu = log(start_alpha - model$d / 2)
      )
      start_values <- c(start_values, starting_values_aux)
    }
  } else {
    start_values <- c(
      sigma_e = log(0.1 * sd(y_resp))
    )
    start_values <- c(start_values, starting_values_aux)
  } 
  
  # Override sigma_e if specified
  if (!is.null(start_sigma_e)) {
    start_values["sigma_e"] <- log(start_sigma_e)
  }
  
  return(start_values)
}


#' Get appropriate auxiliary likelihood function based on model type
#'
#' @param model Model object
#' @return An auxiliary likelihood function appropriate for the model class
#' @noRd
get_aux_likelihood_function <- function(model) {
  if(inherits(model, "intrinsicCBrSPDEobj")) {
    return(aux_lme_intrinsic.loglike)
  } else if(inherits(model, "CBrSPDEobj2d")) {
    return(aux_lme_CBrSPDE.matern2d.loglike)
  } else if(inherits(model, "CBrSPDEobj")) {
    return(aux_lme_CBrSPDE.matern.loglike)
  } else if(inherits(model, "spacetimeobj")) {
    return(aux_lme_spacetime.loglike)
  } else if(inherits(model, "rSPDEobj1d")) {
    return(aux_lme_rSPDE.matern.rational.loglike)
  } else if(inherits(model, "rSPDEobj")) {
    return(aux_lme_rSPDE.matern.loglike)
  } else {
    stop("Unsupported model class")
  }
}

#' Determine which parameters should be estimated
#'
#' @param model Model object
#' @param model_options List of options including fixed parameters
#' @param start_values Named vector of starting values
#' @return Named logical vector indicating which parameters to estimate
#' @noRd
determine_estimate_params <- function(model, model_options, start_values) {
  # Initialize estimate_params with TRUE for all named parameters in start_values
  estimate_params <- setNames(
    rep(TRUE, length(start_values)), 
    names(start_values)
  )
  
  # For each parameter in start_values, check if it should be fixed
  for (param_name in names(start_values)) {
    # Skip unnamed parameters (handled separately below)
    if (param_name == "") {
      next
    }
    
    # Check if there's a corresponding fix_param in model_options
    fix_param_name <- paste0("fix_", param_name)
    if (!is.null(model_options[[fix_param_name]])) {
      estimate_params[param_name] <- FALSE
    }
  }
  
  # Special case for spacetime models - nu is never estimated
  if (inherits(model, "spacetimeobj") && "nu" %in% names(estimate_params)) {
    estimate_params["nu"] <- FALSE
  }

  # Special case for spacetime models - alpha is never estimated
  if (inherits(model, "spacetimeobj") && "alpha" %in% names(estimate_params)) {
    estimate_params["alpha"] <- FALSE
  }  
  
  # Handle overspecification
  if ("alpha" %in% names(estimate_params) && "nu" %in% names(estimate_params)) {
    stop("There was a processing error, alpha and nu were found at the same time at start_values.")
  }
  
  # Handle non-stationary models with unnamed parameters (fix_theta case)
  unnamed_params <- which(names(start_values) == "")
  if (length(unnamed_params) > 0 && !is.null(model_options$fix_theta)) {
    # Only apply for non-stationary models that could have theta parameters
    if (!model$stationary && !inherits(model, "spacetimeobj") && 
        !inherits(model, "CBrSPDEobj2d") && !inherits(model, "intrinsicCBrSPDEobj")) {
      
      # Find theta parameter if it exists in named parameters
      if ("theta" %in% names(estimate_params)) {
        estimate_params["theta"] <- FALSE
      } else {
        # Otherwise set all unnamed parameters to be fixed
        for (i in unnamed_params) {
          if (i <= length(estimate_params)) {
            # Create a new entry (assuming unnamed params are indexed)
            estimate_params[i] <- FALSE
          }
        }
      }
    }
  }
  
  return(estimate_params)
}


#' Extract model update parameters from theta vector
#'
#' @param model Model object
#' @param theta Parameter vector for optimization
#' @param estimate_params Logical vector indicating which parameters to estimate
#' @param model_options Model options with fixed parameter values
#' @param start_values Original starting values vector with names
#' @param n_coeff_nonfixed Number of non-fixed coefficients
#' @param nu_upper_bound Upper bound for nu parameter
#' @return List containing update arguments and additional values
#' @noRd
extract_model_update_args <- function(model, theta, estimate_params, model_options, 
                                      start_values, n_coeff_nonfixed, nu_upper_bound) {
  # Initialize results
  args_list <- list()
  result <- list(
    args_list = args_list,
    sigma_e = NULL,
    beta_cov = NULL,
    early_return = NULL
  )
  
  # Get parameter names from estimate_params
  param_names <- names(estimate_params)
  
  # Initialize index tracker
  index <- 1
  
  # Process parameters based on their names rather than model class
  for (param_name in param_names) {
    # Skip unnamed parameters in start_values (handled specially)
    if (param_name == "") {
      next
    }
    
    # If parameter should be estimated, get from theta
    if (estimate_params[param_name]) {
      # Special handling for each parameter type
      if (param_name == "sigma_e") {
        result$sigma_e <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "nu") {
        nu <- exp(theta[index])
        if (nu %% 1 == 0) nu <- nu - 1e-5
        nu <- min(nu, nu_upper_bound)
        args_list$nu <- nu
        
        # Special check for rSPDEobj1d at upper bound
        if (inherits(model, "rSPDEobj1d") && nu == nu_upper_bound) {
          result$early_return <- -10^100
          return(result)
        }
        
        # Also compute alpha if this is a model that uses alpha derived from nu
        if (!inherits(model, "CBrSPDEobj2d") && !inherits(model, "spacetimeobj")) {
          alpha <- nu + model$d/2
          alpha <- max(1e-5 + model$d/2, alpha)
          args_list$alpha <- alpha
        }
        
        index <- index + 1
      }
      else if (param_name == "alpha") {
        alpha <- exp(theta[index]) + model$d/2
        if (alpha %% 1 == 0) alpha <- alpha - 1e-5
        args_list$alpha <- alpha
        index <- index + 1
      }
      else if (param_name == "beta") {
        beta <- theta2beta(theta[index], model$d)
        if (beta %% 1 == 0) beta <- beta - 1e-5
        args_list$beta <- beta
        index <- index + 1
      }
      else if (param_name == "tau") {
        args_list$tau <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "kappa") {
        args_list$kappa <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "sigma") {
        args_list$sigma <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "gamma") {
        args_list$gamma <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "hx") {
        args_list$hx <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "hy") {
        args_list$hy <- exp(theta[index])
        index <- index + 1
      }
      else if (param_name == "hxy") {
        args_list$hxy <- 2*exp(theta[index])/(1+exp(theta[index])) - 1
        index <- index + 1
      }
      else if (param_name == "rho") {
        # Handle rho parameter (vector) for spacetime models
        if (inherits(model, "spacetimeobj") && model$alpha > 0) {
          if (model$is_bounded_rho) {
            args_list$rho <- model$bound_rho * (2.0 / (1.0 + exp(-theta[index:(index+model$d-1)])) - 1.0)
          } else {
            args_list$rho <- theta[index:(index+model$d-1)]
          }
          index <- index + model$d
        } else {
          args_list$rho <- rep(0, model$d)
        }
      }
      else if (param_name == "theta") {
        # Handle theta parameter for non-stationary models
        if (!model$stationary) {
          theta_length <- length(model$theta)
          args_list$theta <- theta[index:(index + theta_length - 1)]
          index <- index + theta_length
        }
      }
    } else {
      # If parameter is fixed, get from model_options or model
      fix_param_name <- paste0("fix_", param_name)
      
      if (param_name == "sigma_e") {
        result$sigma_e <- model_options[[fix_param_name]]
      }
      else if (param_name == "nu") {
        if (!is.null(model_options[[fix_param_name]])) {
          args_list$nu <- model_options[[fix_param_name]]
        } else if (!is.null(model_options$fix_alpha)) {
          args_list$nu <- model_options$fix_alpha - model$d/2
        } else {
          args_list$nu <- model$nu
        }
        
        # Also compute alpha if this is a model that uses alpha derived from nu
        if (!inherits(model, "CBrSPDEobj2d") && !inherits(model, "spacetimeobj")) {
          args_list$alpha <- args_list$nu + model$d/2
        }
      }
      else if (param_name == "alpha" && !is.null(model_options[[fix_param_name]])) {
        args_list$alpha <- model_options[[fix_param_name]]
      }
      else if (param_name != "sigma_e" && param_name != "nu" && param_name != "alpha") {
        # For other parameters, just get from model_options or model
        if (!is.null(model_options[[fix_param_name]])) {
          args_list[[param_name]] <- model_options[[fix_param_name]]
        } else if (!is.null(model[[param_name]])) {
          args_list[[param_name]] <- model[[param_name]]
        }
      }
    }
  }
  
  # Special handling for non-stationary models with fix_theta
  if (!model$stationary && !inherits(model, "spacetimeobj") && 
      !inherits(model, "CBrSPDEobj2d") && !inherits(model, "intrinsicCBrSPDEobj")) {
    
    # Check if we have unnamed parameters in start_values
    unnamed_indices <- which(names(start_values) == "")
    if (length(unnamed_indices) > 0 && !is.null(model_options$fix_theta)) {
      args_list$theta <- model_options$fix_theta
    }
  }
  
  # Handle special case for intrinsicCBrSPDEobj with alpha=0
  if (inherits(model, "intrinsicCBrSPDEobj")) {
    if (!is.null(args_list$alpha) && args_list$alpha == 0) {
      args_list$kappa <- 0
    }
  }
  
  # Special case for spacetime models with alpha=0
  if (inherits(model, "spacetimeobj")) {
    if (model$alpha == 0) {
      args_list$rho <- rep(0, model$d)
    }
  }
  
  # Extract beta_cov if needed
  n_cov <- length(theta) - n_coeff_nonfixed
  if (n_cov > 0) {
    result$beta_cov <- theta[(n_coeff_nonfixed + 1):length(theta)]
  }
  
  result$args_list <- args_list
  return(result)
}

#' Get arguments for auxiliary likelihood function based on model class
#'
#' @param model Updated model object
#' @param y_resp Response variable
#' @param X_cov Covariate matrix
#' @param repl Replication indicator
#' @param A_list List of observation matrices
#' @param sigma_e Sigma_e parameter
#' @param beta_cov Beta coefficients for covariates
#' @param mean_correction Apply mean correction 
#' @param loc_df Location data frame
#' @return List of arguments for auxiliary likelihood function
#' @noRd
get_aux_lik_fun_args <- function(model, y_resp, X_cov, repl, A_list, 
                                sigma_e, beta_cov, mean_correction = FALSE, 
                                loc_df = NULL) {
  # Common arguments for all auxiliary likelihood functions
  args <- list(
    object = model,
    y = y_resp,
    X_cov = X_cov,
    repl = repl,
    sigma_e = sigma_e,
    beta_cov = beta_cov
  )
  
  # Model-specific arguments
  if(inherits(model, "rSPDEobj1d")) {
    args$loc <- loc_df
  } else {
    args$A_list <- A_list
  }
  
  # Only intrinsicCBrSPDEobj models use mean_correction
  if(inherits(model, "intrinsicCBrSPDEobj")) {
    args$mean_correction <- mean_correction
  }
  
  return(args)
}

#' Create likelihood function for model fitting
#'
#' @param model Original model object
#' @param model_options Model options including fixed/start parameters
#' @param y_resp Response variable
#' @param X_cov Covariate matrix 
#' @param A_list List of observation matrices 
#' @param repl Replication indicator 
#' @param start_values Named vector of starting values
#' @param mean_correction Apply mean correction 
#' @param nu_upper_bound Upper bound for nu parameter 
#' @param loc_df Location data frame
#' @return List containing likelihood function and parameter estimation flags
#' @noRd
create_likelihood <- function(model, model_options, y_resp, 
                             X_cov, A_list,
                             repl, start_values,
                             mean_correction, nu_upper_bound,
                             loc_df) {
  
  # Initialize X_cov if NULL
  if(is.null(X_cov)) {
    X_cov <- matrix(0, nrow = length(y_resp), ncol = 0)
  }
  
  # Get appropriate auxiliary likelihood function
  aux_lik_fun <- get_aux_likelihood_function(model)
  
  # Determine which parameters to estimate
  estimate_params <- determine_estimate_params(model, model_options, start_values)
  
  # Count number of non-fixed coefficients from model parameters
  n_coeff_nonfixed <- sum(estimate_params)
  
  # Create the likelihood function
  likelihood <- function(theta) {
    # Create a working copy of the model
    model_tmp <- model
    
    # Extract model update arguments
    result <- extract_model_update_args(
      model = model_tmp, 
      theta = theta, 
      estimate_params = estimate_params,
      model_options = model_options,
      start_values = start_values,
      n_coeff_nonfixed = n_coeff_nonfixed,
      nu_upper_bound = nu_upper_bound
    )
    
    # Check for early return (e.g., nu at upper bound for rSPDEobj1d)
    if(!is.null(result$early_return)) {
      return(result$early_return)
    }
    
    # Update the model with the extracted parameters
    model_tmp <- do.call(update, c(list(object = model_tmp), result$args_list))
    
    # Get arguments for auxiliary likelihood function
    aux_args <- get_aux_lik_fun_args(
      model = model_tmp,
      y_resp = y_resp,
      X_cov = X_cov,
      repl = repl,
      A_list = A_list,
      sigma_e = result$sigma_e,
      beta_cov = result$beta_cov,
      mean_correction = mean_correction,
      loc_df = loc_df
    )
    
    # Call the auxiliary likelihood function with the appropriate arguments
    loglik <- do.call(aux_lik_fun, aux_args)
    
    return(-loglik)
  }
  
  # Return both the likelihood function and parameter estimation flags
  return(list(
    likelihood = likelihood,
    estimate_params = estimate_params,
    n_coeff_nonfixed = n_coeff_nonfixed
  ))
}


#' Extract parameters from optimization results
#'
#' @param res Optimization result object
#' @param start_values Named vector of starting values
#' @param estimate_params Logical vector indicating which parameters to estimate
#' @param model Original model object
#' @param model_options Model options with fixed parameter values
#' @param X_cov Covariate matrix
#' @param n_coeff_nonfixed Number of non-fixed coefficients
#' @return List containing parameter values and additional information
#' @noRd
extract_parameters_from_optim <- function(res, start_values, estimate_params, model, 
                                         model_options, X_cov, n_coeff_nonfixed) {
  # Initialize result list
  result <- list(
    coeff_random = NULL,
    estimated_alpha = NULL,
    estimated_beta = NULL
  )
  
  # Initialize tracking index for res$par
  index <- 1
  
  # Create coefficient vector with same length as start_values
  param_names <- names(start_values)
  coeff <- numeric(length(start_values))
  names(coeff) <- param_names
  
  # Process each parameter
  for (i in 1:length(start_values)) {
    param_name <- param_names[i]
    
    # If parameter is estimated, get from res$par
    if (estimate_params[i]) {
      if (param_name == "sigma_e" || param_name == "tau" || 
          param_name == "kappa" || param_name == "sigma" || 
          param_name == "gamma" || param_name == "hx" || param_name == "hy") {
        # Parameters with exponential transformation
        coeff[i] <- exp(res$par[index])
        index <- index + 1
      }
      else if (param_name == "nu") {
        # Nu parameter
        nu_value <- exp(res$par[index])
        if (nu_value %% 1 == 0) nu_value <- nu_value - 1e-5
        coeff[i] <- min(nu_value, 10)
        index <- index + 1
        
        # Also calculate estimated_alpha if we have nu
        if (!inherits(model, "spacetimeobj") && !inherits(model, "CBrSPDEobj2d")) {
          result$estimated_alpha <- coeff[i] + model$d/2
        }
      }
      else if (param_name == "alpha") {
        # Alpha parameter
        alpha_value <- exp(res$par[index]) + model$d/2
        if (alpha_value %% 1 == 0) alpha_value <- alpha_value - 1e-5
        coeff[i] <- alpha_value
        index <- index + 1
        
        # Store estimated alpha
        result$estimated_alpha <- alpha_value
      }
      else if (param_name == "beta") {
        # Beta parameter
        beta_value <- theta2beta(res$par[index], model$d)
        if (beta_value %% 1 == 0) beta_value <- beta_value - 1e-5
        coeff[i] <- beta_value
        index <- index + 1
        
        # Store estimated beta
        result$estimated_beta <- beta_value
      }
      else if (param_name == "hxy") {
        # hxy parameter (correlation)
        coeff[i] <- 2*exp(res$par[index])/(1+exp(res$par[index])) - 1
        index <- index + 1
      }
      else if (param_name == "rho") {
        # Handle rho parameter (vector) for spacetime models
        if (inherits(model, "spacetimeobj") && model$alpha > 0) {
          if (model$is_bounded_rho) {
            bound_rho <- model$bound_rho
            coeff[i] <- bound_rho * (2.0 / (1.0 + exp(-res$par[index:(index+model$d-1)])) - 1.0)
          } else {
            coeff[i] <- res$par[index:(index+model$d-1)]
          }
          index <- index + model$d
        } else {
          coeff[i] <- rep(0, model$d)
        }
      }
      else if (param_name == "theta" || param_name == "") {
        # Handle theta parameter for non-stationary models
        if (!model$stationary) {
          theta_length <- length(model$theta)
          coeff[i] <- res$par[index:(index + theta_length - 1)]
          index <- index + theta_length
        }
      }
      else {
        # Generic parameter with no transformation
        coeff[i] <- res$par[index]
        index <- index + 1
      }
    }
    else {
      # Parameter is fixed - get from model_options or model
      fix_param_name <- paste0("fix_", param_name)
      
      if (param_name == "sigma_e") {
        coeff[i] <- model_options[[fix_param_name]]
      }
      else if (param_name == "nu") {
        if (!is.null(model_options[[fix_param_name]])) {
          coeff[i] <- model_options[[fix_param_name]]
        } else if (!is.null(model_options$fix_alpha)) {
          coeff[i] <- model_options$fix_alpha - model$d/2
        } else {
          coeff[i] <- model$nu
        }
      }
      else if (param_name == "alpha") {
        if (!is.null(model_options[[fix_param_name]])) {
          coeff[i] <- model_options[[fix_param_name]]
        } else {
          coeff[i] <- model$alpha
        }
      }
      else if (param_name != "sigma_e" && param_name != "nu" && param_name != "alpha") {
        # For other parameters, just get from model_options or model
        if (!is.null(model_options[[fix_param_name]])) {
          coeff[i] <- model_options[[fix_param_name]]
        } else if (!is.null(model[[param_name]])) {
          coeff[i] <- model[[param_name]]
        }
      }
    }
  }
  
  # Add coefficients for covariates if any
  n_fixed <- ncol(X_cov)
  if (n_fixed > 0) {
    coeff_fixed <- res$par[(n_coeff_nonfixed + 1):(n_coeff_nonfixed + n_fixed)]
    result$coeff_fixed <- coeff_fixed
  }
  
  # Set result coefficient vector
  result$coeff_random <- coeff
  
  return(result)
}

#' Organize extracted parameters into appropriate categories
#'
#' @param coeff Full coefficient vector
#' @param model Model object
#' @param estimate_params Logical vector indicating which parameters were estimated
#' @param X_cov Covariate matrix
#' @return List containing organized parameter vectors
#' @noRd
organize_parameters <- function(coeff, model, estimate_params, X_cov) {
  result <- list(
    coeff_meas = NULL,
    coeff_random = NULL,
    coeff_fixed = NULL,
    par_names = NULL
  )
  
  # Set measurement error coefficient (always the first parameter)
  result$coeff_meas <- coeff[1]
  names(result$coeff_meas) <- "std. dev"
  
  # Determine parameter names based on model type
  if (inherits(model, "intrinsicCBrSPDEobj")) {
    par_names <- c()

    # Check if "alpha" and "beta" are present in names of estimate_params
    has_alpha_param <- "alpha" %in% names(estimate_params)
    has_beta_param <- "beta" %in% names(estimate_params)

    # Add parameters to par_names list based on whether they should be estimated
    if (has_alpha_param && estimate_params["alpha"]) {
      par_names <- c(par_names, "alpha")
    }

    if (has_beta_param && estimate_params["beta"]) {
      par_names <- c(par_names, "beta")
    }

    # Always add tau
    par_names <- c(par_names, "tau")

    # Add kappa only if alpha is not fixed to zero, or if alpha is not in estimate_params
    if ( (has_alpha_param && estimate_params["alpha"]) || (model["alpha"] > 0)) {
      par_names <- c(par_names, "kappa")
    }
  } else if (inherits(model, "CBrSPDEobj2d")) {
    if(!is.na(estimate_params["nu"])){
      estimate_nu <- estimate_params["nu"]
    } else if(!is.na(estimate_params["alpha"])){
      estimate_nu <- estimate_params["alpha"]
    } else{
      estimate_nu <- FALSE
    }

    if (estimate_nu) {
      par_names <- c("nu", "sigma", "hx", "hy", "hxy")
    } else {
      par_names <- c("sigma", "hx", "hy", "hxy")
    }
  } 
  else if ((inherits(model, "CBrSPDEobj") || inherits(model, "rSPDEobj"))) {
    if(!is.na(estimate_params["nu"])){
      estimate_nu <- estimate_params["nu"]
    } else if(!is.na(estimate_params["alpha"])){
      estimate_nu <- estimate_params["alpha"]
    } else{
      estimate_nu <- FALSE
    }
        
    if (model$stationary) {
      if (estimate_nu) {
        par_names <- c("nu", "tau", "kappa")
      } else {
        par_names <- c("tau", "kappa")
      }
    } else {
      # For non-stationary models, use B.tau to determine theta parameter names
      if (estimate_nu) {
        par_names <- c("nu")
      } else {
        par_names <- c()
      }
      
      # Add theta parameter names
      if (estimate_params["theta"] || !is.null(names(coeff)["theta"])) {
        if (ncol(model$B.tau) <= 2) {
          par_names <- c(par_names, "Theta 1")
        } else {
          for (i in 1:(ncol(model$B.tau) - 1)) {
            par_names <- c(par_names, paste("Theta", i))
          }
        }
      }
    }
  } 
  else if (inherits(model, "rSPDEobj1d")) {
    if(!is.na(estimate_params["nu"])){
      estimate_nu <- estimate_params["nu"]
    } else if(!is.na(estimate_params["alpha"])){
      estimate_nu <- estimate_params["alpha"]
    } else{
      estimate_nu <- FALSE
    }
        
    if (estimate_nu) {
      par_names <- c("nu", "tau", "kappa")
    } else {
      par_names <- c("tau", "kappa")
    }
  } 
  else if (inherits(model, "spacetimeobj")) {
    if (model$alpha == 0) {
      par_names <- c("kappa", "sigma", "gamma")
    } else {
      par_names <- c("kappa", "sigma", "gamma", "rho")
    }
  } 
  else {
    # Generic model
    par_names <- names(coeff)[-1]  # All except sigma_e
  }
  
  # Extract random effects coefficients
  # First find which parameters in coeff correspond to par_names
  param_idx <- numeric(0)
  param_names <- names(coeff)

  if(!model$stationary){
    if(!is.na(estimate_params["nu"])){
      estimate_nu <- estimate_params["nu"]
    } else if(!is.na(estimate_params["alpha"])){
      estimate_nu <- estimate_params["alpha"]
    } else{
      estimate_nu <- FALSE
    }     
  }
  
  for (name in par_names) {
    # For parameters like "Theta 1", we need to check for both the named parameter
    # and numeric position
    if (name %in% param_names) {
      param_idx <- c(param_idx, which(param_names == name))
    } else if (name == "theta" && "theta" %in% param_names) {
      param_idx <- c(param_idx, which(param_names == "theta"))
    } else if (startsWith(name, "Theta ") && !model$stationary) {
      # For non-stationary models with unnamed theta parameters
      # This could be in various positions depending on what else is estimated
      if ("theta" %in% param_names) {
        param_idx <- c(param_idx, which(param_names == "theta"))
      } else {
        # If we have unnamed parameters, attempt to find them by position
        theta_num <- as.numeric(substring(name, 7))
        # Start from position 2 (after sigma_e) and count
        if (estimate_nu) {
          pos <- 3 + theta_num - 1  # sigma_e, nu, then theta positions
        } else {
          pos <- 2 + theta_num - 1  # sigma_e, then theta positions
        }
        if (pos <= length(coeff)) {
          param_idx <- c(param_idx, pos)
        }
      }
    }
  }
  
  # If we couldn't match by name, use positions (theta is usually after sigma_e and nu)
  if (length(param_idx) == 0 && !model$stationary) {
    if (estimate_nu) {
      start_pos <- 3  # sigma_e, nu
    } else {
      start_pos <- 2  # sigma_e
    }
    
    param_idx <- start_pos:(start_pos + length(par_names) - 1)
    if (max(param_idx) > length(coeff)) {
      param_idx <- param_idx[param_idx <= length(coeff)]
    }
  }
  
  # If we still don't have indices, just use positions 2 through n+1
  if (length(param_idx) == 0) {
    param_idx <- 2:(1 + length(par_names))
    if (max(param_idx) > length(coeff)) {
      param_idx <- param_idx[param_idx <= length(coeff)]
    }
  }
  
  # Get the random effect coefficients
  result$coeff_random <- coeff[param_idx]

  print(par_names)

  names(result$coeff_random) <- par_names
  
  # Set parameter names
  result$par_names <- par_names
  
  # Extract fixed effects coefficients if any
  n_fixed <- ncol(X_cov)
  if (n_fixed > 0) {
    if (length(coeff) >= length(estimate_params)) {
      # If coeff contains unnamed elements at the end (fixed effects)
      result$coeff_fixed <- coeff[(length(estimate_params) + 1):length(coeff)]
    }
  }
  
  return(result)
}

#' Process optimization results into parameter estimates
#'
#' @param res Optimization result object
#' @param start_values Named vector of starting values
#' @param estimate_params Logical vector indicating which parameters to estimate
#' @param model Original model object
#' @param model_options Model options with fixed parameter values
#' @param X_cov Covariate matrix
#' @param n_coeff_nonfixed Number of non-fixed coefficients
#' @return List containing organized parameter estimates
#' @noRd
extract_and_organize_parameters <- function(res, start_values, estimate_params, model,
                                           model_options, X_cov, n_coeff_nonfixed) {
  
  # Extract parameters from optimization results
  extracted_params <- extract_parameters_from_optim(
    res = res,
    start_values = start_values,
    estimate_params = estimate_params,
    model = model,
    model_options = model_options,
    X_cov = X_cov,
    n_coeff_nonfixed = n_coeff_nonfixed
  )
  
  # Organize parameters into categories
  result <- organize_parameters(
    coeff = c(extracted_params$coeff_random, extracted_params$coeff_fixed),
    model = model,
    estimate_params = estimate_params,
    X_cov = X_cov
  )
  
  # Add estimated alpha/beta if available
  if (!is.null(extracted_params$estimated_alpha)) {
    result$estimated_alpha <- extracted_params$estimated_alpha
  }
  if (!is.null(extracted_params$estimated_beta)) {
    result$estimated_beta <- extracted_params$estimated_beta
  }
    
  return(result)
}

#' Calculate Jacobian of parameter transformation for Fisher information matrix
#'
#' @param res Optimization result object
#' @param estimate_params Logical vector indicating which parameters to estimate
#' @param model Original model object
#' @param model_options Model options with fixed parameter values
#' @param X_cov Covariate matrix
#' @return Matrix representing the parameter transformation Jacobian
#' @noRd
calculate_parameter_jacobian <- function(res, estimate_params, model, model_options, X_cov) {
  # Get dimensions
  n_coeff_nonfixed <- sum(estimate_params)
  n_fixed <- ncol(X_cov)
  n_total <- n_coeff_nonfixed + n_fixed
  
  # Create diagonal matrix for the Jacobian
  par_change <- diag(n_total)
  
  # Track position in res$par
  index <- 1
  
  # Process parameters in order of estimate_params
  param_names <- names(estimate_params)
  for (i in 1:length(estimate_params)) {
    param_name <- param_names[i]
    
    if (estimate_params[i]) {
      # Parameter is estimated, get transformation from res$par
      if (param_name == "sigma_e" || param_name == "tau" || 
          param_name == "kappa" || param_name == "sigma" || 
          param_name == "gamma" || param_name == "hx" || param_name == "hy") {
        # Parameters with exp transformation
        par_change[index, index] <- exp(-res$par[index])
        index <- index + 1
      }
      else if (param_name == "nu" || param_name == "alpha") {
        # Nu/alpha parameters
        par_change[index, index] <- exp(-res$par[index])
        index <- index + 1
      }
      else if (param_name == "beta") {
        # Beta parameter (special transformation)
        par_change[index, index] <- dbetadtheta(res$par[index])
        index <- index + 1
      }
      else if (param_name == "hxy") {
        # hxy parameter (correlation)
        hxy <- 2*exp(res$par[index])/(1+exp(res$par[index]))-1
        hxy_trans <- 2/(2*(hxy+1)-(hxy+1)^2)
        par_change[index, index] <- hxy_trans
        index <- index + 1
      }
      else if (param_name == "rho" && inherits(model, "spacetimeobj")) {
        # Rho parameter for spacetime models
        if (model$alpha > 0) {
          if (model$is_bounded_rho) {
            bound_rho <- model$bound_rho
            if (model$d == 1) {
              rho_trans <- bound_rho * 2.0 * exp(res$par[index])/ ((1.0 + exp(res$par[index]))^2)
              par_change[index, index] <- rho_trans
              index <- index + 1
            } else {
              # For multiple dimensions
              for (j in 1:model$d) {
                rho_trans <- bound_rho * 2.0 * exp(res$par[index])/ ((1.0 + exp(res$par[index]))^2)
                par_change[index, index] <- rho_trans
                index <- index + 1
              }
            }
          } else {
            # Unbounded rho - no transformation
            index <- index + model$d
          }
        }
      }
      else if (param_name == "theta" || param_name == "") {
        # For non-stationary models
        if (!model$stationary) {
          theta_length <- length(model$theta)
          # No transformation for theta parameters
          index <- index + theta_length
        }
      }
      else {
        # Generic parameter with no transformation
        index <- index + 1
      }
    }
  }
    
  return(par_change)
}

#' Calculate standard errors from optimization results
#'
#' @param observed_fisher Observed Fisher information matrix
#' @param res Optimization result object
#' @param estimate_params Logical vector indicating which parameters to estimate
#' @param model Original model object
#' @param model_options Model options with fixed parameter values
#' @param X_cov Covariate matrix
#' @param n_coeff_nonfixed Number of non-fixed coefficients
#' @param param_results Results from parameter extraction
#' @return List containing standard errors and other results
#' @noRd
calculate_standard_errors <- function(observed_fisher, res, estimate_params, model, 
                                     model_options, X_cov, n_coeff_nonfixed, param_results) {
  # Handle edge cases
  all_fixed = all(!estimate_params)
  only_sigma_e = sum(estimate_params) == 1 && estimate_params[1]
  
  # If all parameters are fixed, return NAs for all standard errors
  if (all_fixed) {
    std_err = rep(NA, length(estimate_params) + ncol(X_cov))
    std_meas = NA
    std_random = rep(NA, length(param_results$coeff_random))
    names(std_random) = names(param_results$coeff_random)
    std_fixed = NULL
    if (!is.null(param_results$coeff_fixed)) {
      std_fixed = rep(NA, length(param_results$coeff_fixed))
    }
    
    return(list(
      std_err = std_err,
      std_meas = std_meas,
      std_random = std_random,
      std_fixed = std_fixed,
      inv_fisher = matrix(NA, length(std_err), length(std_err))
    ))
  }
  
  # If only sigma_e is estimated (all latent parameters fixed), 
  # we still need to calculate its standard error
  if (only_sigma_e) {
    # If we only have one parameter, observed_fisher is just a 1x1 matrix
    if (is.null(dim(observed_fisher)) || (nrow(observed_fisher) == 1 && ncol(observed_fisher) == 1)) {
      inv_fisher = 1/observed_fisher
      std_err = sqrt(inv_fisher)
      std_meas = std_err[1]
      
      std_random = rep(NA, length(param_results$coeff_random))
      names(std_random) = names(param_results$coeff_random)
      
      std_fixed = NULL
      if (ncol(X_cov) > 0) {
        std_fixed = rep(NA, ncol(X_cov))
      }
      
      return(list(
        std_err = c(std_meas, rep(NA, length(param_results$coeff_random) + ncol(X_cov))),
        std_meas = std_meas,
        std_random = std_random,
        std_fixed = std_fixed,
        inv_fisher = matrix(inv_fisher, 1, 1)
      ))
    }
  }
  
  # Regular case - calculate parameter Jacobian and standard errors
  
  # Calculate the parameter transformation Jacobian
  par_change <- calculate_parameter_jacobian(
    res = res,
    estimate_params = estimate_params,
    model = model,
    model_options = model_options,
    X_cov = X_cov
  )
  
  # Apply parameter transformation to observed Fisher information
  observed_fisher <- par_change %*% observed_fisher %*% par_change
  
  # Attempt to invert the Fisher information matrix
  inv_fisher <- tryCatch(
    solve(observed_fisher), 
    error = function(e) matrix(NA, nrow(observed_fisher), ncol(observed_fisher))
  )
  
  # Calculate standard errors from inverse Fisher information
  std_err <- sqrt(diag(inv_fisher))
  
  # Create standard error vectors with same structure as parameter vectors
  # Initialize with NA
  std_random <- rep(NA, length(param_results$coeff_random))
  names(std_random) <- names(param_results$coeff_random)
  
  std_fixed <- NULL
  if (!is.null(param_results$coeff_fixed)) {
    std_fixed <- rep(NA, length(param_results$coeff_fixed))
  }
  
  # Get standard error for measurement error (sigma_e) - always first parameter
  if (estimate_params[1]) {
    std_meas <- std_err[1]
  } else {
    std_meas <- NA
  }
  
  # If we have any estimated random effect parameters
  if (sum(estimate_params) > 1 || (!estimate_params[1] && sum(estimate_params) > 0)) {
    # Get parameter names from the parameter results
    par_names <- param_results$par_names
    
    # Fill in standard errors for estimated random effect parameters
    # Track position in std_err vector
    index <- 2  # Start after sigma_e (or position 1 if sigma_e is not estimated)
    if (!estimate_params[1]) {
      index <- 1
    }
    
    # Get indices of estimated parameters (excluding sigma_e if it's estimated)
    est_param_indices <- which(estimate_params)
    if (estimate_params[1]) {
      est_param_indices <- est_param_indices[-1]
    }
    
    # Map estimated parameters to their positions in std_random
    for (i in 1:length(est_param_indices)) {
      param_name <- names(estimate_params)[est_param_indices[i]]
      
      # Find the position of this parameter in coeff_random
      pos <- which(param_name == names(std_random))
      
      # For parameters with special naming (like "Theta 1" vs "theta")
      if (length(pos) == 0 && param_name == "theta") {
        # Find positions that start with "Theta"
        pos <- grep("^Theta", names(std_random))
      }
      
      if (length(pos) > 0) {
        if (length(pos) == 1) {
          # Regular parameter
          std_random[pos] <- std_err[index]
        } else {
          # Vector parameter (like theta or rho)
          if (param_name == "rho" && inherits(model, "spacetimeobj")) {
            for (j in 1:model$d) {
              if (j <= length(pos) && index <= length(std_err)) {
                std_random[pos[j]] <- std_err[index]
                index <- index + 1
              }
            }
            next  # Skip the normal index increment
          } else if (param_name == "theta" || param_name == "") {
            # For theta parameters in non-stationary models
            theta_length <- length(model$theta)
            for (j in 1:min(theta_length, length(pos))) {
              if (index <= length(std_err)) {
                std_random[pos[j]] <- std_err[index]
                index <- index + 1
              }
            }
            next  # Skip the normal index increment
          }
        }
      }
      
      index <- index + 1
    }
  }
  
  # Fill in standard errors for fixed effect parameters (covariates)
  n_fixed <- ncol(X_cov)
  if (n_fixed > 0) {
    start_idx <- length(std_err) - n_fixed + 1
    if (start_idx <= length(std_err)) {
      std_fixed <- std_err[start_idx:length(std_err)]
    } else {
      std_fixed <- rep(NA, n_fixed)
    }
  }
  
  # Return all standard errors
  return(list(
    std_err = std_err,
    std_meas = std_meas,
    std_random = std_random,
    std_fixed = std_fixed,
    inv_fisher = inv_fisher
  ))
}

#' Process model fitting results to extract parameters and standard errors
#'
#' @param res Optimization result object
#' @param observed_fisher Observed Fisher information matrix
#' @param start_values Named starting values vector
#' @param estimate_params Logical vector indicating which parameters to estimate
#' @param model Original model object
#' @param model_options Model options with fixed parameter values
#' @param X_cov Covariate matrix
#' @param n_coeff_nonfixed Number of non-fixed coefficients
#' @return List containing extracted parameters and standard errors
#' @noRd
process_model_results <- function(res, observed_fisher, start_values, estimate_params, 
                                model, model_options, X_cov, n_coeff_nonfixed) {
  
  # Extract parameters
  param_results <- extract_and_organize_parameters(
    res = res,
    start_values = start_values,
    estimate_params = estimate_params,
    model = model,
    model_options = model_options,
    X_cov = X_cov,
    n_coeff_nonfixed = n_coeff_nonfixed
  )
  
  # Calculate standard errors
  se_results <- calculate_standard_errors(
    observed_fisher = observed_fisher,
    res = res,
    estimate_params = estimate_params,
    model = model,
    model_options = model_options,
    X_cov = X_cov,
    n_coeff_nonfixed = n_coeff_nonfixed,
    param_results = param_results
  )
  
  # Add "(fixed)" to parameter names for fixed parameters
  # For measurement error parameter (sigma_e)
  if (!estimate_params[1]) {
    names(param_results$coeff_meas) <- paste0(names(param_results$coeff_meas), " (fixed)")
  }
  
  # For random effects parameters
  if (length(param_results$coeff_random) > 0) {
    # Get parameter names from the parameter results
    param_names <- names(param_results$coeff_random)
    
    # For each parameter in coeff_random, check if it was estimated
    for (i in 1:length(param_names)) {
      param_name <- param_names[i]
      
      # Find if this parameter was estimated
      was_estimated <- FALSE
      
      # Check in standard param names
      if (param_name %in% names(estimate_params)) {
        was_estimated <- estimate_params[param_name]
      } 
      # Check for "Theta N" parameters
      else if (startsWith(param_name, "Theta ")) {
        # If theta was estimated or this was a non-stationary model with unnamed params
        if ("theta" %in% names(estimate_params)) {
          was_estimated <- estimate_params["theta"]
        } else {
          # Check unnamed parameters
          was_estimated <- any(names(estimate_params) == "" & estimate_params)
        }
      }
      
      if (!was_estimated) {
        # Rename the parameter to include "(fixed)"
        names(param_results$coeff_random)[i] <- paste0(param_name, " (fixed)")
        
        # Also rename in param_names for consistency
        param_results$par_names[i] <- paste0(param_name, " (fixed)")
      }
    }
  }
  
  # Combine results
  result <- c(
    param_results,
    list(
      std_err = se_results$std_err,
      std_meas = se_results$std_meas,
      std_random = se_results$std_random,
      std_fixed = se_results$std_fixed
    )
  )
  
  return(result)
}
#' @name intrinsic.operators
#' @title Covariance-based approximations of intrinsic fields
#' @description `intrinsic.matern.operators` is used for computing a
#' covariance-based rational SPDE approximation of intrinsic
#' fields on \eqn{R^d} defined through the SPDE
#' \deqn{(-\Delta)^{\beta/2} (\tau u) = \mathcal{W}}{(-\Delta)^{\beta/2} (\tau u) = \mathcal{W}}
#' @param tau precision parameter
#' @param beta Smoothness parameter
#' @param G The stiffness matrix of a finite element discretization
#' of the domain of interest.
#' @param C The mass matrix of a finite element discretization of
#' the domain of interest.
#' @param d The dimension of the domain.
#' @param mesh An inla mesh.
#' @param graph An optional `metric_graph` object. Replaces `d`, `C` and `G`.
#' @param loc_mesh locations for the mesh for `d=1`.
#' @param m The order of the rational approximation for the intrinsic part,
#' which needs to be a positive integer. The default value is 2.
#' @param compute_higher_order Logical. Should the higher order finite
#' element matrices be computed?
#' @param return_block_list Logical. For `type = "covariance"`,
#' should the block parts of the precision matrix be returned
#' separately as a list?
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are
#' "brasil", "chebfun" or "chebfunLB".
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @param scaling scaling factor, see details. 
#' @param opts options for numerical calulcation of the scaling, see details. 
#' @return `intrinsic.operators` returns an object of
#' class "intrinsicCBrSPDEobj". 
#' @export
#' @details The covariance operator
#' \deqn{\tau^{-2}(-\Delta)^{\beta}}{\tau^{-2}(-\Delta)^{\beta}}
#' is approximated based on a rational approximation. The Laplacian is 
#' equipped with homogeneous Neumann boundary
#' conditions and a zero-mean constraint is additionally imposed to obtained
#' a non-intrinsic model. The scaling is computed as the lowest positive eigenvalue of 
#' sqrt(solve(c0))%*%g1sqrt(solve(c0)). opts provides a list of options for the 
#' numerical calculation of the scaling factor, which is done using `Rspectra::eigs_sym`. 
#' See the help of that function for details. 
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   op <- intrinsic.operators(tau = 1, beta = beta, loc_mesh = x, d = 1)
#'   # Compute and plot the variogram of the model
#'   Sigma <- op$A[,-1] %*% solve(op$Q[-1,-1], t(op$A[,-1]))
#'   One <- rep(1, times = ncol(Sigma))
#'   D <- diag(Sigma)
#'   Gamma <- 0.5 * (One %*% t(D) + D %*% t(One) - 2 * Sigma)
#'   k <- 100
#'   plot(x, Gamma[k, ], type = "l")
#'   lines(x,
#'     variogram.intrinsic.spde(x[k], x, kappa = 0, alpha = 0, 
#'     beta = beta, L = 10, d = 1),
#'     col = 2, lty = 2
#'   )
#' }
intrinsic.operators <- function(tau = NULL,
                                beta = NULL,
                                G = NULL,
                                C = NULL,
                                d = NULL,
                                mesh = NULL,
                                graph = NULL,
                                loc_mesh = NULL,
                                m = 1,
                                compute_higher_order = FALSE,
                                return_block_list = FALSE,
                                type_rational_approximation = c(
                                           "brasil",
                                           "chebfun",
                                           "chebfunLB"
                                           ),
                                fem_mesh_matrices = NULL,
                                scaling = NULL,
                                opts = NULL) {
    
    if(is.null(tau) || is.null(beta)) {
        stop("tau and beta must be provided.")
    }
 return(intrinsic.matern.operators(tau = tau,
                                  kappa = 1, 
                                  alpha = 0, 
                                  beta = beta, 
                                  G = G,
                                  C = C,
                                  d = d,
                                  mesh = mesh,
                                  graph = graph,
                                  loc_mesh = loc_mesh,
                                  m_beta = m,
                                  compute_higher_order =  compute_higher_order,
                                  return_block_list = return_block_list,
                                  type_rational_approximation = type_rational_approximation,
                                  fem_mesh_matrices = fem_mesh_matrices,
                                  scaling = scaling,
                                  opts = opts))
}
#' @name intrinsic.operators.internal
#' @title Covariance-based approximations of intrinsic fields
#' @description `intrinsic.operators.internal` is used for computing a
#' covariance-based rational SPDE approximation of intrinsic
#' fields on \eqn{R^d} defined through the SPDE
#' \deqn{(-\Delta)^{\alpha/2}u = \mathcal{W}}{(-\Delta)^{\alpha/2}u = \mathcal{W}}
#' @param alpha Smoothness parameter
#' @param G The stiffness matrix of a finite element discretization
#' of the domain of interest.
#' @param C The mass matrix of a finite element discretization of
#' the domain of interest.
#' @param mesh An inla mesh.
#' @param d The dimension of the domain.
#' @param m The order of the rational approximation, which needs
#' to be a positive integer.
#' The default value is 2.
#' @param compute_higher_order Logical. Should the higher order finite
#' element matrices be computed?
#' @param return_block_list Logical. For `type = "covariance"`,
#' should the block parts of the precision matrix be returned
#' separately as a list?
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are
#' "brasil", "chebfun" or "chebfunLB".
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @param scaling scaling factor.
#' @param opts options for numerical calculation of scaling factor.
#' @return `intrinsic.operators.internal` returns an object of
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
#' \item{type}{String indicating the type of approximation.}
#' @details We use the covariance-based rational approximation of the
#' fractional operator. It is assumed that a mean-zero contraint is imposed
#' so that the equation has a unique solution. This contraint needs to be
#' imposed while working with the model later.
#' @noRd
intrinsic.operators.internal <- function(C,
                                G,
                                mesh,
                                alpha,
                                m = 2,
                                d,
                                compute_higher_order = FALSE,
                                return_block_list = FALSE,
                                type_rational_approximation = c(
                                  "brasil",
                                  "chebfun",
                                  "chebfunLB"
                                ),
                                fem_mesh_matrices = NULL,
                                scaling = NULL,
                                opts = NULL) {
  type_rational_approximation <- type_rational_approximation[[1]]

  if (is.null(fem_mesh_matrices)) {
    if (!is.null(mesh)) {
      d <- fmesher::fm_manifold_dim(mesh)
      m_alpha <- floor(alpha)
      m_order <- m_alpha + 1

      if (d > 1) {
        if (compute_higher_order) {
          # fem <- INLA::inla.mesh.fem(mesh, order = m_alpha + 1)
          fem <- fmesher::fm_fem(mesh, order = m_alpha + 1)
        } else {
          # fem <- INLA::inla.mesh.fem(mesh)
          fem <- fmesher::fm_fem(mesh)
        }

        C <- fem$c0
        G <- fem$g1
        Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
        GCi <- G %*% Ci
        CiG <- Ci %*% G
        Gk <- list()
        Gk[[1]] <- G
        for (i in 2:m_order) {
          Gk[[i]] <- fem[[paste0("g", i)]]
        }
      } else if (d == 1) {
        # fem <- INLA::inla.mesh.fem(mesh, order = 2)
        fem <- fmesher::fm_fem(mesh)
        C <- fem$c0
        G <- fem$g1
        Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
        GCi <- G %*% Ci
        CiG <- Ci %*% G
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
      m_alpha <- floor(alpha)
      m_order <- m_alpha + 1

      ## get lumped mass matrix
      C <- Matrix::Diagonal(dim(C)[1], rowSums(C))

      ## get G_k matrix: k is up to m_alpha if alpha is integer,
      # k is up to m_alpha + 1 otherwise.
      # inverse lumped mass matrix
      Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))

      GCi <- G %*% Ci
      CiG <- Ci %*% G
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
    CiG <- Ci %*% G
    m_alpha <- floor(alpha)
    m_order <- m_alpha + 1
    if (!is.null(mesh)) {
      d <- fmesher::fm_manifold_dim(mesh)
    }
    Gk <- list()
    Gk[[1]] <- G
    for (i in 2:m_order) {
      Gk[[i]] <- fem_mesh_matrices[[paste0("g", i)]]
    }
  }
  if (is.null(scaling)) {
      Cd <- Diagonal(dim(C)[1], 1/sqrt(diag(C)))
      Gg <- Cd%*%G%*%Cd
      if(is.null(opts)) { 
          opts = list(tol = 1e-10, maxitr = 1e4)
      }
      scaling <- RSpectra::eigs_sym(as(Gg, "CsparseMatrix"), 2, which = "SM",
                                    opts = list(tol = 1e-10, maxitr = 1e4))$values[1]
      if(is.na(scaling)){
          stop("Computation of scaling failed, provide the scaling manually or change opts to allow for higher maxitr or lower tol")
      }
  } else {
      if (!is.numeric(scaling) || length(scaling) != 1 || scaling <= 0) {
          stop("scaling must be a positive numeric value of length 1")
      }
  }

  if(alpha %% 1 == 0) {
      L <- G / scaling
      Q <- G
      if (alpha > 1) {
          for (k in 1:(alpha - 1)) {
              Q <- Q %*% CiG
          }
      }
      Q.frac <- Matrix::Diagonal(dim(G)[1])
      Q.int <- Q
      if (return_block_list) {
          Q <- list(Q)
      } else {
          Q.int <- list(Q.int = Q, order = m_alpha)
      }
  } else { # fractional case 
      L <- G / scaling
      
      CiL <- CiG / scaling
      
      if (m_alpha == 0) {
          aux_mat <- Diagonal(dim(L)[1])
      } else {
          aux_mat <- CiL
      }
      
      
      if (return_block_list) {
          Q.int <- aux_mat
          
          if (alpha %% 1 == 0) {
              Q.frac <- Matrix::Diagonal(dim(L)[1])
              Q <- G
              
              if (alpha > 1) {
                  for (k in 1:(alpha - 1)) {
                      Q <- Q %*% CiG
                  }
              }
              Q.int <- Q
          } else {
              if (m == 0) {
                  stop("Return block list does not work with m = 0, either increase m or set return_block_list to FALSE.")
              }
              Q.frac <- intrinsic.precision(
                  alpha = alpha, rspde.order = m, dim = d,
                  fem_mesh_matrices = fem_mesh_matrices, only_fractional = TRUE,
                  return_block_list = TRUE,
                  type_rational_approx = type_rational_approximation,
                  scaling = scaling
              )
              
              Q <- Q.frac
              
              if (m_alpha > 0) {
                  for (j in seq_len(length(Q))) {
                      for (i in 1:m_alpha) {
                          Q[[j]] <- Q[[j]] %*% Q.int
                      }
                  }
              }
              Q.int <- list(Q.int = Q.int, order = m_alpha)
          }
      } else {
          Q.int <- list(Q.int = kronecker(Diagonal(m + 1), aux_mat), order = m_alpha)
          
          if (alpha %% 1 == 0) {
              Q.frac <- Matrix::Diagonal(dim(G)[1])
              Q <- G
              
              if (alpha > 1) {
                  for (k in 1:(alpha - 1)) {
                      Q <- Q %*% CiL
                  }
              }
              
              Q.int <- list(Q.int = Q, order = m_alpha)
          } else if (m > 0) {
              Q.frac <- intrinsic.precision(
                  alpha = alpha,
                  rspde.order = m, dim = d,
                  fem_mesh_matrices = fem_mesh_matrices, only_fractional = TRUE,
                  type_rational_approx = type_rational_approximation,
                  scaling = scaling
              )
              
              Q <- Q.frac
              
              if (m_alpha > 0) {
                  for (i in 1:m_alpha) {
                      Q <- Q %*% Q.int$Q.int
                  }
              }
          } else {
              stop("m > 0 required for intrinsic fields")
          }
      }    
  }
  

  ## output
  output <- list(
    C = C, 
    G = G, 
    L = L, 
    Ci = Ci, 
    GCi = GCi, 
    Gk = Gk,
    fem_mesh_matrices = fem_mesh_matrices,
    alpha = alpha, 
    m = m, 
    d = d,
    Q.frac = Q.frac, 
    Q.int = Q.int,
    Q = Q, 
    sizeC = dim(C)[1],
    higher_order = compute_higher_order,
    type_rational_approximation = type_rational_approximation,
    return_block_list = return_block_list,
    stationary = TRUE,
    scaling = scaling
  )
  output$type <- "Covariance-Based intrinsic SPDE Approximation"
  class(output) <- "CBrSPDEobj"
  return(output)
}



intrinsic.precision <- function(alpha, rspde.order, dim, fem_mesh_matrices,
                                only_fractional = FALSE, return_block_list = FALSE,
                                type_rational_approx = "chebfun",
                                scaling = NULL) {
  n_m <- rspde.order
  
  mt <- get_rational_coefficients(n_m, type_rational_approx)


  m_alpha <- floor(alpha)

  row_nu <- round(1000 * cut_decimals(alpha))
  r <- unlist(mt[row_nu, 2:(1 + rspde.order)])
  p <- unlist(mt[row_nu, (2 + rspde.order):(1 + 2 * rspde.order)])
  k <- unlist(mt[row_nu, 2 + 2 * rspde.order])

  if (!only_fractional) {
    if (m_alpha == 0) {
      L <- fem_mesh_matrices[["g1"]] / scaling
      Q <- (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]
      if (length(r) > 1) {
        for (i in 2:length(r)) {
          Q <- bdiag(Q, (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i])
        }
      }
    } else {
      Malpha <- fem_mesh_matrices[[paste0("g", m_alpha)]] / scaling^m_alpha
      Malpha2 <- fem_mesh_matrices[[paste0("g", m_alpha + 1)]] / scaling^(m_alpha + 1)

      Q <- 1 / r[1] * (Malpha2 - p[1] * Malpha)

      if (length(r) > 1) {
        for (i in 2:length(r)) {
          Q <- bdiag(Q, 1 / r[i] * (Malpha2 - p[i] * Malpha))
        }
      }
    }


    # add k_part into Q

    if (m_alpha == 0) {
      C <- fem_mesh_matrices[["c0"]]
      Kpart <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
    } else {
      Kpart <- fem_mesh_matrices[[paste0("g", m_alpha)]] / scaling^m_alpha
    }
    Kpart <- Kpart / k

    Q <- bdiag(Q, Kpart)

    Q <- Q * scaling^(alpha)

    return(Q)
  } else {
    L <- fem_mesh_matrices[["g1"]] / scaling

    if (return_block_list) {
      Q <- list()

      Q[[length(Q) + 1]] <- scaling^alpha * (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]

      if (n_m > 1) {
        for (i in 2:(n_m)) {
          Q[[length(Q) + 1]] <- scaling^alpha *
            (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i]
        }
      }
      if (m_alpha == 0) {
        C <- fem_mesh_matrices[["c0"]]
        Kpart <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
      } else {
        Kpart <- fem_mesh_matrices[["c0"]]
      }
      Q[[length(Q) + 1]] <- scaling^alpha * Kpart / k

      return(Q)
    } else {
      Q <- (L - p[1] * fem_mesh_matrices[["c0"]]) / r[1]

      if (n_m > 1) {
        for (i in 2:(n_m)) {
          temp <- (L - p[i] * fem_mesh_matrices[["c0"]]) / r[i]
          Q <- bdiag(Q, temp)
        }
      }
      if (m_alpha == 0) {
        C <- fem_mesh_matrices[["c0"]]
        Kpart <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
      } else {
        Kpart <- fem_mesh_matrices[["c0"]]
      }

      Q <- bdiag(Q, Kpart / k)

      Q <- Q * scaling^alpha

      return(Q)
    }
  }
}

#' @name intrinsic.matern.operators
#' @title Covariance-based approximations of intrinsic fields
#' @description `intrinsic.matern.operators` is used for computing a
#' covariance-based rational SPDE approximation of intrinsic
#' fields on \eqn{R^d} defined through the SPDE
#' \deqn{(-\Delta)^{\beta/2}(\kappa^2-\Delta)^{\alpha/2} (\tau u) = \mathcal{W}}{(-\Delta)^{\beta/2}(\kappa^2-\Delta)^{\alpha/2} (\tau u) = \mathcal{W}}
#' @param kappa range parameter
#' @param tau precision parameter
#' @param alpha Smoothness parameter
#' @param beta Smoothness parameter
#' @param G The stiffness matrix of a finite element discretization
#' of the domain of interest.
#' @param C The mass matrix of a finite element discretization of
#' the domain of interest.
#' @param d The dimension of the domain.
#' @param mesh An inla mesh.
#' @param graph An optional `metric_graph` object. Replaces `d`, `C` and `G`.
#' @param loc_mesh locations for the mesh for `d=1`.
#' @param m_alpha The order of the rational approximation for the Matérn part,
#' which needs to be a positive integer. The default value is 2.
#' @param m_beta The order of the rational approximation for the intrinsic part,
#' which needs to be a positive integer. The default value is 2.
#' @param compute_higher_order Logical. Should the higher order finite
#' element matrices be computed?
#' @param return_block_list Logical. For `type = "covariance"`,
#' should the block parts of the precision matrix be returned
#' separately as a list?
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are
#' "brasil", "chebfun" or "chebfunLB".
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @param scaling scaling factor, see details. 
#' @param opts options for numerical calulcation of the scaling, see details. 
#' @return `intrinsic.matern.operators` returns an object of
#' class "intrinsicCBrSPDEobj". This object is a list including the
#' following quantities:
#' \item{C}{The mass lumped mass matrix.}
#' \item{Ci}{The inverse of `C`.}
#' \item{GCi}{The stiffness matrix G times `Ci`}
#' \item{Q}{The precision matrix.}
#' \item{Q_list}{A list containing the blocks required to assemble the precision matrix.}
#' \item{alpha}{The fractional power of the Matérn part of the operator.}
#' \item{beta}{The fractional power of the intrinsic part of the operator.}
#' \item{kappa}{Range parameter of the covariance function}
#' \item{tau}{Scale parameter of the covariance function.}
#' \item{m_alpha}{The order of the rational approximation for the Matérn part.}
#' \item{m_beta}{The order of the rational approximation for the intrinsic part.}
#' \item{m}{The total number of blocks in the precision matrix.}
#' \item{n}{The number of mesh nodes.}
#' \item{d}{The dimension of the domain.}
#' \item{type_rational_approximation}{String indicating the type of rational approximation.}
#' \item{higher_order}{Boolean indicating if higher order FEM-related matrices are computed.}
#' \item{return_block_list}{Boolean indicating if the precision matrix is returned as 
#' a list with the blocks.}
#' \item{fem_mesh_matrices}{A list containing the mass lumped mass
#' matrix, the stiffness matrix and
#' the higher-order FEM-related matrices.}
#' \item{make_A}{A function to compute the projection matrix which links the field to observation locations.}
#' \item{variogram}{A function to compute the variogram of the model at a specified node.}
#' \item{A}{Matrix that sums the components in the approximation to the mesh nodes.}
#' \item{scaling}{The scaling used in the intrinsic part of the model.}
#' @export
#' @details The covariance operator
#' \deqn{\tau^{-2}(-\Delta)^{\beta}(\kappa^2-\Delta)^{\alpha}}{\tau^{-2}(-\Delta)^{\beta}(\kappa^2-\Delta)^{\alpha}}
#' is approximated based on rational approximations of the two fractional
#' components. The Laplacians are equipped with homogeneous Neumann boundary
#' conditions. Unless supplied, the scaling is computed as the lowest positive eigenvalue of 
#' `sqrt(solve(c0))%*%g1%*%sqrt(solve(c0))`. opts provides a list of options for the 
#' numerical calculation of the scaling factor, which is done using `Rspectra::eigs_sym`. 
#' See the help of that function for details. 
#' 
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   kappa <- 1
#'   op <- intrinsic.matern.operators(
#'     kappa = kappa, tau = 1, alpha = alpha,
#'     beta = beta, loc_mesh = x, d = 1
#'   )
#'   # Compute and plot the variogram of the model
#'   Sigma <- op$A[,-1] %*% solve(op$Q[-1,-1], t(op$A[,-1]))
#'   One <- rep(1, times = ncol(Sigma))
#'   D <- diag(Sigma)
#'   Gamma <- 0.5 * (One %*% t(D) + D %*% t(One) - 2 * Sigma)
#'   k <- 100
#'   plot(x, Gamma[k, ], type = "l")
#'   lines(x,
#'     variogram.intrinsic.spde(x[k], x, kappa, alpha, beta, L = 10, d = 1),
#'     col = 2, lty = 2
#'   )
#' }
intrinsic.matern.operators <- function(kappa,
                                       tau,
                                       alpha,
                                       beta = 1,
                                       G = NULL,
                                       C = NULL,
                                       d = NULL,
                                       mesh = NULL,
                                       graph = NULL,
                                       loc_mesh = NULL,
                                       m_alpha = 2,
                                       m_beta = 2,
                                       compute_higher_order = FALSE,
                                       return_block_list = FALSE,
                                       type_rational_approximation = c(
                                           "brasil", "chebfun",
                                         "chebfunLB"
                                       ),
                                       fem_mesh_matrices = NULL,
                                       scaling = NULL,
                                       opts = NULL) {
  if (is.null(d) && is.null(mesh) && is.null(graph)) {
    stop("You should give either the dimension d, the mesh or graph!")
  }

  if ((is.null(C) || is.null(G)) && is.null(mesh) && is.null(graph) && (is.null(loc_mesh) || d != 1)) {
    stop("You should either provide mesh, graph, or provide both C *and* G!")
  }

  if ((is.null(C) || is.null(G)) && (is.null(graph)) && (!is.null(loc_mesh) && d == 1)) {
    fem <- rSPDE.fem1d(loc_mesh)
    C <- fem$C
    G <- fem$G

    fem_mesh_matrices <- list()
    fem_mesh_matrices[["c0"]] <- C
    fem_mesh_matrices[["g1"]] <- G

    Gk <- list()
    Gk[[1]] <- G

    m_order <- floor(alpha) + 1

    if (compute_higher_order) {
      for (i in 1:m_order) {
        fem_mesh_matrices[[paste0("g", i)]] <- Gk[[i]]
      }
    }
  }

  has_mesh <- FALSE
  has_graph <- FALSE

  if (!is.null(loc_mesh)) {
    if (!is.numeric(loc_mesh)) {
      stop("loc_mesh must be numerical.")
    }
  }

  if (!is.null(graph)) {
    if (!inherits(graph, "metric_graph")) {
      stop("graph should be a metric_graph object!")
    }
    d <- 1
    if (is.null(graph$mesh)) {
      warning("The graph object did not contain a mesh, one was created with h = 0.01. Use the build_mesh() method to replace this mesh with a new one.")
      graph$build_mesh(h = 0.01)
    }
    graph$compute_fem()
    C <- graph$mesh$C
    G <- graph$mesh$G
    has_graph <- TRUE
  }

  if (!is.null(mesh)) {
    d <- fmesher::fm_manifold_dim(mesh)
    # fem <- INLA::inla.mesh.fem(mesh)
    fem <- fmesher::fm_fem(mesh)
    C <- fem$c0
    G <- fem$g1
    has_mesh <- TRUE
  }



  kappa <- rspde_check_user_input(kappa, "kappa", 0)
  tau <- rspde_check_user_input(tau, "tau", 0)

  alpha <- rspde_check_user_input(alpha, "alpha", 0)
  alpha <- min(alpha, 10)

  #beta <- rspde_check_user_input(beta, "beta", lower_bound = 0, upper_bound = 1+d/2)
  beta <- rspde_check_user_input(beta, "beta", lower_bound = 0)

  if (alpha + beta < d / 2) {
    stop("One must have alpha + beta > d/2")
  }

  if (alpha > 0 && beta > 0 && kappa > 0) {
    op1 <- CBrSPDE.matern.operators(
      C = C, G = G, mesh = mesh, nu = alpha - d / 2, kappa = kappa, tau = tau,
      fem_mesh_matrices = fem_mesh_matrices,
      m = m_alpha, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE,
      type_rational_approximation = type_rational_approximation[[1]]
    )
    
    if (is.list(op1$Q)) {
        Q.list1 <- op1$Q
    } else {
        Q.list1 <- list(op1$Q)
    }
    
    op2 <- intrinsic.operators.internal(
      C = C, G = G, mesh = mesh, alpha = beta,
      m = m_beta, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE, fem_mesh_matrices = fem_mesh_matrices,
      type_rational_approximation = type_rational_approximation[[1]],
      scaling = scaling
    )
    
    if (is.list(op2$Q)) {
        Q.list2 <- op2$Q
    } else {
        Q.list2 <- list(op2$Q)
    }
    
    scaling <- op2$scaling
    block_list <- list()
    
    
    m1 <- length(Q.list1)
    m2 <- length(Q.list2)
    k <- 1
    if (return_block_list) {
      Q <- list()
    }
    Q_list <- list(Qproper = Q.list1,
                   Qintrinsic = Q.list2)
    for (i in 1:m1) {
      for (j in 1:m2) {
        if (return_block_list) {
          Q[[k]] <- Q.list1[[i]] %*% op1$Ci %*% Q.list2[[j]]
        } else {
          if (i == 1 && j == 1) {
            Q <- Q.list1[[i]] %*% op1$Ci %*% Q.list2[[j]]
          } else {
            Q <- bdiag(Q, Q.list1[[i]] %*% op1$Ci %*% Q.list2[[j]])
          }
        }
        k <- k + 1
      }
    }
    n <- dim(op1$C)[1]
    A <- kronecker(matrix(rep(1, m1 * m2), 1, m1 * m2), Diagonal(n))    
  } else if (alpha > 0 && kappa > 0) {
    op1 <- CBrSPDE.matern.operators(
      C = C, G = G, mesh = mesh, nu = alpha - d / 2, kappa = kappa, tau = tau,
      fem_mesh_matrices = fem_mesh_matrices,
      m = m_alpha, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE,
      type_rational_approximation = type_rational_approximation[[1]]
    )
    if (is.list(op1$Q)) {
      Q.list1 <- op1$Q
    } else {
      Q.list1 <- list(op1$Q)
    }
    m1 <- length(Q.list1)
    if (!return_block_list) {
      for (i in 1:m1) {
        if (i == 1) {
          Q <- Q.list1[[i]]
        } else {
          Q <- bdiag(Q, Q.list1[[i]])
        }
      }
    } else {
      Q <- Q.list1
    }
    Q_list <- list(Qproper = Q.list1,
                   Qintrinsic = NULL)
    
    n <- dim(op1$C)[1]
    A <- kronecker(matrix(rep(1, m1), 1, m1), Diagonal(n))
  } else if (beta > 0) {
    if (kappa == 0) {
      alpha_beta <- alpha + beta
    } else {
      alpha_beta <- beta
    }
    op1 <- intrinsic.operators.internal(
      C = C, G = G, mesh = mesh, alpha = alpha_beta,
      m = m_beta, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE, fem_mesh_matrices = fem_mesh_matrices,
      type_rational_approximation = type_rational_approximation[[1]],
      scaling = scaling
    )
    scaling <- op1$scaling
    if (is.list(op1$Q)) {
      Q.list1 <- lapply(op1$Q, function(x) { x * tau^2})
    } else {
      Q.list1 <- list(op1$Q*tau^2)
    }
    m1 <- length(Q.list1)
    if (!return_block_list) {
      for (i in 1:m1) {
        if (i == 1) {
          Q <- Q.list1[[i]]
        } else {
          Q <- bdiag(Q, Q.list1[[i]])
        }
      }
    } else {
      Q <- Q.list1
    }
    Q_list <- list(Qintrinsic = Q.list1,
                   Qproper = NULL)
    n <- dim(op1$C)[1]

    A <- kronecker(matrix(rep(1, m1), 1, m1), Diagonal(n))    
  }

  n <- dim(op1$C)[1]
  m1 <- max(c(length(Q_list$Qproper),1))
  m2 <- max(c(length(Q_list$Qintrinsic),1))
  m <- m1*m2
  if (!is.null(mesh)) {
      make_A <- function(loc) {
        Ai <- fmesher::fm_basis(x = mesh, loc = loc)
        A <- kronecker(matrix(rep(1, m), 1, m), Ai)   
        return(A)
      }
  } else if (!is.null(graph)) {
      make_A <- function(loc) {
          Ai <- graph$fem_basis(loc)
          A <- kronecker(matrix(rep(1, m), 1, m), Ai)   
          return(A)
      }
  } else if (!is.null(loc_mesh) && d == 1) {
      make_A <- function(loc) {
          Ai <- rSPDE::rSPDE.A1d(x = loc_mesh, loc = loc)
          A <- kronecker(matrix(rep(1, m), 1, m), Ai)   
          return(A)
      }
  } else {
      make_A <- NULL
  }
  
  variogram <- function(loc, semi = FALSE) {
    if(return_block_list) { 
        QQ <- Q[[1]]
        if(m>1) {
            for(i in 2:m){
                QQ <- bdiag(QQ,Q[[i]])
            }
        }
    } else {
        QQ <- Q
    }
    if(beta < 1) {
      Sigma <- A%*% solve(QQ, t(A))
    } else {
      ind.fix <- 1 + seq(from=0,to= m * n, by = n)
      Sigma <- A[,-ind.fix] %*% solve(QQ[-ind.fix,-ind.fix], t(A[,-ind.fix]))    
    }
    
    One <- rep(1, times = ncol(Sigma))
    D <- diag(Sigma)
    
    if(!is.null(mesh)) {
      Ai <- fmesher::fm_basis(x = mesh, loc = loc)    
    } else if (!is.null(graph)) {
    Ai <- graph$fem_basis(loc)
    } else if (!is.null(loc_mesh) && d == 1) {
    Ai <- rSPDE.A1d(x = loc_mesh, loc = loc)
    }
    Gamma <- Ai %*% (One %*% t(D) + D %*% t(One) - 2 * Sigma)
    
    if(semi) {
      return(0.5 * Gamma) 
    } else {
      return(Gamma)
    }
  }
  mean_correction <- function(full = FALSE, index = 1) {
      if(index < 1 || index > n) {
          stop("Index out of bounds.")
      }
      if(!full) {
          out <- rep(0,n)    
      }
      
      for(i in 1:m) {
          if(return_block_list) { 
              QQ <- Q[[i]][-index,-index]
          } else {
              ind <- setdiff((1+n*(i-1)) : (n*i), n*(i-1) + index)
              QQ <- Q[ind,ind]
          }
          
          
          if(full) {
              vec <- rep(0,n)
              tryCatch(
                  expr = {
                      vec[-index] <- -diag(MetricGraph::selected_inv(QQ))/2},
                  error = function(e) {
                      vec[-index] = -diag(solve(QQ))/2
                  }
              ) 
              if(i == 1) {
                  out <- vec 
              } else {
                  out <- c(out,vec)    
              }
          } else {
              tryCatch(
                  expr = {out[-index] <- out[-index] - diag(MetricGraph::selected_inv(QQ))/2},
                  error = function(e) {
                      out[-index] <- out[-index] - diag(solve(QQ))/2
                  }
              )    
          }
          
          
      }
      return(out)
  }
  out <- list(
    C = op1$C, 
    G = op1$G, 
    Ci = op1$Ci, 
    GCi = op1$GCi,
    Q = Q,
    Q_list = Q_list,
    alpha = alpha, 
    beta = beta, 
    kappa = kappa, 
    tau = tau,
    m_alpha = m_alpha, 
    m_beta = m_beta, 
    m = m,
    n = dim(C)[1],
    d = d,
    type_rational_approximation = type_rational_approximation[[1]],
    higher_order = compute_higher_order,
    return_block_list = return_block_list,
    stationary = TRUE,
    has_mesh = has_mesh,
    has_graph = has_graph,
    make_A = make_A,
    variogram = variogram,
    mean_correction = mean_correction,
    A = A,
    mesh = mesh,
    graph = graph,
    loc_mesh = loc_mesh,
    scaling = scaling
  )
  out$type <- "Covariance-Based intrinsic Matern SPDE Approximation"
  class(out) <- "intrinsicCBrSPDEobj"

  return(out)
}



#' @name simulate.intrinsicCBrSPDEobj
#' @title Simulation of a fractional intrinsic SPDE using the
#' covariance-based rational SPDE approximation
#' @description The function samples a Gaussian random field based using the
#' covariance-based rational SPDE approximation.
#' @param object The covariance-based rational SPDE approximation,
#' computed using [intrinsic.matern.operators()]
#' @param nsim The number of simulations.
#' @param seed An object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param kappa new value of kappa to use
#' @param tau new value of tau to use
#' @param alpha new value of alpha to use
#' @param beta new value of beta to use
#' @param integral.constraint Should the contraint on the integral be done?
#' @param use_kl Simulate based on a KL expansion?
#' @param ... Currently not used.
#' @return A matrix with the `nsim` samples as columns.
#' @method simulate intrinsicCBrSPDEobj
#' @export
simulate.intrinsicCBrSPDEobj <- function(object, nsim = 1, 
                                         seed = NULL,
                                         kappa = NULL,
                                         tau = NULL,
                                         alpha = NULL,
                                         beta = NULL,
                                         integral.constraint = TRUE, 
                                         use_kl = NULL, ...) {
    
    if(!is.null(kappa) || !is.null(tau) || !is.null(alpha) || !is.null(beta)) {
        object <- update(object, kappa = kappa, tau = tau, 
                         alpha = alpha, beta = beta)    
    }
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if(is.null(use_kl)) {
        if(object$alpha == 0) {
            use_kl <- FALSE
        } else {
            use_kl <- TRUE
        }
    }
    
    if (object$return_block_list) {
        n <- object$n
        m <- object$m
        
        X <- matrix(0,n,nsim)
        for(i in 1:m){
            if(object$beta < 1 && object$alpha == 0) {
                Q <- object$Q[[i]]
                Z <- rnorm(n * nsim)
                dim(Z) <- c(n, nsim)
                LQ <-  Matrix::Cholesky(forceSymmetric(Q))
                X <- X + as.matrix(solve(LQ, Z))  
            } else {
                if(!use_kl) {
                    Q <- object$Q[[i]]
                    R <- Cholesky(Q, LDL = TRUE)
                    tmp <- expand2(R)
                    Z <- rnorm(n * nsim)
                    dim(Z) <- c(n, nsim)
                    Di <- Diagonal(n,c(1/sqrt(diag(tmp$D)[-n]),0))
                    X <- X + as.matrix(solve(tmp$P1, solve(t(tmp$L1), Di%*%Z)))
                } else {
                    ev <- eigen(object$Q[[i]])
                    V <- ev$vectors[,-n]
                    lambda <- ev$values[-n]
                    for(j in 1:(n-1)) {
                        Z <- kronecker(matrix(rnorm(nsim),1,nsim), rep(1,n))
                        X <- X + Z * kronecker(matrix(rep(1,nsim),1,nsim), V[,j]) / sqrt(lambda[j])
                    }
                }
            }
        }
    } else {
      n <- object$n
      m <- object$m
      
      X <- matrix(0,n,nsim)
      for(i in 1:m){
          if(object$beta < 1 && object$alpha == 0) {
              ind <- (1+n*(i-1)) : (n*i)
              Q <- object$Q[ind,ind]
              Z <- rnorm(n * nsim)
              dim(Z) <- c(n, nsim)
              LQ <-  chol(forceSymmetric(Q))
              X <- X + as.matrix(solve(LQ, Z))      
          } else {
              if(!use_kl) {
                  ind <- (1+n*(i-1)) : (n*i)
                  Q <- object$Q[ind,ind]
                  R <- Cholesky(Q, LDL = TRUE)
                  tmp <- expand2(R)
                  Z <- rnorm(n * nsim)
                  dim(Z) <- c(n, nsim)
                  Di <- Diagonal(n,c(1/sqrt(diag(tmp$D)[-n]),0))
                  X <- X + as.matrix(solve(tmp$P1, solve(t(tmp$L1), Di%*%Z)))
              } else {
                  ind <- (1+n*(i-1)) : (n*i)
                  ev <- eigen(object$Q[ind,ind])
                  V <- ev$vectors[,-n]
                  lambda <- ev$values[-n]
                  for(j in 1:(n-1)) {
                      Z <- kronecker(matrix(rnorm(nsim),1,nsim), rep(1,n))
                      X <- X + Z * kronecker(matrix(rep(1,nsim),1,nsim), V[,j]) / sqrt(lambda[j])
                  }
              }
          }  
        }
    }
    #Add zero integral constraint
    if(integral.constraint){
        h <- diag(object$C)
        for(i in 1:nsim){
            X[,i] <- X[,i] - sum(h*X[,i])/sum(h)
        }    
    }
    
    return(X)
}

#' @name update.intrinsicCBrSPDEobj
#' @title Update parameters of intrinsicCBrSPDEobj objects
#' @description Function to change the parameters of a intrinsicCBrSPDEobj object
#' @param object Model object created by [intrinsic.matern.operators()]
#' @param kappa kappa value to be updated.
#' @param tau tau value to be update.
#' @param alpha alpha value to be updated.
#' @param beta beta value to be updated. .
#' @param ... currently not used.
#'
#' @return An object of type intrinsicCBrSPDEobj with updated parameters.
#' @export
#' @method update intrinsicCBrSPDEobj
#'
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   kappa <- 1
#'   op <- intrinsic.matern.operators(
#'     kappa = kappa, tau = 1, alpha = alpha,
#'     beta = beta, loc_mesh = x, d = 1
#'   )
#'op <- update(op, beta = 1.1, alpha = 0.9) 
#'}
update.intrinsicCBrSPDEobj <- function(object, 
                                       kappa = NULL,
                                       tau = NULL,
                                       alpha = NULL,
                                       beta = NULL,
                                       ...) {
    
    if (!is.null(kappa)) {
        kappa <- rspde_check_user_input(kappa, "kappa", 0)
    } else {
        kappa <- object$kappa
    }
    
    if (!is.null(tau)) {
        tau <- rspde_check_user_input(tau, "tau", 0)
    } else {
        tau <- object$tau
    }
    
    if(!is.null(beta)){
        beta <- rspde_check_user_input(beta, "beta", lower_bound = 0)
        
    } else {
        beta <- object$beta
    }
    if(!is.null(alpha)){
        alpha <- rspde_check_user_input(alpha, "alpha", lower_bound = 0)
    } else {
        alpha <- object$alpha
    }
    
    if(alpha + beta < object$d/2){
        stop("One must have alpha + beta > d/2")
    }
    
    return(intrinsic.matern.operators(kappa = kappa,
                                      tau = tau,
                                      alpha = alpha,
                                      beta = beta,
                                      G = object$G,
                                      C = object$C,
                                      d = object$d,
                                      mesh = object$mesh,
                                      graph = object$graph,
                                      loc_mesh = object$loc_mesh,
                                      m_alpha = object$m_alpha,
                                      m_beta = object$m_beta,
                                      compute_higher_order = object$higher_order,
                                      return_block_list = object$return_block_list,
                                      type_rational_approximation = object$type_rational_approximation,
                                      scaling = object$scaling))
}


#' @name precision.intrinsicCBrSPDEobj
#' @title Get the precision matrix of intrinsicCBrSPDEobj objects
#' @description Function to get the precision matrix of a intrinsicCBrSPDEobj object
#' @param object The model object computed using [intrinsic.matern.operators()]
#' @param kappa If non-null, update the range parameter.
#' @param tau If non-null, update the precision parameter.
#' @param alpha If non-null, update the alpha parameter.
#' @param beta If non-null, update the beta parameter.
#' @param ld If TRUE, return the log determinant of the precision matrix instead 
#' of the precision matrix. By default FALSE.
#' @param ... Currently not used.
#' @return The precision matrix.
#' @method precision intrinsicCBrSPDEobj
#' @seealso [simulate.intrinsicCBrSPDEobj()], [intrinsic.matern.operators()]
#' @export
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   kappa <- 1
#'   op <- intrinsic.matern.operators(
#'     kappa = kappa, tau = 1, alpha = alpha,
#'     beta = beta, loc_mesh = x, d = 1
#'   )
#' Q <- precision(op) 
#'}
precision.intrinsicCBrSPDEobj <- function(object,
                                   kappa = NULL,
                                   tau = NULL,
                                   alpha = NULL,
                                   beta = NULL,
                                   ld = FALSE,
                                   ...) {
    object <- update.intrinsicCBrSPDEobj(
        object = object,
        kappa = kappa,
        tau = tau,
        alpha = alpha,
        beta = beta
    )
    if(ld) {
        n <- dim(object$C)[1]
        if(!is.null(object$Q_list$Qproper) && !is.null(object$Q_list$Qintrinsic)) {
            m1 <- length(object$Q_list$Qproper)
            m2 <- length(object$Q_list$Qintrinsic)
            detCi <- sum(log(diag(object$Ci)))
            det1 <- rep(1,m1)
            det2 <- rep(1,m2)
            for(i in 1:m1) {
                Q.R <- Matrix::Cholesky(object$Q_list$Qproper[[i]])
                det1[i] <- 2 * c(determinant(Q.R, logarithm = TRUE, sqrt = TRUE)$modulus)
            }
            for(i in 1:m2) {
                Q.R <- Matrix::Cholesky(object$Q_list$Qintrinsic[[i]][-1,-1])
                det2[i] <- 2 * c(determinant(Q.R, logarithm = TRUE, sqrt = TRUE)$modulus)
            }
            logQ <- m1*m2*detCi + m1*(sum(det2) + m2*log(n)) + m2*sum(det1)
        } else if (!is.null(object$Q_list$Qintrinsic)){
            m2 <- length(object$Q_list$Qintrinsic)
            det2 <- rep(1,m2)
            for(i in 1:m2) {
                if(object$beta < 1 && object$alpha == 0) {
                #if(0) {
                    Q.R <- Matrix::Cholesky(object$Q_list$Qintrinsic[[i]])
                    det2[i] <- 2 * c(determinant(Q.R, logarithm = TRUE, sqrt = TRUE)$modulus)
                } else {
                    Q.R <- Matrix::Cholesky(object$Q_list$Qintrinsic[[i]][-1,-1])    
                    det2[i] <- 2 * c(determinant(Q.R, logarithm = TRUE, sqrt = TRUE)$modulus) + log(n)
                }
            }
            logQ <- sum(det2) 
        } else {
            m1 <- length(object$Q_list$Qproper)
            det1 <- rep(1,m1)
            for(i in 1:m1) {
                Q.R <- Matrix::Cholesky(object$Q_list$Qproper[[i]])
                det1[i] <- 2 * c(determinant(Q.R, logarithm = TRUE, sqrt = TRUE)$modulus)
            }
            logQ <- sum(det1)
        }
        return(logQ)
    } else {
        return(object$Q)    
    }
}


#' @name predict.intrinsicCBrSPDEobj
#' @title Prediction of an intrinsic Whittle-Matern model
#' @description The function is used for computing kriging predictions based
#' on data \eqn{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon}{\epsilon}
#' is mean-zero Gaussian measurement noise and \eqn{u(s)}{u(s)} is defined by
#' an intrinsic SPDE as described in [intrinsic.matern.operators()].
#' @param object The covariance-based rational SPDE approximation,
#' computed using [intrinsic.matern.operators()]
#' @param A A matrix linking the measurement locations to the basis of the FEM
#' approximation of the latent model.
#' @param Aprd A matrix linking the prediction locations to the basis of the
#' FEM approximation of the latent model.
#' @param Y A vector with the observed data, can also be a matrix where the
#' columns are observations of independent replicates of \eqn{u}.
#' @param sigma.e The standard deviation of the Gaussian measurement noise.
#' Put to zero if the model does not have measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param compute.variances Set to also TRUE to compute the kriging variances.
#' @param posterior_samples If `TRUE`, posterior samples will be returned.
#' @param n_samples Number of samples to be returned. Will only be used if `sampling` is `TRUE`.
#' @param only_latent Should the posterior samples be only given to the laten model?
#' @param ... further arguments passed to or from other methods.
#' @return A list with elements
#' \item{mean }{The kriging predictor (the posterior mean of u|Y).}
#' \item{variance }{The posterior variances (if computed).}
#' @export
#' @method predict intrinsicCBrSPDEobj
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   kappa <- 1
#'   op <- intrinsic.matern.operators(
#'     kappa = kappa, tau = 1, alpha = alpha,
#'     beta = beta, loc_mesh = x, d = 1
#'   )
#' # Create some data
#' u <-  simulate(op)
#' sigma.e <- 0.1
#' obs.loc <- runif(n = 20, min = 0, max = 10)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y <- as.vector(A %*% u + sigma.e * rnorm(20))
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
#'}
predict.intrinsicCBrSPDEobj <- function(object, 
                                        A, 
                                        Aprd, 
                                        Y, 
                                        sigma.e, 
                                        mu = 0,
                                        compute.variances = FALSE, 
                                        posterior_samples = FALSE,
                                        n_samples = 100, only_latent = FALSE,
                                        ...) {
    Y <- as.matrix(Y)
    if (dim(Y)[1] != dim(A)[1]) {
        stop("the dimensions of A does not match the number of observations")
    }
    
    n <- dim(Y)[1]
    out <- list()
    
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
    
    if(length(mu) == 1) {
        mu <- rep(mu, dim(object$Q)[1])
    } else {
        if(length(mu) == dim(object$C)[1] && dim(object$C)[1] < dim(object$Q)[1]) {
            mu <- rep(mu,object$m) / object$m
        } else if (length(mu) != dim(object$Q)[1]) {
            stop("the length of mu is wrong.")
        }
    }
    
    if (!no_nugget) {
        ## construct Q
        Q <- object$Q
        ## compute Q_x|y
        Q_xgiveny <- (t(A) %*% Q.e %*% A) + Q
        
        ## construct mu_x|y
        mu_xgiveny <- t(A) %*% Q.e %*% (Y - A%*%mu)    
        
        R <- Matrix::Cholesky(forceSymmetric(Q_xgiveny))
        mu_xgiveny <- solve(R, mu_xgiveny, system = "A")
        
        mu_xgiveny <- mu + mu_xgiveny
        out$mean <- Aprd %*% mu_xgiveny
        
        if (compute.variances) {
            out$variance <- diag(Aprd %*% solve(R, t(Aprd), system = "A"))
        }
    } else {
        Q <- object$Q
        
        QiAt <- solve(Q, t(A))
        AQiA <- A %*% QiAt
        xhat <- solve(Q, t(A) %*% solve(AQiA, Y))
        
        out$mean <- as.vector(Aprd %*% xhat)
        if (compute.variances) {
            M <- Q - QiAt %*% solve(AQiA, t(QiAt))
            out$variance <- diag(Aprd %*% M %*% t(Aprd))
        }
    }
    
    
    if (posterior_samples) {
        if (!no_nugget) {
            Z <- rnorm(dim(object$Q)[1] * n_samples)
            dim(Z) <- c(dim(object$Q)[1], n_samples)
            LQ <-  chol(forceSymmetric(Q_xgiveny))
            X <- as.matrix(solve(LQ, Z)) + kronecker(as.matrix(mu_xgiveny), 
                                                     matrix(rep(1,n_samples),1,n_samples))
            X <- Aprd %*% X
            if (!only_latent) {
                X <- X + matrix(rnorm(n_samples * dim(Aprd)[1], sd = sigma.e), nrow = dim(Aprd)[1])
            }
            return(X)
        } else {
            M <- Q - QiAt %*% solve(AQiA, t(QiAt))
            post_cov <- Aprd %*% M %*% t(Aprd)
            Y_tmp <- as.matrix(Y)
            mean_tmp <- as.matrix(out$mean)
            out$samples <- lapply(1:ncol(Y_tmp), function(i) {
                Z <- rnorm(dim(post_cov)[1] * n_samples)
                dim(Z) <- c(dim(post_cov)[1], n_samples)
                LQ <-  Matrix::Cholesky(forceSymmetric(post_cov))
                X <- LQ %*% Z
                X <- X + mean_tmp[, i]
                if (!only_latent) {
                    X <- X + matrix(rnorm(n_samples * dim(Aprd)[1], sd = sigma.e), nrow = dim(Aprd)[1])
                }
                return(X)
            })
        }
    }
    return(out)
}

#' @name intrinsic.loglike
#' @title Object-based log-likelihood function for latent intrinsic SPDE model 
#' @description This function evaluates the log-likelihood function for an
#' intrinsic Whittle-Matern SPDE model, that is observed under
#' Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables.
#' @param object The model object computed using [intrinsic.matern.operators()]
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location
#' to the finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param kappa If non-null, update the range parameter.
#' @param tau If non-null, update the precision parameter.
#' @param alpha If non-null, update the alpha parameter.
#' @param beta If non-null, update the beta parameter.
#' @return The log-likelihood value.
#' @noRd
#' @seealso [intrinsic.matern.operators()], [predict.intrinsicCBrSPDEobj()]
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   kappa <- 1
#'   op <- intrinsic.matern.operators(
#'     kappa = kappa, tau = 1, alpha = alpha,
#'     beta = beta, loc_mesh = x, d = 1
#'   )
#' # Create some data
#' obs.loc <- runif(n = 20, min = 0, max = 10)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y <- as.vector(A %*% u + sigma.e * rnorm(20))
#' loglik <- intrinsic.loglike(object, Y, A, sigma.e)
#'}
intrinsic.loglike <- function(object, 
                              Y, 
                              A, 
                              sigma.e, 
                              mu = 0,
                              kappa = NULL,
                              tau = NULL,
                              alpha = NULL,
                              beta = NULL) {
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
    
    if(!is.null(kappa) || !is.null(tau) || !is.null(alpha) || !is.null(beta)) {
        object <- update.intrinsicCBrSPDEobj(
            object = object,
            kappa = kappa,
            tau = tau,
            alpha = alpha,
            beta = beta)    
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
    
    Q <- object$Q
    #compute log determinant of Q
    logQ <- precision(object, ld = TRUE)
    
    ## compute Q_x|y
    Q <- object$Q
    
    Q_xgiveny <- t(A) %*% Q.e %*% A + Q
    ## construct mu_x|y
    
    mu_xgiveny <- t(A) %*% Q.e %*% Y
    # upper triangle with reordering
    
    
    R <- Matrix::Cholesky(Q_xgiveny)
    
    mu_xgiveny <- solve(R, mu_xgiveny, system = "A")
    
    mu_xgiveny <- mu + mu_xgiveny
    
    ## compute log|Q_xgiveny|
    ind <- 1 + seq(from = 0, to = (object$m-1)*object$n, by = object$n)
    R <- Matrix::Cholesky(forceSymmetric(Q_xgiveny[-ind,-ind]))
    log_Q_xgiveny <- 2 * determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus
    
    ## compute mu_x|y*Q*mu_x|y
    if (n.rep > 1) {
        mu_part <- sum(colSums((mu_xgiveny - mu) * (Q %*% (mu_xgiveny - mu))))
    } else {
        mu_part <- t(mu_xgiveny - mu) %*% Q %*% (mu_xgiveny - mu)
    }
    ## compute central part
    if (n.rep > 1) {
        central_part <- sum(colSums((Y - A %*% mu_xgiveny) * (Q.e %*% (Y - A %*% mu_xgiveny))))
    } else {
        central_part <- t(Y - A %*% mu_xgiveny) %*% Q.e %*% (Y - A %*% mu_xgiveny)
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

#' @noRd
aux_lme_intrinsic.loglike <- function(object, y, X_cov, repl, A_list, sigma_e, 
                                      beta_cov, mean_correction) {
    l_tmp <- tryCatch(
        aux2_lme_intrinsic.loglike(
            object = object,
            y = y, X_cov = X_cov, repl = repl, A_list = A_list,
            sigma_e = sigma_e, beta_cov = beta_cov, mean_correction
        ),
        error = function(e) {
            return(NULL)
        }
    )
    if (is.null(l_tmp)) {
        return(-10^100)
    }
    return(l_tmp)
}

#' @noRd
aux2_lme_intrinsic.loglike <- function(object, y, X_cov, repl, A_list, sigma_e, 
                                       beta_cov, mean_correction) {
    
    #compute prior log determinant
    Q <- object$Q
    prior.ld <- 0.5*precision(object, ld = TRUE)
    if(mean_correction) {
        mean_latent <- object$mean_correction()
    } else {
        mean_latent <- rep(0,dim(Q)[1])
    }
    repl_val <- unique(repl)
    
    l <- 0
    
    for (i in repl_val) {
        ind_tmp <- (repl %in% i)
        y_tmp <- y[ind_tmp]
        
        if (ncol(X_cov) == 0) {
            X_cov_tmp <- 0
        } else {
            X_cov_tmp <- X_cov[ind_tmp, , drop = FALSE]
        }
        
        na_obs <- is.na(y_tmp)
        
        y_ <- y_tmp[!na_obs]
        
        n.o <- length(y_)
        A_tmp <- A_list[[as.character(i)]]
        Q.p <- Q + t(A_tmp) %*% A_tmp / sigma_e^2
        
        if(object$beta < 1 && object$alpha == 0) {
        #if(0) {
            posterior.ld <- 0.5*c(determinant(Q.p, logarithm = TRUE)$modulus)
        } else {
            if(object$m == 1) {
            #    if(0) {
                posterior.ld <- 0.5*c(determinant(Q.p, logarithm = TRUE)$modulus)
            } else {
                ind <- 1 + seq(from = object$n, 
                               to = (object$m-1)*object$n, by = object$n)
                posterior.ld <- 0.5*(c(determinant(Q.p[-ind,-ind], logarithm = TRUE)$modulus) + (object$m-1)*log(dim(Q)[1]) - log(object$m))  
                l <- l - log(2*pi)    
            }
        }
        
        l <- l + prior.ld - posterior.ld - n.o * log(sigma_e)
        
        v <- y_
        
        if (ncol(X_cov) > 0) {
            X_cov_tmp <- X_cov_tmp[!na_obs, , drop = FALSE]
            # X_cov_tmp <- X_cov_list[[as.character(i)]]
            v <- v - X_cov_tmp %*% beta_cov - A_tmp %*% mean_latent
        }
        R.p <- Matrix::Cholesky(Q.p)
        mu.p <- solve(R.p, as.vector(t(A_tmp) %*% v / sigma_e^2), system = "A")
        
        v <- v - A_tmp %*% mu.p
        
        l <- l - 0.5 * (t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
            0.5 * n.o * log(2 * pi)
    }
    
    return(as.double(l))
}



#' @name variogram.intrinsic.spde
#' @title Variogram of intrinsic SPDE model
#' @description Variogram \eqn{\gamma(s_0,s)}{\gamma(s_0,s)} of intrinsic SPDE
#' model
#' \deqn{(-\Delta)^{\beta/2}(\kappa^2-\Delta)^{\alpha/2} (\tau u) = \mathcal{W}}{(-\Delta)^{\beta/2}(\kappa^2-\Delta)^{\alpha/2} (\tau u) = \mathcal{W}}
#' with Neumann boundary conditions and a mean-zero constraint on a
#' square \eqn{[0,L]^d}{[0,L]^d} for \eqn{d=1}{d=1} or \eqn{d=2}{d=2}.
#' @param s0 The location where the variogram should be evaluated, either
#' a double for 1d or a vector for 2d
#' @param s A vector (in 1d) or matrix (in 2d) with all locations where the
#' variogram is computed
#' @param kappa Range parameter.
#' @param alpha Smoothness parameter.
#' @param beta Smoothness parameter.
#' @param tau Precision parameter.
#' @param L The side length of the square domain.
#' @param N The number of terms in the Karhunen-Loeve expansion.
#' @param d The dimension (1 or 2).
#' @param semi Compute the semi variogram? Default FALSE.
#' @details The variogram is computed based on a Karhunen-Loeve expansion of the
#' covariance function.
#'
#' @export
#' @seealso [intrinsic.matern.operators()]
#'
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)) {
#'   x <- seq(from = 0, to = 10, length.out = 201)
#'   beta <- 1
#'   alpha <- 1
#'   kappa <- 1
#'   op <- intrinsic.matern.operators(
#'     kappa = kappa, tau = 1, alpha = alpha,
#'     beta = beta, loc_mesh = x, d = 1
#'   )
#'   # Compute and plot the variogram of the model
#'   Sigma <- op$A[,-1] %*% solve(op$Q[-1,-1], t(op$A[,-1]))
#'   One <- rep(1, times = ncol(Sigma))
#'   D <- diag(Sigma)
#'   Gamma <- 0.5 * (One %*% t(D) + D %*% t(One) - 2 * Sigma)
#'   k <- 100
#'   plot(x, Gamma[k, ], type = "l")
#'   lines(x,
#'     variogram.intrinsic.spde(x[k], x, kappa, alpha, beta, L = 10, d = 1),
#'     col = 2, lty = 2
#'   )
#' }
variogram.intrinsic.spde <- function(s0 = NULL,
                                     s = NULL,
                                     kappa = 0,
                                     alpha = 0,
                                     beta = NULL,
                                     tau = 1,
                                     L = NULL,
                                     N = 100,
                                     d = NULL,
                                     semi = FALSE) {
  if (is.null(kappa) || is.null(alpha) || is.null(beta)) {
    stop("All model parameters must be provided.")
  }
  if (is.null(s0) || is.null(s) || is.null(d) || is.null(L)) {
    stop("s0, s, L and d must be provided.")
  }

  if (d == 1) {
    if (is.matrix(s)) {
      n <- max(dim(s))
      if (min(dim(s)) > 1) {
        stop("s has wrong dimensions for d = 1")
      }
    } else {
      n <- length(s)
    }
    vario <- rep(0, n)
    for (i in 1:N) {
      lambda <- (i * pi / L)^(-2 * beta) * ((i * pi / L)^2 + kappa^2)^(-alpha)
      vario <- vario + 0.5 * (2 / L) * lambda * (cos(i * pi * s / L) - cos(i * pi * s0 / L))^2
    }
  } else if (d == 2) {
    if (!is.matrix(s)) {
      stop("s should be a matrix if d=2")
    }
    vario <- rep(0, dim(s)[1])
    for (i in 1:N) {
      f <- i^2 * pi^2 / L^2
      lambda <- f^(-beta) * (f + kappa^2)^(-alpha)
      e1 <- (sqrt(2) / L) * cos(i * pi * s[, 1] / L)
      e2 <- (sqrt(2) / L) * cos(i * pi * s0[1] / L)
      vario <- vario + 0.5 * lambda * (e1 - e2)^2
    }
    for (i in 1:N) {
      f <- i^2 * pi^2 / L^2
      lambda <- f^(-beta) * (f + kappa^2)^(-alpha)
      e1 <- (sqrt(2) / L) * cos(i * pi * s[, 2] / L)
      e2 <- (sqrt(2) / L) * cos(i * pi * s0[2] / L)
      vario <- vario + 0.5 * lambda * (e1 - e2)^2
    }
    for (i in 1:N) {
      for (j in 1:N) {
        f <- (i^2 + j^2) * pi^2 / L^2
        lambda <- f^(-beta) * (f + kappa^2)^(-alpha)
        e1 <- (2 / L) * cos(i * pi * s[, 1] / L) * cos(j * pi * s[, 2] / L)
        e2 <- (2 / L) * cos(i * pi * s0[1] / L) * cos(j * pi * s0[2] / L)
        vario <- vario + 0.5 * lambda * (e1 - e2)^2
      }
    }
  } else {
    stop("d should be 1 or 2.")
  }
  if(semi) {
      return(vario / tau^2)    
  } else {
      return(2 * vario / tau^2)    
  }
}

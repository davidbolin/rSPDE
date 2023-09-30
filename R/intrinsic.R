#' @name intrinsic.operators
#' @title Covariance-based approximations of intrinsic fields
#' @description `intrinsic.operators` is used for computing a
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
#' "chebfun", "brasil" or "chebfunLB".
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @param scaling second lowest eigenvalue of g1
#' @return `intrinsic.operators` returns an object of
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
intrinsic.operators <- function(C, 
                                G,
                                mesh,
                                alpha,
                                m = 2,
                                d,
                                compute_higher_order = FALSE,
                                return_block_list = FALSE,
                                type_rational_approximation = c("chebfun",
                                                                "brasil", 
                                                                "chebfunLB"),
                                fem_mesh_matrices = NULL,
                                scaling = NULL) {
  type_rational_approximation <- type_rational_approximation[[1]]
  
  if (is.null(fem_mesh_matrices)) {
    if (!is.null(mesh)) {
      d <- get_inla_mesh_dimension(inla_mesh = mesh)
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
      d <- get_inla_mesh_dimension(inla_mesh = mesh)
    }
    Gk <- list()
    Gk[[1]] <- G
    for (i in 2:m_order) {
      Gk[[i]] <- fem_mesh_matrices[[paste0("g", i)]]
    }
  }
  if(is.null(scaling)) {
    scaling <- RSpectra::eigs(as(G,"CsparseMatrix"),2, which = "SM")$values[1]
  }
  
  L <- G / scaling
  
  CiL <- CiG / scaling 
  
  if (m_alpha == 0) {
    aux_mat <- Diagonal(dim(L)[1])
  } else {
    aux_mat <- CiL
  }
  
  
  if (return_block_list) {
    Q.int <- aux_mat
    
    if(alpha %% 1 == 0){
      Q.frac <- Matrix::Diagonal(dim(L)[1])
      Q <- G
      
      if(alpha > 1){
        for(k in 1:(alpha-1)){
          Q <- Q %*% CiG
        }
      } 
      
      Q.int <- Q
    } else {
      if(m == 0){
        stop("Return block list does not work with m = 0, either increase m or set return_block_list to FALSE.")
      }
      Q.frac <- intrinsic.precision(alpha = alpha, rspde.order = m, dim = d,
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
    
    if(alpha %% 1 == 0){
      Q.frac <- Matrix::Diagonal(dim(L)[1])
      Q <- G
      
      if(alpha > 1){
        for(k in 1:(alpha-1)){
          Q <- Q %*% CiL
        }
      } 
      
      Q.int <- list(Q.int = Q, order = m_alpha)
    } else if (m > 0){
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
    } else{
      stop("m > 0 required for intrinsic fields")
    }
  }
  
  ## output
  output <- list(
    C = C, G = G, L = L, Ci = Ci, GCi = GCi, Gk = Gk,
    fem_mesh_matrices = fem_mesh_matrices,
    alpha = alpha, m = m, d = d,
    Q.frac = Q.frac, Q.int = Q.int,
    Q = Q, sizeC = dim(C)[1],
    higher_order = compute_higher_order,
    type_rational_approximation = type_rational_approximation,
    return_block_list = return_block_list,
    stationary = TRUE
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
  
  row_nu <- round(1000*cut_decimals(alpha))
  r <- unlist(mt[row_nu, 2:(1+rspde.order)])
  p <- unlist(mt[row_nu, (2+rspde.order):(1+2*rspde.order)])
  k <- unlist(mt[row_nu, 2+2*rspde.order])
  
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
      Malpha2 <- fem_mesh_matrices[[paste0("g", m_alpha+1)]] / scaling^(m_alpha+1)
      
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
    
    Q <- Q * scaling^alpha
    
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
      if(m_alpha==0) {
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
      if(m_alpha==0) {
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
#' "chebfun", "brasil" or "chebfunLB".
#' @param fem_mesh_matrices A list containing FEM-related matrices.
#' The list should contain elements c0, g1, g2, g3, etc.
#' @param scaling second lowest eigenvalue of g1
#' @return `intrinsic.matern.operators` returns an object of
#' class "intrinsicCBrSPDEobj". This object is a list containing the
#' following quantities:
#' \item{C}{The mass lumped mass matrix.}
#' \item{Ci}{The inverse of `C`.}
#' \item{GCi}{The stiffness matrix G times `Ci`}
#' \item{Gk}{The stiffness matrix G along with the higher-order
#' FEM-related matrices G2, G3, etc.}
#' \item{fem_mesh_matrices}{A list containing the mass lumped mass
#' matrix, the stiffness matrix and
#' the higher-order FEM-related matrices.}
#' \item{m_alpha}{The order of the rational approximation for the Matérn part.}
#' \item{m_beta}{The order of the rational approximation for the intrinsic part.}
#' \item{alpha}{The fractional power of the Matérn part of the operator.}
#' \item{beta}{The fractional power of the intrinsic part of the operator.}
#' \item{type}{String indicating the type of approximation.}
#' \item{d}{The dimension of the domain.}
#' \item{A}{Matrix that sums the components in the approximation to the mesh nodes.}
#' \item{kappa}{Range parameter of the covariance function}
#' \item{tau}{Scale parameter of the covariance function.}
#' \item{type}{String indicating the type of approximation.}
#' @export
#' @details The covariance operator 
#' \deqn{\tau^{-2}(-\Delta)^{\beta}(\kappa^2-\Delta)^{\alpha}}{\tau^{-2}(-\Delta)^{\beta}(\kappa^2-\Delta)^{\alpha}}
#' is approximated based on rational approximations of the two fractional 
#' components. The Laplacians are equipped with homogeneous Neumann boundary
#' conditions and a zero-mean constraint is additionally imposed to obtained 
#' a non-intrinsic model. 
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)){
#'  x <- seq(from = 0, to = 10, length.out = 201)
#'  beta <- 1
#'  alpha <- 1
#'  kappa <- 1
#'  op <- intrinsic.matern.operators(kappa = kappa, tau = 1, alpha = alpha, 
#'                                  beta = beta, loc_mesh = x, d=1) 
#'  # Compute and plot the variogram of the model
#'  Sigma <- op$A %*% solve(op$Q,t(op$A))
#'  One <- rep(1, times = ncol(Sigma))
#'  D <- diag(Sigma)
#'  Gamma <- 0.5*(One %*% t(D) + D %*% t(One) - 2 * Sigma)
#'  k <- 100
#'  plot(x, Gamma[k, ], type = "l")
#'  lines(x, 
#'       variogram.intrinsic.spde(x[k], x, kappa, alpha, beta, L = 10, d = 1),
#'       col=2, lty = 2)
#' }
intrinsic.matern.operators <- function(kappa,
                                       tau,
                                       alpha ,
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
                                       type_rational_approximation = c("chebfun",
                                                                       "brasil", "chebfunLB"),
                                       fem_mesh_matrices = NULL,
                                       scaling = NULL) {
  
  if (is.null(d) && is.null(mesh) && is.null(graph)) {
    stop("You should give either the dimension d, the mesh or graph!")
  }
  
  if ((is.null(C) || is.null(G)) && is.null(mesh) && is.null(graph) &&(is.null(loc_mesh) || d != 1)) {
    stop("You should either provide mesh, graph, or provide both C *and* G!")
  }
  
  if( (is.null(C) || is.null(G)) && (is.null(graph)) && (!is.null(loc_mesh) && d==1)){
    fem <- rSPDE.fem1d(loc_mesh)
    C <- fem$C
    G <- fem$G

    fem_mesh_matrices <- list()
    fem_mesh_matrices[["c0"]] <- C
    fem_mesh_matrices[["g1"]] <- G

    Gk <- list()
    Gk[[1]] <- G

    m_alpha <- floor(alpha)
    m_order <- m_alpha + 1

    if (compute_higher_order) {
      for (i in 1:m_order) {
        fem_mesh_matrices[[paste0("g", i)]] <- Gk[[i]]
      }
    }    
  }
  
  has_mesh <- FALSE
  has_graph <- FALSE
  
  if(!is.null(loc_mesh)){
    if(!is.numeric(loc_mesh)){
      stop("loc_mesh must be numerical.")
    }
  }
  
  if(!is.null(graph)){
    if(!inherits(graph, "metric_graph")){
      stop("graph should be a metric_graph object!")
    }
    d <- 1
    if(is.null(graph$mesh)){
      warning("The graph object did not contain a mesh, one was created with h = 0.01. Use the build_mesh() method to replace this mesh with a new one.")
      graph$build_mesh(h = 0.01)
    }
    graph$compute_fem()
    C <- graph$mesh$C
    G <- graph$mesh$G
    has_graph <- TRUE
  }
  
  if (!is.null(mesh)) {
    d <- get_inla_mesh_dimension(inla_mesh = mesh)
    # fem <- INLA::inla.mesh.fem(mesh)
    fem <- fmesher::fm_fem(mesh)
    C <- fem$c0
    G <- fem$g1
    has_mesh <- TRUE
  }
  
  
  
  kappa <- rspde_check_user_input(kappa, "kappa" , 0)
  tau <- rspde_check_user_input(tau, "tau" , 0)
  
  alpha <- rspde_check_user_input(alpha, "alpha" , 0)
  alpha <- min(alpha, 10)
  
  beta <- rspde_check_user_input(beta, "beta" , 0)
  beta <- min(beta, 10)
  
  if(alpha + beta < d/2) {
    stop("One must have alpha + beta > d/2")
  }
  
  if(!is.null(mesh)){
    make_A <- function(loc){
      # return(INLA::inla.spde.make.A(mesh = mesh, loc = loc))
      return(fmesher::fm_basis(x = mesh, loc = loc))
    }
  } else if(!is.null(graph)){
    make_A <- function(loc){
      return(graph$mesh_A(loc))
    }
  } else if(!is.null(loc_mesh) && d == 1){
    make_A <- function(loc){
      return(rSPDE::rSPDE.A1d(x = loc_mesh, loc = loc))
    }   
  } else {
    make_A <- NULL
  }
  
  if(alpha>0 && beta>0 && kappa > 0) {
    op1 <- CBrSPDE.matern.operators(
      C = C, G = G, mesh = mesh, nu = alpha - d/2, kappa = kappa, tau = tau,
      fem_mesh_matrices = fem_mesh_matrices,
      m = m_alpha, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE,
      type_rational_approximation = type_rational_approximation[[1]]
    )
    op2 <-intrinsic.operators(
      C = C, G = G, mesh = mesh, alpha = beta, 
      m = m_beta, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE, fem_mesh_matrices = fem_mesh_matrices,
      type_rational_approximation = type_rational_approximation[[1]],
      scaling = scaling
    )
    block_list <- list()
    if(is.list(op1$Q)) {
      Q.list1 <- op1$Q  
    } else {
      Q.list1 <- list(op1$Q)
    }
    if(is.list(op2$Q)) {
      Q.list2 <- op2$Q  
    } else {
      Q.list2 <- list(op2$Q)
    }
    m1 <- length(Q.list1)
    m2 <- length(Q.list2)
    k <- 1
    if(return_block_list) {
      Q <- list()
    } 
    for(i in 1:m1) {
      for(j in 1:m2) {
        if(return_block_list) {
          Q[[k]] <- Q.list1[[i]]%*%op1$Ci%*%Q.list2[[j]]  
        } else {
          if(i == 1 && j == 1) {
            Q <- Q.list1[[i]]%*%op1$Ci%*%Q.list2[[j]]  
          } else {
            Q <- bdiag(Q, Q.list1[[i]]%*%op1$Ci%*%Q.list2[[j]])
          }  
        }
      }
    }
    h <- rep(rowSums(op1$C),m1*m2)
    if(!return_block_list) {
      Q <- rbind(cbind(Q, h), c(h,0))  
    }  
    n <- dim(op1$C)[1]
    A <- cbind(kronecker(matrix(rep(1,m1*m2) , 1, m1*m2), Diagonal(n)), 
               Matrix(0,ncol=1,nrow=n))
  } else if (alpha > 0 && kappa > 0) {
    op1 <- CBrSPDE.matern.operators(
      C = C, G = G, mesh = mesh, nu = alpha - d/2, kappa = kappa, tau = tau,
      fem_mesh_matrices = fem_mesh_matrices,
      m = m_alpha, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE,
      type_rational_approximation = type_rational_approximation[[1]]
    )
    if(is.list(op1$Q)) {
      Q.list1 <- op1$Q  
    } else {
      Q.list1 <- list(op1$Q)
    }
    m1 <- length(Q.list1)
    if(!return_block_list) {
      for(i in 1:m1) {
        if(i == 1) {
          Q <- Q.list1[[i]]
        } else {
          Q <- bdiag(Q, Q.list1[[i]])
        }  
      }  
    } else {
      Q <- Q.list1 
    }
    h <- rep(rowSums(op1$C),m1)
    if(!return_block_list) {
      Q <- rbind(cbind(Q, h), c(h,0))  
    } 
    n <- dim(op1$C)[1]
    A <- cbind(kronecker(matrix(rep(1,m1), 1, m1), Diagonal(n)), 
               Matrix(0,ncol=1,nrow=n))
  } else if (beta > 0) {
    if(kappa == 0) {
      alpha_beta = alpha + beta
    } else {
      alpha_beta = beta
    }
    op1 <-intrinsic.operators(
      C = C, G = G, mesh = mesh, alpha = alpha_beta, 
      m = m_beta, d = d, compute_higher_order = compute_higher_order,
      return_block_list = TRUE, fem_mesh_matrices = fem_mesh_matrices,
      type_rational_approximation = type_rational_approximation[[1]]
    )
    if(is.list(op1$Q)) {
      Q.list1 <- op1$Q  
    } else {
      Q.list1 <- list(op1$Q)
    }
    m1 <- length(Q.list1)
    if(!return_block_list) {
      for(i in 1:m1) {
        if(i == 1) {
          Q <- Q.list1[[i]]
        } else {
          Q <- bdiag(Q, Q.list1[[i]])
        }  
      }  
    } else {
      Q <- Q.list1 
    }
    h <- rep(rowSums(op1$C),m1)
    if(!return_block_list) {
      Q <- rbind(cbind(Q, h), c(h,0))  
    }
    n <- dim(op1$C)[1]
    A <- cbind(kronecker(matrix(rep(1,m1), 1, m1), Diagonal(n)), 
               Matrix(0,ncol=1,nrow=n))
  } 
  
  out <- list(C = op1$C, G = op1$G, Ci = op1$Ci, GCi = op1$GCi,
              Q = Q, 
              alpha = alpha, beta = beta, kappa = kappa, tau = tau,
              m_alpha = m_alpha, m_beta = m_beta, d = d,
              type_rational_approximation = type_rational_approximation[[1]],
              higher_order = compute_higher_order,
              return_block_list = return_block_list,
              stationary = TRUE, 
              has_mesh = has_mesh,
              has_graph = has_graph,
              make_A = make_A,
              A = A,
              mesh = mesh,
              graph = graph
  )
  out$type <- "Covariance-Based intrinsic Matern SPDE Approximation"
  class(out) <- "intrinsicCBrSPDEobj"
  
  return(out)
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
#' @details The variogram is computed based on a Karhunen-Loeve expansion of the 
#' covariance function. 
#' 
#' @export
#' @seealso [intrinsic.matern.operators()]
#'
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE)){
#'  x <- seq(from = 0, to = 10, length.out = 201)
#'  beta <- 1
#'  alpha <- 1
#'  kappa <- 1
#'  op <- intrinsic.matern.operators(kappa = kappa, tau = 1, alpha = alpha, 
#'                                  beta = beta, loc_mesh = x, d=1) 
#'  # Compute and plot the variogram of the model
#'  Sigma <- op$A %*% solve(op$Q,t(op$A))
#'  One <- rep(1, times = ncol(Sigma))
#'  D <- diag(Sigma)
#'  Gamma <- 0.5*(One %*% t(D) + D %*% t(One) - 2 * Sigma)
#'  k <- 100
#'  plot(x, Gamma[k, ], type = "l")
#'  lines(x, 
#'       variogram.intrinsic.spde(x[k], x, kappa, alpha, beta, L = 10, d = 1),
#'       col=2, lty = 2)
#' }
variogram.intrinsic.spde <- function(s0 = NULL, 
                                     s = NULL, 
                                     kappa = NULL, 
                                     alpha = NULL, 
                                     beta = NULL,
                                     tau = 1,
                                     L = NULL,
                                     N = 100,
                                     d = NULL) {
  if(is.null(kappa) || is.null(alpha) || is.null(beta)) {
    stop("All model parameters must be provided.")
  }
  if(is.null(s0) || is.null(s) || is.null(d) || is.null(L)) {
    stop("s0, s, L and d must be provided.")
  }
  
  if(d==1) {
    if(is.matrix(s)) {
      n = max(dim(s))
      if(min(dim(s))>1) {
        stop("s has wrong dimensions for d = 1")
      }
    } else {
      n = length(s)
    }
    vario <- rep(0,n)
    for(i in 1:N) {
      lambda <- (i*pi/L)^(-2*beta)*((i*pi/L)^2+kappa^2)^(-alpha)
      vario <- vario + 0.5*(2/L)*lambda*(cos(i*pi*s/L)-cos(i*pi*s0/L))^2
    }
  } else if (d == 2) {
    if(!is.matrix(s)) {
      stop("s should be a matrix if d=2")
    }
    vario <- rep(0,dim(s)[1])
    for(i in 1:N) {
      f <- i^2*pi^2/L^2
      lambda <- f^(-beta)*(f+kappa^2)^(-alpha)    
      e1 <- (sqrt(2)/L)*cos(i*pi*s[,1]/L)
      e2 <- (sqrt(2)/L)*cos(i*pi*s0[1]/L)
      vario <- vario + 0.5*lambda*(e1-e2)^2
    }
    for(i in 1:N) {
      f <- i^2*pi^2/L^2
      lambda <- f^(-beta)*(f+kappa^2)^(-alpha)    
      e1 <- (sqrt(2)/L)*cos(i*pi*s[,2]/L)
      e2 <- (sqrt(2)/L)*cos(i*pi*s0[2]/L)
      vario <- vario + 0.5*lambda*(e1-e2)^2
    }
    for(i in 1:N) {
      for(j in 1:N) {
        f <- (i^2+j^2)*pi^2/L^2
        lambda <- f^(-beta)*(f+kappa^2)^(-alpha)    
        e1 <- (2/L)*cos(i*pi*s[,1]/L)*cos(j*pi*s[,2]/L)
        e2 <- (2/L)*cos(i*pi*s0[1]/L)*cos(j*pi*s0[2]/L)
        vario <- vario + 0.5*lambda*(e1-e2)^2
      }
    }  
  } else {
    stop("d should be 1 or 2.")
  }
  return(vario/tau^2)
}


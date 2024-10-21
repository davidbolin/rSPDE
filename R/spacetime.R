
#' Space-time random fields 
#' 
#' `spacetime.operators` is used for computing a FEM approximation of a Gaussian 
#' random field defined as a solution to the SPDE
#' \deqn{d u + \gamma(\kappa^2 + \rho\cdot\Nabla - \Delta)^\alpha u = \sigma dW_C.}
#'where C is a Whittle-Matern covariance operator with smoothness parameter 
#'\eqn{\beta} and range parameter \eqn{\kappa}
#' @param mesh Spatial mesh for FEM approximation
#' @param time Time points for the temporal discretization
#' @param graph An optional `metric_graph` object. Replaces `mesh` for models on 
#' metric graphs.
#' @param kappa Positive spatial range parameter
#' @param sigma Positive variance parameter
#' @param rho Drift parameter. Real number on metric graphs and 
#' one-dimensional spatial domains, a vector with two number on 2d domains.
#' @param gamma Temporal range parameter.  
#' @param alpha Integer smoothness parameter alpha.
#' @param beta Integer smoothness parameter beta.
#'
#' @return An object of type spacetimeobj. 
#' @export
#'
#' @examples
#'s <- seq(from = 0, to = 20, length.out = 101)
#'mesh <- fm_mesh_1d(s)
#'
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(mesh = mesh, time = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#'Q <- op_cov$Q
#'v <- rep(0,dim(Q)[1])
#'v[1565] <- 1
#'Sigma <- solve(Q,v)
#'
#'image(matrix(Sigma, nrow=length(s), ncol = length(t)))
spacetime.operators <- function(mesh = NULL,
                                time = NULL,
                                graph = NULL,
                                kappa = NULL,
                                sigma = NULL,
                                gamma = NULL,
                                rho = NULL,
                                alpha = NULL,
                                beta = NULL) {
    
    if (is.null(mesh) && is.null(graph)) {
        stop("You should provide mesh_space or graph.")
    }
    
    if (!is.null(mesh) && !is.null(graph)) {
        stop("You should provide either mesh_space or graph.")
    }
    
    if (is.null(time)) {
        stop("You should provide time.")
    }
    
    if (is.null(alpha)) {
        stop("alpha must be provided")
    } else if(alpha%%1 != 0){
        stop("alpha must be integer")
    } 
    
    if (is.null(beta)) {
        stop("beta must be provided")
    } else if(beta%%1 != 0){
        stop("beta must be integer")
    } 
    
    
    nt <- length(time)
    d <- c(Inf, diff(time))
    dm1 <- c(d[2:nt], Inf)
    Gt <- -bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                      diagonals = cbind(1 / dm1, -(1 / dm1 + 1 / d), 1 / dm1))
    Ct <- bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                     diagonals = cbind(dm1 / 6, (dm1 + d) / 3, c(d[2:nt],Inf) / 6))
    Ct[1, 1:2] <- c(d[2], d[2] / 2) / 3
    Ct[nt, (nt - 1):nt] <- c(d[nt] / 2, d[nt]) / 3
    Bt <- bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                     diagonals = cbind(rep(0.5,nt), rep(0,nt), rep(-0.5,nt)))
    Bt[1,1] = -0.5
    Bt[nt,nt] = 0.5
    B0 <- Diagonal(n=nt,0)
    B0[1,1] <- B0[nt,nt] <- 1/2
    
    has_mesh <- FALSE
    has_graph <- FALSE
    
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
        Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
        G <- graph$mesh$G
        B <- graph$mesh$B
        has_graph <- TRUE
    }
    
    if (!is.null(mesh)) {
        d <- get_inla_mesh_dimension(inla_mesh = mesh)
        if(d==2){
            P <- mesh$loc[,1:2]
            FV <- mesh$graph$tv
            fem <- rSPDE.fem2d(FV, P)
            C <- fem$Cd
            Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
            G <- fem$G
            Bx <- fem$Bx
            By <- fem$By
        } else if(d==1){
            fem <- rSPDE.fem1d(mesh$loc)
            C <- fem$Cd
            Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
            G <- fem$G
            B <- fem$B
        } else {
            stop("Only supported for 1d and 2d meshes.")
        }
        
        has_mesh <- TRUE
    } 
    if(alpha+beta<d/2){
        stop("You must have alpha+beta >= d, where d is the dimension of the spatial domain")
    }
    
    if (is.null(kappa) || is.null(sigma)) {
        if (has_mesh) {
            param <- get.initial.values.rSPDE(dim = d, parameterization = "SPDE", 
                                              mesh = mesh, nu = alpha + beta - d/2)
        } else if (has_graph) {
            param <- get.initial.values.rSPDE(graph.obj = graph, parameterization = "SPDE", 
                                              nu = alpha + beta - d/2)
        }
    }
        
    if (is.null(kappa)) {
        kappa <- exp(param[2])
    } else {
        kappa <- rspde_check_user_input(kappa, "kappa", 0, 1)
    }
    if (is.null(sigma)) {
        sigma <- exp(-param[1])
    } else {
        sigma <- rspde_check_user_input(sigma, "sigma", 0, 1)
    }
    
    if(is.null(rho)){
        rho <- 0
    } else {
        if(has_graph) {
            rho <- rspde_check_user_input(rho, "rho", dim = 1) 
        } else {
            if(d==1){
                rho <- rspde_check_user_input(rho, "rho", dim = 1)     
            } else {
                rho <- rspde_check_user_input(rho, "rho", dim = 2) 
            }
            
        }
    }
    
    if (is.null(gamma)) {
        gamma <- 2/range(time)
    } else {
        gamma <- rspde_check_user_input(gamma, "gamma", 0, 1)
    }
    
    Glist <- make.Glist(beta+2*alpha, C, G)
    Ctlist <- kron.Glist(Ct,Glist, left = TRUE)
    Gtlist <- kron.Glist(Gt,Glist, left = TRUE)
    B0list <- kron.Glist(B0,Glist, left = TRUE)
    M1list <- list()
    M2list <- list()
    M2list2 <- list()
    for(k in 0:alpha) {
        M1list.tmp <- mult.Glist(Ci%*%Glist[[k+1]], Glist, left = FALSE)
        M1list[[k+1]] <- kron.Glist(Ct, M1list.tmp)
        
        if(d==2){
            M2list.tmp <- mult.Glist(Ci%*%Glist[[floor(k/2)+1]]%*%Ci%*%Bx, Glist, left = FALSE)
            M2list[[k+1]] <- kron.Glist(t(Bt), M2list.tmp)
            M2list.tmp <- mult.Glist(Ci%*%Glist[[floor(k/2)+1]]%*%Ci%*%By, Glist, left = FALSE)
            M2list2[[k+1]] <- kron.Glist(t(Bt), M2list.tmp)
        } else {
            M2list.tmp <- mult.Glist(Ci%*%Glist[[floor(k/2)+1]]%*%Ci%*%B, Glist, left = FALSE)
            M2list[[k+1]] <- kron.Glist(t(Bt), M2list.tmp)
        }
    }
    if(alpha == 0) {
        Q <- make.L(beta,kappa,Gtlist) + gamma^2*make.L(beta,kappa, Ctlist)
    } else if (alpha==1) {
        Q <- make.L(beta,kappa,Gtlist) + 2*gamma*make.L(beta+alpha,kappa, B0list)
        Q <- Q + gamma^2*make.L(beta+2*alpha, kappa, Ctlist) - gamma^2*rho^2*Ctlist[[3]]
        Q <- Q + 6*gamma^2*rho^2*make.L(beta+2,kappa,M1list[[2]])
        M2 <- make.L(beta+1,kappa,M2list[[1]])
        Q <- Q - 2*rho*gamma*(M2 + t(M2))
    } else {
        Q <- make.L(beta,kappa,Gtlist) + 2*gamma*make.L(beta+alpha,kappa, B0list)
        
        for(k in 0:alpha) {
            Q <- Q + gamma^2*choose(alpha,k)*rho^(2*k)*make.L(beta+2*(alpha-k),kappa, M1list[[k+1]])
            if(d==2){
                M2x <- make.L(beta+alpha-k,kappa,M2list[[k+1]])
                M2y <- make.L(beta+alpha-k,kappa,M2list2[[k+1]])
                Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho[1]^(k)*(M2x + t(M2x))
                Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho[2]^(k)*(M2y + t(M2y))
            } else {
                M2 <- make.L(beta+alpha-k,kappa,M2list[[k+1]])
                Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho^(k)*(M2 + t(M2))    
            }
        }
    }
    
    out <- list()
    out$Q <- Q/sigma^2
    out$Gtlist <- Gtlist
    out$Ctlist <- Ctlist
    out$B0list <- B0list
    out$M1list <- M1list
    out$M2list <- M2list
    out$M2list2 <- M2list2
    out$kappa <- kappa
    out$sigma <- sigma
    out$alpha <- alpha
    out$beta <- beta
    out$gamma <- gamma
    out$rho <- rho
    out$has_mesh <- has_mesh
    out$has_graph <- has_graph
    out$time <- time
    out$mesh <- mesh
    out$graph <- graph
    class(out) <- "spacetimeobj"
    return(out)

}



#' Update space-time operator object with new parameters
#'
#' @param object Space-time object created by [spacetime.operators()]
#' @param user_kappa kappa value to be updated.
#' @param user_sigma sigma value to be updated.
#' @param user_gamma gamma value to be updated.
#' @param user_rho rho value to be updated.
#'
#' @return An object of type spacetimeobj with updated parameters.
#' @export
#'
#' @examples
#'s <- seq(from = 0, to = 20, length.out = 101)
#'mesh <- fm_mesh_1d(s)
#'
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(mesh = mesh, time = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#'op_cov <- update.spacetimeobj(op_cov, user_kappa = 4, 
#'                              user_sigma = 2, user_gamma = 0.1)                              
update.spacetimeobj <- function(object, 
                                user_kappa = NULL,
                                user_sigma = NULL,
                                user_gamma = NULL,
                                user_rho = NULL) {
    new_object <- object
    
    if (!is.null(user_kappa)) {
        kappa <- rspde_check_user_input(user_kappa, "kappa", 0, 1)
    } else {
        kappa <- object$kappa
    }
    
    if (!is.null(user_sigma)) {
        sigma <- rspde_check_user_input(user_sigma, "sigma", 0, 1)
    } else {
        sigma <- object$sigma
    }
    
    if(!is.null(user_rho)){
        if(object$has_graph) {
            rho <- rspde_check_user_input(user_rho, "rho", dim = 1) 
        } else {
            d <- get_inla_mesh_dimension(inla_mesh = mesh)
            rho <- rspde_check_user_input(user_rho, "rho", dim = d) 
        }
    } else {
        rho <- object$rho
    }
    
    if (!is.null(user_gamma)) {
        gamma <- rspde_check_user_input(user_gamma, "gamma", 0, 1)
    } else {
        gamma <- object$gamma
    }
    
    alpha <- object$alpha
    beta <- object$beta 
    
    if(alpha == 0) {
        Q <- make.L(beta,kappa,object$Gtlist) + gamma^2*make.L(beta,kappa, object$Ctlist)
    } else if (alpha==1) {
        Q <- make.L(beta,kappa,object$Gtlist) + 2*gamma*make.L(beta+alpha,kappa, object$B0list)
        Q <- Q + gamma^2*make.L(beta+2*alpha, kappa, object$Ctlist) - gamma^2*rho^2*object$Ctlist[[3]]
        Q <- Q + 6*gamma^2*rho^2*make.L(beta+2,kappa,object$M1list[[2]])
        M2 <- make.L(beta+1,kappa,object$M2list[[1]])
        Q <- Q - 2*rho*gamma*(M2 + t(M2))
    } else {
        Q <- make.L(beta,kappa,object$Gtlist) + 2*gamma*make.L(beta+alpha,kappa, object$B0list)
        
        for(k in 0:alpha) {
            Q <- Q + gamma^2*choose(alpha,k)*rho^(2*k)*make.L(beta+2*(alpha-k),kappa, object$M1list[[k+1]])
            if(d==2){
                M2x <- make.L(beta+alpha-k,kappa,object$M2list[[k+1]])
                M2y <- make.L(beta+alpha-k,kappa,object$M2list2[[k+1]])
                Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho[1]^(k)*(M2x + t(M2x))
                Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho[2]^(k)*(M2y + t(M2y))
            } else {
                M2 <- make.L(beta+alpha-k,kappa,object$M2list[[k+1]])
                Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho^(k)*(M2 + t(M2))    
            }
        }
    }
    
    new_object$Q <- Q/sigma^2
    new_object$kappa <- kappa
    new_object$sigma <- sigma
    new_object$gamma <- gamma
    new_object$rho <- rho
    return(new_object)
}

#' Simulation of space-time models
#'
#' @param object Space-time object created by [spacetime.operators()]
#' @param nsim The number of simulations.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param user_kappa kappa parameter if it should be updated
#' @param user_sigma sigma parameter if it should be updated
#' @param user_gamma gamma parameter if it should be updated
#' @param user_rho rho parameter if it should be updated
#' @param ... Currently not used.
#'
#' @return A matrix with the simulations as columns.
#' @export
#'
#' @examples
#'s <- seq(from = 0, to = 20, length.out = 101)
#'mesh <- fm_mesh_1d(s)
#'
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(mesh = mesh, time = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#'x <- simulate(op_cov, nsim = 1) 
#'image(matrix(x, nrow = length(s), ncol = length(t)))                              
simulate.spacetimeobj <- function(object, nsim = 1,
                                seed = NULL,
                                user_kappa = NULL,
                                user_sigma = NULL,
                                user_gamma = NULL,
                                user_rho = NULL,
                                ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    object <- update.spacetimeobj(object, user_kappa = user_kappa,
                                  user_sigma = user_sigma,
                                  user_gamma = user_gamma,
                                  user_rho = user_rho)
        
    sizeQ <- dim(object$Q)[1]
    Z <- rnorm(sizeQ * nsim)
    dim(Z) <- c(sizeQ, nsim)
        
    LQ <- chol(forceSymmetric(object$Q))
    X <- solve(LQ, Z)
    
    return(X)
}

## Util functions below 

# Build the operator (kappa^2 C + G)^n from a Glist
#' @noRd
make.L <- function(n,kappa,Glist){
    if(n > length(Glist)){
        stop("Glist too short.")
    }
    L <- 0
    for(k in 0:n){
        L <- L + choose(n,k)*kappa^(2*(n-k))*Glist[[k+1]]
    }
    return(L)
}

#make a list with elements G[[1]] = C, G[[2]] = G, G[[k]] = G%*%solve(C,G[[k-1]])
#' @noRd
make.Glist <- function(n,C, G){
    Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
    Glist <- list() 
    
    for(k in 0:n){
        if(k==0){
            Gk <- C        
        } else if(k==1) {
            Gk <- G
        } else {
            Gk <- Gk%*%Ci%*%G
        }
        Glist[[k+1]] <- Gk
    }
    return(Glist)
}

# make a list with elements M%*%Glist[[k]] (if left = TRUE) or Glist[[k]]%*%M
#' @noRd
mult.Glist <- function(M,Glist, left = TRUE){
    Mlist <- list()
    for(k in 1:length(Glist)){
        if(left) {
            Mlist[[k]] <- M%*%Glist[[k]]    
        } else {
            Mlist[[k]] <- Glist[[k]]%*%M
        }
    }
    return(Mlist)
}

# make a list with elements kronecker(M,Glist[[k]]) (if left = TRUE) or kronecker(Glist[[k]],M)
#' @noRd
kron.Glist <- function(M,Glist, left = TRUE){
    Mlist <- list()
    for(k in 1:length(Glist)){
        if(left) {
            Mlist[[k]] <- kronecker(M,Glist[[k]])
        } else {
            Mlist[[k]] <- kronecker(Glist[[k]],M)
        }
    }
    return(Mlist)
}


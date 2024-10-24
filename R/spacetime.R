
#' Space-time random fields 
#' 
#' `spacetime.operators` is used for computing a FEM approximation of a Gaussian 
#' random field defined as a solution to the SPDE
#' \deqn{d u + \gamma(\kappa^2 + \rho\cdot\Nabla - \Delta)^\alpha u = \sigma dW_C.}
#'where C is a Whittle-Matern covariance operator with smoothness parameter 
#'\eqn{\beta} and range parameter \eqn{\kappa}
#' @param mesh_space Spatial mesh for FEM approximation
#' @param mesh_time Temporal mesh for FEM approximation
#' @param space_loc Locations of mesh nodes for spatial mesh for 1d models.
#' @param time_loc Locations of temporal mesh nodes.
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
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(space_loc = s, time_loc = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#'Q <- op_cov$Q
#'v <- rep(0,dim(Q)[1])
#'v[1565] <- 1
#'Sigma <- solve(Q,v)
#'
#'image(matrix(Sigma, nrow=length(s), ncol = length(t)))
spacetime.operators <- function(mesh_space = NULL,
                                mesh_time = NULL,
                                space_loc = NULL,
                                time_loc = NULL,
                                graph = NULL,
                                kappa = NULL,
                                sigma = NULL,
                                gamma = NULL,
                                rho = NULL,
                                alpha = NULL,
                                beta = NULL) {
    
    
    if ((!is.null(mesh_space) && !is.null(graph)) || (!is.null(mesh_space) && !is.null(space_loc)) || (!is.null(graph) && !is.null(space_loc))){
        stop("You should provide only one of mesh_space, space_loc or graph.")
    }
    
    if (!is.null(mesh_space) && !is.null(graph) && !is.null(space_loc)) {
        stop("You should provide mesh_space, space_loc or graph.")
    }
    
    if (is.null(mesh_time) && is.null(time_loc)) {
        stop("You should provide mesh_time or time_loc.")
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
    
    if(!is.null(mesh_time)){
        d_time <- get_inla_mesh_dimension(inla_mesh = mesh_time)
        if(d_time != 1) {
            stop("mesh_time should be a 1d mesh")
        }
        time <- mesh_time$loc
    } else {
        time <- time_loc
        mesh_time <- mesh <- fm_mesh_1d(time)
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
    } else if (!is.null(mesh_space)) {
        mesh <- mesh_space
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
    } else {
        space_loc <- as.matrix(space_loc)
        if(min(dim(space_loc))>1) {
            stop("For 2d domains, please provide mesh_space instead of space_loc")
        }
        fem <- rSPDE.fem1d(space_loc)
        C <- fem$Cd
        Ci <- Diagonal(1/rowSums(C),n=dim(C)[1])
        G <- fem$G
        B <- fem$B
        d <- 1
        mesh_space <- fm_mesh_1d(space_loc)
    }
    if(alpha+beta<d/2){
        stop("You must have alpha+beta >= d, where d is the dimension of the spatial domain")
    }
    
    if (is.null(kappa) || is.null(sigma)) {
        if (has_mesh) {
            param <- get.initial.values.rSPDE(dim = d, parameterization = "spde", 
                                              mesh = mesh_space, nu = alpha + beta - d/2)
        } else if (has_graph) {
            param <- get.initial.values.rSPDE(graph.obj = graph, parameterization = "spde", 
                                              nu = alpha + beta - d/2)
        } else {
            param <- get.initial.values.rSPDE(dim = 1, parameterization = "spde", 
                                              mesh.range = diff(range(space_loc)), 
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
        gamma <- 2/diff(range(time))
    } else {
        gamma <- rspde_check_user_input(gamma, "gamma", 0, 1)
    }
    
    Glist <- make.Glist(beta+3*alpha, C, G)
    
    Ctlist <- kron.Glist(Ct,make.Glist(beta+3*alpha, C, G), left = TRUE)
    Gtlist <- kron.Glist(Gt,make.Glist(1+beta, C, G), left = TRUE)
    B0list <- kron.Glist(B0,make.Glist(beta+alpha, C, G), left = TRUE)
    M2list <- list()
    M2list2 <- list()
    for(k in 0:alpha) {
        Glist <- make.Glist(1+beta+alpha-k, C, G)
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
    
    Q <- make.L(beta,kappa,Gtlist) + 2*gamma*make.L(beta+alpha,kappa, B0list)
        
    for(k in 0:alpha) {
        Q <- Q + gamma^2*choose(alpha,k)*rho^(2*k)*make.L(beta+2*(alpha-k),kappa,Ctlist[(k+1):length(Ctlist)])

        if(d==2){
            M2x <- make.L(beta+alpha-k,kappa,M2list[[k+1]])
            M2y <- make.L(beta+alpha-k,kappa,M2list2[[k+1]])
            Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho[1]^(k)*(M2x + t(M2x))
            Q <- Q - 0.5*gamma*choose(alpha,k)*(1-(-1)^k)*rho[2]^(k)*(M2y + t(M2y))
        } else {
            M2 <- make.L(beta+alpha-k,kappa,M2list[[k+1]])
            Q <- Q - 0.5*(-1)^(floor(k/2))*gamma*choose(alpha,k)*(1-(-1)^k)*rho^(k)*(M2 + t(M2))    
        }
    }

    if (!is.null(graph)) {
        make_A <- function(loc, time) {
            return(rSPDE.Ast(graph = graph, mesh_time = mesh_time, 
                             obs.s = loc, obs.t = time))
        }
    } else {
        make_A <- function(loc, time) {
            return(rSPDE.Ast(mesh_space = mesh_space, mesh_time = mesh_time, 
                             obs.s = loc, obs.t = time))
        }
    }
    
    if(has_graph) { 
        plot_covariances <- function(t.ind, s.ind, t.shift=0, show.temporal = TRUE) {
            
            N <- dim(Q)[1]
        
            n <- N/length(mesh_time$loc)
            
            T <- N/n
            if(length(t.ind)>4)
                stop("max 4 curves allowed")
            if(s.ind > n)
                stop("too large space index")
            if(max(t.ind)>length(mesh_time$loc))
                stop("too large time index")
            
            cols <- c("green", "cyan","blue","red")[(4-length(t.ind)+1):4]
            
            time.index <- n*(0:(T-1)) + s.ind
            ct <- matrix(0,nrow = length(t.ind),ncol = T)
            for(i in 1:length(t.ind)) {
                v <- rep(0,N)
                v[(t.ind[i]-1)*n+s.ind] <- 1
                
                tmp <- solve(Q,v)
        
                ct[i,] <- tmp[time.index]
                for(j in 1:length(t.shift)) {
                    ind <- ((t.ind[i]-t.shift[j]-1)*n+1):((t.ind[i]-t.shift[j])*n)
                    c <- tmp[ind]
                    if(length(t.shift)>1) {
                        col <- cols[j]
                    } else {
                        col <- cols[i]
                    }
                    if(i == 1) {
                        p <- graph$plot_function(as.vector(c), 
                                                 plotly = TRUE, 
                                                 support_width = 0, 
                                                 line_color = col)
                    } else {
                        p <- graph$plot_function(as.vector(c), 
                                                 plotly = TRUE, 
                                                 p = p, 
                                                 support_width = 0, 
                                                 line_color = col)
                    }
                    
                }
            }
            if(show.temporal){
                df <- data.frame(t=rep(mesh_time$loc,length(t.ind)),
                                 y=c(t(ct)), 
                                 i=rep(1:length(t.ind), each=length(mesh_time$loc)))
                pt <- plotly::plot_ly(df, x = ~t, y = ~y, split = ~i, 
                                      type = 'scatter', mode = 'lines')
                fig <- plotly::layout(plotly::subplot(p,pt), 
                                      title = "Marginal covariances",
                                      scene = list(domain=list(x=c(0,0.5),
                                                               y=c(0,1))),
                                      scene2 = list(domain=list(x=c(0.5,1),
                                                                y=c(0,1))))
                fig$x$layout <- fig$x$layout[grep('NA', names(fig$x$layout), 
                                                  invert = TRUE)]
            } else {
                fig <- p
            }
            
            print(fig)
            return(fig)
        }
    } else {
        plot_covariances <- NULL
    }
    out <- list()
    out$Q <- Q/sigma^2
    out$Gtlist <- Gtlist
    out$Ctlist <- Ctlist
    out$B0list <- B0list
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
    out$mesh_time <- mesh_time
    out$mesh_space <- mesh_space
    out$graph <- graph
    out$d <- d
    out$make_A <- make_A
    out$plot_covariances <- plot_covariances
    out$stationary <- TRUE
    
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
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(space_loc = s, time_loc = t,
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
            rho <- rspde_check_user_input(user_rho, "rho", dim = object$d) 
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
        Q <- Q + 6*gamma^2*rho^2*make.L(beta+2,kappa,object$Ctlist[2:length(object$Ctlist)])
        M2 <- make.L(beta+1,kappa,object$M2list[[1]])
        Q <- Q - 2*rho*gamma*(M2 + t(M2))
    } else {
        Q <- make.L(beta,kappa,object$Gtlist) + 2*gamma*make.L(beta+alpha,kappa, object$B0list)
        
        for(k in 0:alpha) {
            Q <- Q + gamma^2*choose(alpha,k)*rho^(2*k)*make.L(beta+2*(alpha-k),kappa,
                                                              object$Ctlist[(k+1):length(object$Ctlist)])
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
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(space_loc = s, time_loc = t,
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
    
    return(as.matrix(X))
}





#' @name precision.spacetimeobj
#' @title Get the precision matrix of spacetimeobj objects
#' @description Function to get the precision matrix of a spacetimeobj object
#' @param object The model object computed using [spacetime.operators()]
#' @param user_kappa If non-null, update the range parameter of
#' the covariance function.
#' @param user_sigma If non-null, update the standard deviation of
#' the covariance function.
#' @param user_gamma If non-null, update the temporal range parameter
#' of the covariance function.
#' @param user_rho If non-null, update the drift parameter of the
#' covariance function.
#' @param ... Currently not used.
#' @return The precision matrix.
#' @method precision spacetimeobj
#' @seealso [simulate.spacetimeobj()], [spacetime.operators()]
#' @export
#' @examples
#'s <- seq(from = 0, to = 20, length.out = 101)
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(space_loc = s, time_loc = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#' prec_matrix <- precision(op_cov)
precision.spacetimeobj <- function(object,
                                   user_kappa = NULL,
                                   user_sigma = NULL,
                                   user_gamma = NULL,
                                   user_rho = NULL,
                                   ...) {
    object <- update.spacetimeobj(
        object = object,
        user_kappa = user_kappa,
        user_sigma = user_sigma,
        user_gamma = user_gamma,
        user_rho = user_rho
    )
    
    Q <- object$Q
    return(Q)
}


#' @name predict.spacetimeobj
#' @title Prediction of a space-time SPDE
#' @description The function is used for computing kriging predictions based
#' on data \eqn{Y_i = u(s_i,t_i) + \epsilon_i}, where \eqn{\epsilon}{\epsilon}
#' is mean-zero Gaussian measurement noise and \eqn{u(s,t)}{u(s,t)} is defined by
#' a spatio-temporal SPDE as described in [spacetime.operators()].
#' @param object The covariance-based rational SPDE approximation,
#' computed using [spacetime.operators()]
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
#' @param posterior_samples If `TRUE`, posterior samples will be returned.
#' @param n_samples Number of samples to be returned. Will only be used if `sampling` is `TRUE`.
#' @param only_latent Should the posterior samples be only given to the laten model?
#' @param ... further arguments passed to or from other methods.
#' @return A list with elements
#' \item{mean }{The kriging predictor (the posterior mean of u|Y).}
#' \item{variance }{The posterior variances (if computed).}
#' @export
#' @method predict spacetimeobj
#' @examples
#'s <- seq(from = 0, to = 20, length.out = 101)
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(space_loc = s, time_loc = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#'# generate data
#'sigma.e <- 0.01
#'n.obs <- 500
#'obs.loc <- data.frame(x = max(s)*runif(n.obs), 
#'                      t = max(t)*runif(n.obs))
#'A <- rSPDE.Ast(space_loc = s, time_loc = t, obs.s = obs.loc$x, obs.t = obs.loc$t)
#'Aprd <- Diagonal(dim(A)[2])
#'x <- simulate(op_cov, nsim = 1) 
#'Y <- A%*%x + sigma.e*rnorm(n.obs)
#'u.krig <- predict(op_cov, A, Aprd, Y, sigma.e)
predict.spacetimeobj <- function(object, A, Aprd, Y, sigma.e, mu = 0,
                                 compute.variances = FALSE, posterior_samples = FALSE,
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
    
   
    if (!no_nugget) {
        ## construct Q
        Q <- object$Q
        ## compute Q_x|y
        Q_xgiveny <- (t(A) %*% Q.e %*% A) + Q
        ## construct mu_x|y
        mu_xgiveny <- t(A) %*% Q.e %*% Y
            
        R <- Matrix::Cholesky(forceSymmetric(Q_xgiveny))
        mu_xgiveny <- solve(R, mu_xgiveny, system = "A")
            
        mu_xgiveny <- mu + mu_xgiveny
        out$mean <- Aprd %*% mu_xgiveny
            
        if (compute.variances) {
            out$variance <- diag(Aprd %*% solve(Q_xgiveny, t(Aprd)))
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
            post_cov <- Aprd %*% solve(Q_xgiveny, t(Aprd))
        } else {
            M <- Q - QiAt %*% solve(AQiA, t(QiAt))
            post_cov <- Aprd %*% M %*% t(Aprd)
        }
        Y_tmp <- as.matrix(Y)
        mean_tmp <- as.matrix(out$mean)
        out$samples <- lapply(1:ncol(Y_tmp), function(i) {
            Z <- rnorm(dim(post_cov)[1] * n_samples)
            dim(Z) <- c(dim(post_cov)[1], n_samples)
            LQ <- chol(forceSymmetric(post_cov))
            X <- LQ %*% Z
            X <- X + mean_tmp[, i]
            if (!only_latent) {
                X <- X + matrix(rnorm(n_samples * dim(Aprd)[1], sd = sigma.e), nrow = dim(Aprd)[1])
            }
            return(X)
        })
    }
    return(out)
}


#' @name spacetime.loglike
#' @title Object-based log-likelihood function for latent spatio-temporal
#' SPDE model 
#' @description This function evaluates the log-likelihood function for a
#' Gaussian process with a space-time SPDE covariance function, that is observed under
#' Gaussian measurement noise:
#' \eqn{Y_i = u(s_i,t_i) + \epsilon_i}{Y_i = u(s_i,t_i) + \epsilon_i}, where
#' \eqn{\epsilon_i}{\epsilon_i} are iid mean-zero Gaussian variables.
#' @param object The model object computed using [spacetime.operators()]
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location
#' to the finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0).
#' @param user_kappa If non-null, update the range parameter of the
#' covariance function.
#' @param user_sigma If non-null, update the standard deviation of
#' the covariance function.
#' @param user_gamma If non-null, update the temporal range parameter
#' of the covariance function.
#' @param user_rho If non-null, update the drift parameter of the
#' covariance function.
#' @return The log-likelihood value.
#' @noRd
#' @seealso [spacetime.operators()], [predict.spacetimeobj()]
#' @examples
#'s <- seq(from = 0, to = 20, length.out = 101)
#'t <- seq(from = 0, to = 20, length.out = 31)
#'
#'op_cov <- spacetime.operators(space_loc = s, time_loc = t,
#'                              kappa = 5, sigma = 10, alpha = 1,
#'                              beta = 2, rho = 1, gamma = 0.05)
#'# generate data
#'sigma.e <- 0.01
#'n.obs <- 500
#'obs.loc <- data.frame(x = max(s)*runif(n.obs), 
#'                      t = max(t)*runif(n.obs))
#'A <- rSPDE.Ast(space_loc = s, time_loc = t, obs.s = obs.loc$x, obs.t = obs.loc$t)
#'Y <- A%*%x + sigma.e*rnorm(n.obs)
#'spacetime.loglike(object, Y, A, sigma.e)
spacetime.loglike <- function(object, Y, A, sigma.e, mu = 0,
                              user_kappa = NULL,
                              user_sigma = NULL,
                              user_gamma = NULL,
                              user_rho = NULL) {
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
    
    object <- update.spacetimeobj(
        object = object,
        user_kappa = user_kappa,
        user_sigma = user_sigma,
        user_gamma = user_gamma,
        user_rho = user_rho)
    
    
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
    Q.R <- Matrix::Cholesky(Q)
        
    logQ <- 2 * c(determinant(Q.R, logarithm = TRUE, sqrt = TRUE)$modulus)
    
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
aux_lme_spacetime.loglike <- function(object, y, X_cov, repl, A_list, sigma_e, beta_cov) {
    l_tmp <- tryCatch(
        aux2_lme_spacetime.loglike(
            object = object,
            y = y, X_cov = X_cov, repl = repl, A_list = A_list,
            sigma_e = sigma_e, beta_cov = beta_cov
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
aux2_lme_spacetime.loglike <- function(object, y, X_cov, repl, A_list, sigma_e, beta_cov) {
    
    Q <- object$Q
    
    R <- Matrix::Cholesky(Q)
    
    prior.ld <- c(determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus)
    
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
        R.p <- Matrix::Cholesky(Q.p)
        
        posterior.ld <- c(determinant(R.p, logarithm = TRUE, sqrt = TRUE)$modulus)
 
        l <- l + prior.ld - posterior.ld - n.o * log(sigma_e)
        
        v <- y_
        
        if (ncol(X_cov) > 0) {
            X_cov_tmp <- X_cov_tmp[!na_obs, , drop = FALSE]
            # X_cov_tmp <- X_cov_list[[as.character(i)]]
            v <- v - X_cov_tmp %*% beta_cov
        }
        
        mu.p <- solve(R.p, as.vector(t(A_tmp) %*% v / sigma_e^2), system = "A")
        
        v <- v - A_tmp %*% mu.p
        
        l <- l - 0.5 * (t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
            0.5 * n.o * log(2 * pi)
    }
    
    return(as.double(l))
}



#' Observation matrix for space-time models
#'
#' @param mesh_space mesh object for models on 1d or 2d domains
#' @param space_loc mesh locations for models on 1d domains
#' @param mesh_time mesh object for time discretization
#' @param time_loc mesh locations for time discretization
#' @param graph MetricGraph object for models on metric graphs
#' @param obs.s spatial locations of observations
#' @param obs.t time points for observations
#'
#' @return Observation matrix linking observation locations to mesh nodes
#' @export
#'
#' @examples
#' s <- seq(from = 0, to = 20, length.out = 11)
#' t <- seq(from = 0, to = 20, length.out = 5)
#' n.obs <- 10
#' obs.loc <- data.frame(x = max(s)*runif(n.obs), 
#' t = max(t)*runif(n.obs))
#' A <- rSPDE.Ast(space_loc = s,time_loc = t, 
#'                obs.s = obs.loc$x, obs.t = obs.loc$t)
rSPDE.Ast <- function(mesh_space = NULL,
                      space_loc = NULL,
                      mesh_time = NULL,
                      time_loc = NULL,
                      graph = NULL,
                      obs.s = NULL, 
                      obs.t = NULL) {
    
    if ((!is.null(mesh_space) && !is.null(graph)) || (!is.null(mesh_space) && !is.null(space_loc)) || (!is.null(graph) && !is.null(space_loc))){
        stop("You should provide only one of mesh_space, space_loc or graph.")
    }
    
    if (is.null(mesh_space) && is.null(graph) && is.null(space_loc)) {
        stop("You must provide one of mesh_space, space_loc or graph.")
    }
    
    if (is.null(mesh_time) && is.null(time_loc)) {
        stop("You should provide mesh_time or time_loc.")
    }
    
    if(!is.null(mesh_time)){
        d_time <- get_inla_mesh_dimension(inla_mesh = mesh_time)
        if(d_time != 1) {
            stop("mesh_time should be a 1d mesh")
        }
        time <- mesh_time$loc
    } else {
        time <- time_loc
    }
    
    if(is.null(obs.s) || is.null(obs.t)) {
        stop("obs.s and obs.t must be provided")
    } 
    
    
    At <- rSPDE.A1d(time, obs.t)
    if(!is.null(graph)) {
        As <- graph$fem_basis(obs.s)
    } else if (!is.null(mesh_space)){
        As <- fm_basis(mesh_space, obs.s)    
    } else {
        space_loc <- as.matrix(space_loc)
        if(min(dim(space_loc))>1) {
            stop("For 2d domains, please provide mesh_space instead of space_loc")
        }
        mesh_space <- fm_mesh_1d(s)
        As <- fm_basis(mesh_space, obs.s)    
    }
    
    As.bar <- kronecker(matrix(1,nrow=1, ncol=length(time)),As)
    At.bar <- kronecker(At,matrix(1,nrow=1, ncol=dim(As)[2]))
    return(As.bar*At.bar)
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


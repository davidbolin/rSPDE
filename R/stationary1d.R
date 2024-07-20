#' Rational approximation of the Matern fields on intervals and metric graphs
#' 
#' The function is used for computing an approximation,
#' which can be used for inference and simulation, of the fractional SPDE
#' \deqn{(\kappa^2 - \Delta)^{\alpha/2} (\tau u(s)) = W}
#' on intervals or metric graphs. Here \eqn{W} is Gaussian white noise, 
#' \eqn{\kappa} controls the range, \eqn{\alpha = \nu + 1/2} with \eqn{\nu>0} 
#' controls the smoothness and \eqn{\tau} is related to the marginal variances
#' through 
#' \deqn{\sigma^2 = \frac{\Gamma(\nu)}{\tau^2\Gamma(\alpha)2\sqrt{\pi}\kappa^{2\nu}}.}
#' @param graph Metric graph object. The default is NULL, which means that a stationary
#' Matern model on the line is created. 
#' @param loc Locations where to evaluate the model.
#' @param bc Specifies the boundary conditions. The default is "free" which gives
#' stationary Matern models on intervals. Other options are "Neumann" or "Dirichlet".
#' @param m The order of the approximation
#' @param parameterization Which parameterization to use? `matern` uses range, std. deviation and nu 
#' (smoothness). `spde` uses kappa, tau and alpha. The default is `matern`.
#' @param kappa Range parameter
#' @param range practical correlation range
#' @param nu Smoothness parameter
#' @param sigma Standard deviation
#' @param tau Precision parameter
#' @param alpha Smoothness parameter
#' @param type_rational_approximation Method used to compute the coefficients of the rational approximation.
#' @param type_interp Interpolation method for the rational coefficients. 
#'
#' @return A model object for the the approximation
#' @export
#' @examples
#' s <- seq(from = 0, to = 1, length.out = 101)
#' kappa <- 20
#' sigma <- 2
#' nu <- 0.8
#' r <- sqrt(8*nu)/kappa #range parameter
#' op_cov <- matern.rational(loc = s, nu = nu, range = r, sigma = sigma, m = 2, 
#' parameterization = "matern")
#' cov.true <- matern.covariance(abs(s-s[1]), kappa = kappa, sigma = sigma, nu = nu)
#' cov.approx <- op_cov$covariance(ind = 1)
#' 
#' plot(s, cov.true)
#' lines(s, cov.approx, col = 2)
matern.rational = function(graph = NULL,
                           loc = NULL,
                           bc = c("free", "Neumann", "Dirichlet"),
                           kappa = NULL,
                           range = NULL,
                           nu = NULL, 
                           sigma = NULL,
                           tau = NULL,
                           alpha = NULL,
                           m = 2,
                           parameterization = c("matern", "spde"),
                           type_rational_approximation = "brasil", 
                           type_interp = "spline") {
    bc <- bc[[1]]
    has_graph <- FALSE
    if(!is.null(graph)) {
        has_graph <- TRUE
    }
    parameterization <- parameterization[[1]]
    
    if (!parameterization %in% c("matern", "spde")) {
        stop("parameterization should be either 'matern' or 'spde'!")
    }
    
    if (parameterization == "spde") {
        if (!is.null(nu)) {
            stop("For 'spde' parameterization, you should NOT supply 'nu'. You need to provide 'alpha'!")
        }
        if (!is.null(sigma)) {
            stop("For 'spde' parameterization, you should NOT supply 'sigma'. You need to provide 'tau'!")
        }
        if (!is.null(range)) {
            stop("For 'spde' parameterization, you should NOT supply 'range'. You need to provide 'kappa'!")
        }
    } else {
        if (!is.null(alpha)) {
            stop("For 'matern' parameterization, you should NOT supply 'alpha'. You need to provide 'nu'!")
        }
        if (!is.null(tau)) {
            stop("For 'matern' parameterization, you should NOT supply 'tau'. You need to provide 'sigma'!")
        }
        if (!is.null(kappa)) {
            stop("For 'matern' parameterization, you should NOT supply 'kappa'. You need to provide 'range'!")
        }
    }
    
    if (!is.null(graph)) {
        if (!inherits(graph, "metric_graph")) {
            stop("graph should be a metric_graph object!")
        }
    }
    
    if(is.null(graph) && !(bc == "free") ){
        stop("If boundary conditions are used, the domain must be specified.")
    }
    
    if (parameterization == "spde") {
        if (is.null(kappa) || is.null(tau)) {
            if (is.null(graph) && is.null(loc)) {
                stop("You should either provide all the parameters, or you should provide a graph or loc")
            }
            if (!is.null(loc)) {
                range_mesh <- max(loc) - min(loc)
                param <- get.initial.values.rSPDE(mesh.range = range_mesh, dim = 1, 
                                                  parameterization = parameterization, 
                                                  nu = nu)
            } else if (!is.null(graph)) {
                param <- get.initial.values.rSPDE(graph.obj = graph, 
                                                  parameterization = parameterization, 
                                                  nu = nu)
            }
        }
        
        if (is.null(kappa)) {
            kappa <- exp(param[2])
        } else {
            kappa <- rspde_check_user_input(kappa, "kappa", 0)
        }
        if (is.null(tau)) {
            tau <- exp(param[1])
        } else {
            tau <- rspde_check_user_input(tau, "tau", 0)
        }
        
        if (is.null(alpha)) {
            alpha <- 1 + 0.5
        } else {
            alpha <- rspde_check_user_input(alpha, "alpha", 0.5)
            #alpha <- min(alpha, 10)
        }
        
        nu <- alpha - 0.5
        
        range <- sqrt(8 * nu) / kappa
        sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2 * nu) *
                                       (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
    } else if (parameterization == "matern") {
        if (is.null(sigma) || is.null(range)) {
            if (is.null(graph) && is.null(loc)) {
                stop("You should either provide all parameters, or you should provide graph or loc.")
            }
            if (!is.null(loc)) {
                range_mesh <- max(loc) - min(loc)
                param <- get.initial.values.rSPDE(mesh.range = range_mesh, dim = 1, 
                                                  parameterization = parameterization, 
                                                  nu = nu)
            } else if (!is.null(graph)) {
                param <- get.initial.values.rSPDE(graph.obj = graph, 
                                                  parameterization = parameterization, 
                                                  nu = nu)
            }
        }
        
        if (is.null(range)) {
            range <- exp(param[2])
        } else {
            range <- rspde_check_user_input(range, "range", 0)
        }
        if (is.null(sigma)) {
            sigma <- exp(param[1])
        } else {
            sigma <- rspde_check_user_input(sigma, "sigma", 0)
        }
        
        if (is.null(nu)) {
            nu <- 1
        } else {
            nu <- rspde_check_user_input(nu, "nu", 0)
        }
        kappa <- sqrt(8 * nu) / range
        tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
                                     (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
        alpha <- nu + 1 / 2
    }
    
    
    output <- list(
        graph = graph,
        has_graph = has_graph,
        loc = loc,
        bc = bc,
        range = range,
        sigma = sigma,
        nu = nu, 
        kappa = kappa,
        tau = tau,
        alpha = alpha,
        m = m,
        d = 1,
        type_rational_approximation = type_rational_approximation, 
        type_interp = type_interp,
        parameterization = parameterization,
        stationary = TRUE
    )
    output$covariance <- function(ind = NULL) {
        
        if(is.null(loc)) {
            stop("The locations where to simulate must be specified in loc")
        }
        
        if(!is.null(graph)) {
            stop("Not implemented for graphs yet.")
        }
        
        if(!is.null(ind)) {
            h <- abs(loc - loc[ind])
        } else {
            h <- as.matrix(dist(loc))
        }
        
        Sigma <- matern.rational.cov(h = h, order = m,
                                     kappa = kappa,
                                     nu = nu, 
                                     sigma = sigma,
                                     type_rational = type_rational_approximation, 
                                     type_interp = type_interp)
        
        return(Sigma)
    }
    class(output) <- "rSPDEobj1d"
    return(output)
}



#' @name update.rSPDEobj1d
#' @title Update parameters of rSPDEobj1d objects
#' @description Function to change the parameters of a rSPDEobj1d object
#' @param object The covariance-based rational SPDE approximation,
#' computed using [matern.rational()]
#' @param user_kappa If non-null, update the parameter kappa of the SPDE. Will be used if parameterization is 'spde'.
#' @param user_tau If non-null, update the parameter tau of the SPDE. Will be used if parameterization is 'spde'.
#' @param user_sigma If non-null, update the standard deviation of
#' the covariance function. Will be used if parameterization is 'matern'.
#' @param user_range If non-null, update the range parameter
#' of the covariance function. Will be used if parameterization is 'matern'.
#' @param user_theta For non-stationary models. If non-null, update the vector of parameters.
#' @param user_nu If non-null, update the shape parameter of the
#' covariance function. Will be used if parameterization is 'matern'.
#' @param user_alpha If non-null, update the fractional SPDE order parameter. Will be used if parameterization is 'spde'.
#' @param user_m If non-null, update the order of the rational
#' approximation, which needs to be a positive integer.
#' @param graph An optional `metric_graph` object. 
#' @param loc The locations of interest for evaluating the model. 
#' @param parameterization If non-null, update the parameterization. 
#' @param type_rational_approximation Which type of rational
#' approximation should be used? The current types are "chebfun",
#' "brasil" or "chebfunLB".
#' @param ... Currently not used.
#' @return It returns an object of class "rSPDEobj1d". This object contains the
#' same quantities listed in the output of [matern.rational()].
#' @method update rSPDEobj1d
#' @seealso [simulate.rSPDEobj1d()], [matern.rational()]
#' @export
#' @examples
#' 
#' s <- seq(from = 0, to = 1, length.out = 101)
#' kappa <- 20
#' sigma <- 2
#' nu <- 0.8
#' r <- sqrt(8*nu)/kappa #range parameter
#' op_cov <- matern.rational(loc = s, nu = nu, range = r, sigma = sigma, m = 2, 
#' parameterization = "matern")
#' cov1 <- op_cov$covariance(ind = 1)
#' op_cov <- update(op_cov, user_range = 0.2)
#' cov2 <- op_cov$covariance(ind = 1)
#' plot(s, cov1, type = "l")
#' lines(s, cov2, col = 2)
update.rSPDEobj1d <- function(object, 
                              user_nu = NULL, 
                              user_alpha = NULL,
                              user_kappa = NULL,
                              user_tau = NULL,
                              user_sigma = NULL,
                              user_range = NULL,
                              user_theta = NULL,
                              user_m = NULL,
                              loc = NULL,
                              graph = NULL,
                              parameterization = NULL,
                              type_rational_approximation =
                                  object$type_rational_approximation,
                              ...) {
    new_object <- object
    d <- object$d
    
    if (is.null(parameterization)) {
        parameterization <- new_object$parameterization
    } else {
        parameterization <- parameterization[[1]]
        if (!parameterization %in% c("matern", "spde")) {
            stop("parameterization should be either 'matern' or 'spde'!")
        }
    }
    
    if (parameterization == "spde") {
        if (!is.null(user_kappa)) {
            new_object$kappa <- rspde_check_user_input(user_kappa, "kappa", 0)
            new_object$range <- NULL
            new_object$sigma <- NULL
        }
        
        if (!is.null(user_tau)) {
            new_object$tau <- rspde_check_user_input(user_tau, "tau", 0)
            new_object$sigma <- NULL
        }
        
        if (!is.null(user_alpha)) {
            alpha <- rspde_check_user_input(user_alpha, "alpha", d / 2)
            user_nu <- alpha - d / 2
            new_object$nu <- user_nu
            new_object$alpha <- alpha
        }
    } else if (parameterization == "matern") {
        if (!is.null(user_range)) {
            new_object$range <- rspde_check_user_input(user_range, "range", 0)
            new_object$kappa <- NULL
            new_object$tau <- NULL
        }
        
        if (!is.null(user_sigma)) {
            new_object$sigma <- rspde_check_user_input(user_sigma, "sigma", 0)
            new_object$tau <- NULL
        }
        if (!is.null(user_nu)) {
            new_object$nu <- rspde_check_user_input(user_nu, "nu")
        }
        alpha <- new_object$nu + d / 2
        new_object$alpha <- alpha
    }

    ## get parameters
    alpha <- new_object$alpha

    if (!is.null(user_m)) {
        new_object$m <- as.integer(rspde_check_user_input(user_m, "m", 0))
    }
    
    
    if (is.null(loc)) {
        loc <- new_object[["loc"]]
    }
    if (is.null(graph)) {
        graph <- new_object$graph
    }
    
    if (parameterization == "spde") {
        new_object <- matern.rational(
            kappa = new_object$kappa,
            tau = new_object$tau,
            alpha = new_object$alpha,
            m = new_object$m,
            loc = loc,
            graph = graph,
            parameterization = parameterization,
            type_rational_approximation = type_rational_approximation
        )
    } else {
        new_object <- matern.rational(
            sigma = new_object$sigma,
            range = new_object$range,
            nu = new_object$nu,
            m = new_object$m,
            loc = loc,
            graph = graph,
            parameterization = parameterization,
            type_rational_approximation = type_rational_approximation
        )
    }

    return(new_object)
}

#' @title Simulation of a Matern field using a rational SPDE approximation
#'
#' @description The function samples a Gaussian random field based on a
#' pre-computed rational SPDE approximation.
#'
#' @param object The rational SPDE approximation, computed
#' using [matern.rational()].
#' @param nsim The number of simulations.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param ... Currently not used.
#'
#' @return A matrix with the `n` samples as columns.
#' @seealso [matern.rational()]
#' @export
#' @method simulate rSPDEobj1d
#'
#' @examples
#' # Sample a Gaussian Matern process on R using a rational approximation
#' range <- 0.2
#' sigma <- 1
#' nu <- 0.8
#'
#' # compute rational approximation
#' x <- seq(from = 0, to = 1, length.out = 100)
#' op <- matern.rational(
#'   range = range, sigma = sigma,
#'   nu = nu, loc = x
#' )
#'
#' # Sample the model and plot the result
#' Y <- simulate(op)
#' plot(x, Y, type = "l", ylab = "u(x)", xlab = "x")
#'
simulate.rSPDEobj1d <- function(object,
                              nsim = 1,
                              seed = NULL,
                              ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    if (!inherits(object, "rSPDEobj1d")) {
        stop("input object is not of class rSPDEobj1d")
    }
    
    if(is.null(object$loc)) {
        stop("The locations where to simulate must be specified in loc")
    }
    
    if(!is.null(object$graph)) {
        stop("Simulation not implemented for graphs yet.")
    }
    tmp <- sort(object$loc, index.return = TRUE)
    reo <- tmp$ix
    loc_sort <- tmp$x
    ireo <- 1:length(loc_sort)
    ireo[reo] <- 1:length(loc_sort)
    
    ldl <- matern.rational.ldl(loc = loc_sort, order = object$m,
                               nu = object$nu, kappa = object$kappa,
                               sigma = object$sigma,
                               type_rational = object$type_rational,
                               type_interp = object$type_interp)
    m <- dim(ldl$L)[1]
    z <- rnorm(nsim * m)
    dim(z) <- c(m, nsim)
    x <- solve(sqrt(ldl$D), z)
    x <- ldl$A %*% solve(ldl$L,x)
    return(x[ireo,])
}


#' @noRd
aux2_lme_rSPDE.matern.rational.loglike <- function(object, y, X_cov, repl, 
                                                   loc, sigma_e, beta_cov) {
    
    repl_val <- unique(repl)
    l <- 0
    
    for (i in repl_val) {
        ind_tmp <- (repl %in% i)
        y_tmp <- y[ind_tmp]
        loc_ <- loc[ind_tmp]
        
        if (ncol(X_cov) == 0) {
            X_cov_tmp <- 0
        } else {
            X_cov_tmp <- X_cov[ind_tmp, , drop = FALSE]
        }
        na_obs <- is.na(y_tmp)
        
        y_ <- y_tmp[!na_obs]
        loc_ <- loc_[!na_obs]
        n.o <- length(y_)
        
        ls <- sort(loc_, index.return = TRUE)
        loc_ <- ls$x
        ## Construct prior Q
        tmp <- matern.rational.ldl(loc = ls$x, order = object$m, nu = object$nu, 
                                   kappa = object$kappa, sigma = object$sigma, 
                                   type_rational = object$type_rational_approx, 
                                   type_interp =  object$type_interp)    

        prior.ld <- 0.5 * sum(log(diag(tmp$D)))
        
        Q.post <- t(tmp$L)%*%tmp$D%*%tmp$L + t(tmp$A) %*% tmp$A / sigma_e^2
        
        R.post <- tryCatch(Matrix::Cholesky(Q.post), error = function(e) {
            return(NULL)
        })
        if (is.null(R.post)) {
            return(-10^100)
        }
        
        posterior.ld <- c(determinant(R.post, logarithm = TRUE, sqrt = TRUE)$modulus)
        
        v <- y_
        
        if (ncol(X_cov) != 0) {
            X_cov_tmp <- X_cov_tmp[!na_obs, , drop = FALSE]
            v <- v - X_cov_tmp %*% beta_cov
        }
        v <- v[ls$ix]
        AtY <- t(tmp$A) %*% v / sigma_e^2
        
        mu.post <- solve(R.post, AtY, system = "A")
        
        v1 <- tmp$L %*% mu.post
        v2 <- v - tmp$A %*% mu.post
        l <- l + prior.ld - posterior.ld - n.o * log(sigma_e)
        l <- l -  0.5 * (t(v1) %*% tmp$D %*% v1 + t(v2) %*% (v2) / sigma_e^2) -
            0.5 * n.o * log(2 * pi)
    }
    
    return(as.double(l))
}

#' @noRd
aux_lme_rSPDE.matern.rational.loglike <- function(object, y, X_cov, repl, loc, sigma_e, beta_cov) {
    
    l_tmp <- tryCatch(
        aux2_lme_rSPDE.matern.rational.loglike(
            object = object,
            y = y, X_cov = X_cov, repl = repl, loc = loc,
            sigma_e = sigma_e, beta_cov = beta_cov
        ),
        error = function(e) {
            return(NULL)
        }
    )
    if (is.null(l_tmp)) {
        return(-10^100)
    }
    #cat("ll =", l_tmp, "beta = ", beta_cov, ", sigma_e = ", sigma_e, ", kappa = ", object$kappa, ", tau = ", object$tau, ", nu = ", object$nu, "\n")
    return(l_tmp)
}


#' Precision matrix of rational approximation of Matern covariance
#'
#' Computes the precision matrix for a rational approximation of the Matern covariance on intervals.
#' @param loc Locations at which the precision is evaluated
#' @param order Order of the rational approximation
#' @param nu Smoothness parameter
#' @param kappa Range parameter
#' @param sigma Standard deviation
#' @param type_rational Method used to compute the coefficients of the rational approximation.
#' @param type_interp Interpolation method for the rational coefficients. 
#' @param equally_spaced Are the locations equally spaced?
#' @param ordering Return the matrices ordered by field or by location? 
#' @return A list containing the precision matrix `Q` of the process and its derivatives if they exist, and
#' a matrix `A` that extracts the elements corresponding to the process. 
#' @export 
#'
#' @examples
#' h <- seq(from = 0, to = 1, length.out = 100)
#' cov.true <- matern.covariance(h, kappa = 10, sigma = 1, nu = 0.8)
#' Q.approx <- matern.rational.precision(h, kappa = 10, sigma = 1, nu = 0.8, order = 2)
#' cov.approx <- Q.approx$A%*%solve(Q.approx$Q,Q.approx$A[1,])
#' 
#' plot(h, cov.true)
#' lines(h, cov.approx, col = 2)
matern.rational.precision <- function(loc,
                                      order,
                                      nu,
                                      kappa,
                                      sigma,
                                      type_rational = "brasil",
                                      type_interp = "spline",
                                      equally_spaced = FALSE,
                                      ordering = c("field", "location")) {
    ordering <- ordering[[1]]
    if(!(ordering %in% c("field", "location"))) {
        stop("Ordering must be 'field' or 'location'.")
    }
    if(is.matrix(loc) && min(dim(loc)) > 1) {
        stop("Only one dimensional locations supported.")
    }
    alpha=nu+1/2
    n <- length(loc)
    if (alpha%%1 == 0) {
        tmp <- matern.p.precision(loc = loc,kappa = kappa,p = 0,
                                  equally_spaced = equally_spaced, alpha = alpha) 
        Q <- tmp$Q
        A <- tmp$A
    } else {
        coeff <- interp_rational_coefficients(order = order, 
                                              type_rational_approx = type_rational, 
                                              type_interp = type_interp, 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
        
        ## k part
        tmp <- matern.k.precision(loc = loc,kappa,
                                  equally_spaced = equally_spaced, 
                                  alpha = alpha)
        Q <- tmp$Q/(k*sigma^2)
        A <- tmp$A
        
        ## p part
        for(i in 1:length(p)){
            tmp <- matern.p.precision(loc = loc, kappa = kappa, p =p[i],
                                      equally_spaced = equally_spaced, 
                                      alpha = alpha)
            Q <- bdiag(Q, tmp$Q/(r[i]*sigma^2))
            A = cbind(A,tmp$A)
        }
    }
    if(ordering == "location") {
        reo <- compute.reordering(n,order,alpha)
        Q <- Q[reo,reo]
        A <- A[,reo]
    } 
    return(list(Q = Q, A = A))    
}


#' LDL factorization of rational approximation of Matern covariance
#'
#' Computes the LDL factorization for a rational approximation of the Matern covariance on intervals.
#' @param loc Locations at which the precision is evaluated
#' @param order Order of the rational approximation
#' @param nu Smoothness parameter
#' @param kappa Range parameter
#' @param sigma Standard deviation
#' @param type_rational Method used to compute the coefficients of the rational approximation.
#' @param type_interp Interpolation method for the rational coefficients. 
#' @param equally_spaced Are the locations equally spaced?
#' @param ordering Return the matrices ordered by field or by location? 
#' @return A list containing the precision matrix `Q` of the process and its derivatives if they exist, and
#' a matrix `A` that extracts the elements corresponding to the process. 
#' @export 
#'
#' @examples
#' h <- seq(from = 0, to = 1, length.out = 100)
#' cov.true <- matern.covariance(h, kappa = 10, sigma = 1, nu = 0.8)
#' Q.approx <- matern.rational.precision(h, kappa = 10, sigma = 1, nu = 0.8, order = 2)
#' cov.approx <- Q.approx$A%*%solve(Q.approx$Q,Q.approx$A[1,])
#' 
#' plot(h, cov.true)
#' lines(h, cov.approx, col = 2)
matern.rational.ldl <- function(loc,
                                order,
                                nu,
                                kappa,
                                sigma,
                                type_rational = "brasil",
                                type_interp = "spline",
                                equally_spaced = FALSE,
                                ordering = c("field", "location")) {
    ordering <- ordering[[1]]
    if(!(ordering %in% c("field", "location"))) {
        stop("Ordering must be 'field' or 'location'.")
    }
    if(is.matrix(loc) && min(dim(loc)) > 1) {
        stop("Only one dimensional locations supported.")
    }
    alpha=nu+1/2
    n <- length(loc)
    if (alpha%%1 == 0) {
        tmp <- matern.p.chol(loc = loc,kappa = kappa,p = 0,
                             equally_spaced = equally_spaced, alpha = alpha) 
        L <- tmp$Bs
        D <- tmp$Fsi
        A <- tmp$A
    } else {
        coeff <- interp_rational_coefficients(order = order, 
                                              type_rational_approx = type_rational, 
                                              type_interp = type_interp, 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
        
        ## k part
        tmp <- matern.k.chol(loc = loc,kappa,equally_spaced = equally_spaced, 
                             alpha = alpha)
        L <- tmp$Bs
        D <- tmp$Fsi/(k*sigma^2)
        A <- tmp$A
        
        ## p part
        t.p <- rep(0,length(p))
        for(i in 1:length(p)){
            tmp <- matern.p.chol(loc = loc, kappa = kappa, p =p[i],
                                 equally_spaced = equally_spaced, 
                                 alpha = alpha)
            
            L <- bdiag(L, tmp$Bs)
            D <- bdiag(D, tmp$Fsi/(r[i]*sigma^2))
            A = cbind(A,tmp$A)
        }
    }
    if(ordering == "location") {
        stop("not implemented")
        #reo <- compute.reordering(n,order,alpha)
        #Q <- Q[reo,reo]
        #A <- A[,reo]
    } 
    return(list(L = L, D=D, A = A))    
}

### Utility below

# Reorder matern
#' @noRd
compute.reordering <- function(n,m,alpha) {
    if(alpha <0 || alpha > 2) {
        stop("only 0<alpha<2 supported")
    }
    if(alpha < 1) {
        return(as.vector(rep(seq(from=1,to=(m+1)*n,by=n),n) + kronecker(0:(n-1),rep(1,m+1))))
    } else {
        tmp <- rep(c(1,c(rbind(1+n*seq(from=1,to=2*m,by=2), 2+n*seq(from=1,to=2*m,by=2)))),n) 
        return(tmp + as.vector(rbind(0:(n-1), matrix(rep(1,2*m*n),2*m,n)%*%Diagonal(n,2*c(0:(n-1))))))
    }
}

# Derivatives of the matern covariance
matern.derivative = function(h, kappa, nu, sigma, deriv = 1) 
{
    if(deriv == 0) {
        C = matern.covariance(h, kappa = kappa, nu = nu, sigma = sigma)
        return(C)
    } else if (deriv == 1) {
        C = h * matern.covariance(h, kappa = kappa, nu = nu - 1, sigma = sigma)
        C[h == 0] = 0
        return(-(kappa^2/(2 * (nu - 1))) * C)
    }
    else if (deriv == 2) {
        C = matern.covariance(h, kappa = kappa, nu = nu - 1, sigma = sigma) + 
            h * matern.derivative(h, kappa = kappa, nu = nu - 1, 
                                  sigma = sigma, deriv = 1)
        return(-(kappa^2/(2 * (nu - 1))) * C)
    }
    else {
        C = (deriv - 1) * matern.derivative(h, kappa = kappa, 
                                            nu = nu - 1, sigma = sigma, 
                                            deriv = deriv - 2) + 
            h * matern.derivative(h, kappa = kappa, nu = nu - 
                                      1, sigma = sigma, deriv = deriv - 1)
        return(-(kappa^2/(2 * (nu - 1))) * C)
    }
    
}


# Precision for exponential covariance
#' @noRd
exp_precision <- function(loc, kappa, boundary = "free") {
    n <- length(loc)
    l <- diff(loc)
    
    # Precompute reusable terms
    a <- exp(-2 * kappa * l)
    b <- 1 - a
    
    # Initialize matrix Q efficiently
    Q <- Matrix(0, nrow = n, ncol = n)
    
    # Set diagonal elements
    diag(Q)[1:(n-1)] <- 1/2 + a/b
    diag(Q)[2:n] <- diag(Q)[2:n] + 1/2 + a/b
    
    # Set off-diagonal elements
    #  indices <- c(which(row(Q) == col(Q) + 1), which(row(Q) == col(Q) - 1))
    # off_diag_values <- -exp(-kappa * l) /b
    # Q[indices] <- off_diag_values
    
    
    row_indices <- seq_len(n - 1)
    col_indices <- row_indices + 1  # Offsets for off-diagonal elements
    
    # Compute the off-diagonal values
    off_diag_values <- -exp(-kappa * l)/b
    
    # Set the off-diagonal elements in matrix Q
    Q[cbind(row_indices, col_indices)] <- off_diag_values
    Q[cbind(col_indices, row_indices)] <- off_diag_values
    
    # Adjust boundary conditions if required
    if (boundary == "free") {
        Q[1, 1] <- Q[1, 1] + 0.5
        Q[n, n] <- Q[n, n] + 0.5
    }
    
    return(2 * kappa * Q[])
}







#Joint covariance of process and derivative for shifted Matern
matern.p.joint <- function(s,t,kappa,p, alpha = 1){
    
    if(alpha%%1 == 0) {
        fa <- alpha
    } else {
        fa <- floor(alpha) + 1    
    }
    mat <- matrix(0, nrow = fa, ncol = fa)
    for(i in 1:fa) {
        for(j in i:fa) {
            if(i==j) {
                mat[i,i] <- ((-1)^(i-1))*matern.p.deriv(s,t,kappa,p,alpha, deriv = 2*(i-1))
            } else {
                tmp <- matern.p.deriv(s,t,kappa,p,alpha, deriv = i-1 + j - 1)
                mat[i,j] <- (-1)^(j-1)*tmp
                mat[j,i] <- (-1)^(i-1)*tmp
            }
        }
    }
    
    return(mat)
}


matern.p <- function(s,t,kappa,p,alpha){
    h <- s-t
    if(p==0){
        return(matern.covariance(h, kappa = kappa, nu = alpha - 1/2, sigma = 1))
    } else {
        ca <- gamma(alpha)/gamma(alpha-0.5)
        fa <- floor(alpha)
        kp <- kappa*sqrt(1-p)
        out <- matern.covariance(h, kappa = kp, nu = 1/2, 
                                 sigma = sqrt(ca*sqrt(pi)/sqrt(1-p)))
        if(alpha < 1) {
            return(out)
        } else {
            
            for(j in 1:fa) {
                out <- out - matern.covariance(h, kappa = kappa, nu = j-1/2, 
                                               sigma = sqrt(ca*gamma(j-0.5)/gamma(j)))/p^(1 - j)
            }
            out <- out/p^fa
            return(out)    
        }
    }
}

matern.p.deriv <- function(s,t,kappa,p,alpha,deriv = 0){
    h <- s-t
    if(deriv ==0){
        return(matern.p(s,t,kappa,p,alpha))
    } else {
        if(p==0){
            return(matern.derivative(h, kappa = kappa, nu = alpha - 1/2, 
                                     sigma = 1, deriv = deriv))
        } else {
            ca <- gamma(alpha)/gamma(alpha-0.5)
            fa <- floor(alpha)
            kp <- kappa*sqrt(1-p)
            out <- matern.derivative(h, kappa = kp, nu =  1/2, 
                                     sigma = sqrt(ca*sqrt(pi)/sqrt(1-p)), deriv = deriv)
            if(alpha < 1) {
                return(out)
            } else {
                out <- out/p^fa
                for(j in 1:fa) {
                    out <- out - matern.derivative(h, kappa = kappa, nu = j-1/2, 
                                                   sigma = sqrt(ca*gamma(j-0.5)/gamma(j)), 
                                                   deriv = deriv)/p^(fa + 1 - j)
                }
            }
            return(out)
        }    
    }
    
}

matern.p.precision <- function(loc,kappa,p,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    if(alpha==1) {
        Q <- exp_precision(loc,kappa)/(2*kappa)
        A <- Diagonal(n)
        
        return(list(Q=Q,A=A))
    } else {
        if(alpha%%1 == 0) {
            fa <- alpha
        } else {
            fa <- floor(alpha) + 1    
        }
        
        if(fa == 1) {
            N <- n  + n - 1 
        } else {
            N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
        }
        
        ii <- numeric(N)
        jj <- numeric(N)
        val <- numeric(N)
        
        if(equally_spaced){
            Q.1 <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                                     matern.p.joint(loc[1],loc[1+1],kappa,p,alpha)),
                               cbind(matern.p.joint(loc[1+1],loc[1],kappa,p,alpha), 
                                     matern.p.joint(loc[1+1],loc[1+1],kappa,p,alpha))))
            
            Q.i <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                                     matern.p.joint(loc[1],loc[2],kappa,p,alpha),
                                     matern.p.joint(loc[1],loc[3],kappa,p,alpha)),
                               cbind(matern.p.joint(loc[2],loc[1],kappa,p,alpha),
                                     matern.p.joint(loc[2],loc[2],kappa,p,alpha),
                                     matern.p.joint(loc[2],loc[3],kappa,p,alpha)),
                               cbind(matern.p.joint(loc[3],loc[1],kappa,p,alpha),
                                     matern.p.joint(loc[3],loc[2],kappa,p,alpha),
                                     matern.p.joint(loc[3],loc[3],kappa,p,alpha))))[-c(1:fa),-c(1:fa)]
        }
        
        
        for(i in 1:max((n-1),1)){
            if(i==1){
                if(!equally_spaced){
                    Q.1 <- solve(rbind(cbind(matern.p.joint(loc[i],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i],loc[i+1],kappa,p,alpha)),
                                       cbind(matern.p.joint(loc[i+1],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i+1],loc[i+1],kappa,p,alpha))))
                    
                } 
                counter <- 1
                for(ki in 1:fa) {
                    for(kj in ki:(2*fa)) {
                        ii[counter] <- ki
                        jj[counter] <- kj
                        val[counter] <- Q.1[ki,kj]
                        counter <- counter + 1
                    }
                }
            } else {
                if(!equally_spaced){
                    Q.i <- solve(rbind(cbind(matern.p.joint(loc[i-1],loc[i-1],kappa,p,alpha),
                                             matern.p.joint(loc[i-1],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i-1],loc[i+1],kappa,p,alpha)),
                                       cbind(matern.p.joint(loc[i],loc[i-1],kappa,p,alpha),
                                             matern.p.joint(loc[i],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i],loc[i+1],kappa,p,alpha)),
                                       cbind(matern.p.joint(loc[i+1],loc[i-1],kappa,p,alpha),
                                             matern.p.joint(loc[i+1],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i+1],loc[i+1],kappa,p,alpha))))[-c(1:fa),-c(1:fa)]
                } 
                # Q[(2*n-1):(2*n), (2*n-3):(2*n)] = Q.i[3:4,]
                for(ki in 1:fa){
                    for(kj in ki:(2*fa)){
                        ii[counter] <- fa*(i-1) + ki
                        jj[counter] <- fa*(i-1) + kj
                        val[counter] <-Q.i[ki,kj]
                        counter <- counter + 1
                    }
                }     
            }
        }
        if(n<=2){
            Q.i <- Q.1
        }
        for(ki in 1:fa){
            for(kj in ki:fa){
                ii[counter] <- fa*(n-1) + ki
                jj[counter] <- fa*(n-1) + kj
                val[counter] <-Q.i[ki+fa,kj+fa]
                counter <- counter + 1
            }
        }    
        Q <- Matrix::sparseMatrix(i   = ii,
                                  j    = jj,
                                  x    = val,
                                  dims = c(fa*n, fa*n), symmetric=TRUE)
        
        A <-  Matrix::sparseMatrix(i   = 1:n,
                                   j    = seq(from=1,to=n*fa,by=fa),
                                   x    = rep(1,n),
                                   dims = c(n, fa*n))
        return(list(Q=Q,A=A))    
    }
}


matern.p.chol <- function(loc,kappa,p,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    
    if(alpha%%1 == 0) {
        fa <- alpha
    } else {
        fa <- floor(alpha) + 1    
    }
    
    #Bs <- Fs <- Fsi <- Diagonal(fa*n)
    if(fa == 1) {
        N <- n  + n - 1 
    } else {
        N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
    }
    
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    Sigma <- matrix(0,nrow=2*fa, ncol = 2*fa)
    Stransp <- outer((-1)^(0:(fa-1)),(-1)^(0:(fa-1)))
    Sdiag <- matern.p.joint(0,0,kappa,p,alpha)
    
    Sod <- matern.p.joint(loc[1],loc[2],kappa,p,alpha)
    di <- abs(loc[2]-loc[1])
    
    Sigma[1:fa,1:fa] <- Sdiag
    Sigma[(fa+1):(2*fa),(fa+1):(2*fa)] <- Sdiag
    Sigma[1:fa,(fa+1):(2*fa)] <- Sod
    Sigma[(fa+1):(2*fa),1:fa] <- Stransp*Sod
    
    
    Fs.d <- Fsi.d <- numeric(fa*n)
    
    if(equally_spaced){
        Sigma <- rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                             matern.p.joint(loc[1],loc[2],kappa,p,alpha)),
                       cbind(matern.p.joint(loc[2],loc[1],kappa,p,alpha),
                             matern.p.joint(loc[2],loc[2],kappa,p,alpha)))
    }
    
    for(i in 1:n){
        if(i==1){
            Fs.d[1] <- Sdiag[1,1]
            Fsi.d[1] <- 1/Sdiag[1,1]
            val[1] <- 1
            ii[1] <- 1
            jj[1] <- 1
            counter <- 2
            counter2 <- 1
            if(fa > 1) {
                for(k in 2:fa) {
                    if(0) { #k > 2
                        ss <- Sdiag[1:(k-1),1:(k-1)]
                        prec <- diag(1/sqrt(diag(ss)))
                        tmp <- prec%*%solve(prec%*%ss%*%prec, prec%*%Sdiag[1:(k-1),k])  
                    } else {
                        tmp <- solve(Sdiag[1:(k-1),1:(k-1)], Sdiag[1:(k-1),k])   
                    }
                    
                    
                    val[counter2 + (1:k)] <- c(-t(tmp),1)
                    ii[counter2 + (1:k)] <- rep(counter,k)
                    jj[counter2 + (1:k)] <- (counter-k+1):counter
                    Fs.d[k] <- Sdiag[k,k] - t(Sdiag[1:(k-1),k])%*%tmp
                    Fsi.d[k] <- 1/Fs.d[k]
                    counter <- counter + 1
                    counter2 <- counter2 + k
                }    
            }
        } else {
            if(!equally_spaced){
                #Sigma <- rbind(cbind(matern.p.joint(loc[i-1],loc[i-1],kappa,p,alpha),
                #                     matern.p.joint(loc[i-1],loc[i],kappa,p,alpha)),
                #               cbind(matern.p.joint(loc[i],loc[i-1],kappa,p,alpha),
                #                     matern.p.joint(loc[i],loc[i],kappa,p,alpha)))
                if(!(di==abs(loc[i]-loc[i-1]))) {
                    di=abs(loc[i]-loc[i-1])
                    Sod <- matern.p.joint(loc[i-1],loc[i],kappa,p,alpha)
                    Sigma[1:fa,(fa+1):(2*fa)] <- Sod
                    Sigma[(fa+1):(2*fa),1:fa] <- Stransp*Sod
                }

            } 
            for(k in (fa+1):(2*fa)) {
                if(0) { #k > 2
                    ss <- Sigma[1:(k-1),1:(k-1)]
                    prec <- diag(1/sqrt(diag(ss)))
                    tmp <- prec%*%solve(prec%*%ss%*%prec, prec%*%Sigma[1:(k-1),k])    
                } else {
                    tmp <- solve(Sigma[1:(k-1),1:(k-1)], Sigma[1:(k-1),k])   
                }
                
                val[counter2 + (1:k)] <- c(-tmp,1)
                ii[counter2 + (1:k)] <- rep(counter,k)
                jj[counter2 + (1:k)] <- (counter-k+1):counter
                Fs.d[counter] <- Sigma[k,k] - sum(Sigma[1:(k-1),k]*tmp)#t(Sigma[1:(k-1),k])%*%tmp
                Fsi.d[counter] <- 1/Fs.d[counter]
                counter <- counter + 1
                counter2 <- counter2 + k
            }    
        }
    }
    
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(fa*n, fa*n))
    Fs <-  Matrix::Diagonal(fa*n,Fs.d)
    Fsi <-  Matrix::Diagonal(fa*n,Fsi.d)
    A <-  Matrix::sparseMatrix(i   = 1:n,
                               j    = seq(from=1,to=n*fa,by=fa),
                               x    = rep(1,n),
                               dims = c(n, fa*n))
    return(list(Bs=Bs, Fs = Fs, Fsi = Fsi, A=A))    
}

matern.k.chol <- function(loc,kappa,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    
    fa <- floor(alpha)
    if(fa == 0) {
        N <- n 
    } else {
        N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
    }
    fa <- max(fa,1)
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    Sigma <- matrix(0, nrow = 2*fa,ncol = 2*fa)
    Stransp <- outer((-1)^(0:(fa-1)),(-1)^(0:(fa-1)))
    Sdiag <- matern.k.joint(0,0,kappa,alpha)
    Sigma[1:fa,1:fa] <- Sdiag
    Sigma[(fa+1):(2*fa),(fa+1):(2*fa)] <- Sdiag
    Fs.d <- Fsi.d <- numeric(fa*n)
    
    if(equally_spaced){
        Sigma <- rbind(cbind(matern.k.joint(loc[1],loc[1],kappa,alpha), 
                             matern.k.joint(loc[1],loc[2],kappa,alpha)),
                       cbind(matern.k.joint(loc[2],loc[1],kappa,alpha),
                             matern.k.joint(loc[2],loc[2],kappa,alpha)))
    }
    
    for(i in 1:n){
        if(i==1){
            
            Fs.d[1] <- Sdiag[1,1]
            Fsi.d[1] <- 1/Sdiag[1,1]
            val[1] <- 1
            ii[1] <- 1
            jj[1] <- 1
            counter <- 2
            counter2 <- 1
            if(fa > 1) {
                for(k in 2:fa) {
                    #tmp <- solve(Sigma.1[1:(k-1),1:(k-1)], Sigma.1[1:(k-1),k])
                    if(0) { #k > 2
                        ss <- Sdiag[1:(k-1),1:(k-1)]
                        prec <- diag(1/sqrt(diag(ss)))
                        tmp <- prec%*%solve(prec%*%ss%*%prec, prec%*%Sdiag[1:(k-1),k])  
                    } else {
                        tmp <- solve(Sdiag[1:(k-1),1:(k-1)], Sdiag[1:(k-1),k])   
                    }
                    
                    val[counter2 + (1:k)] <- c(-t(tmp),1)
                    ii[counter2 + (1:k)] <- rep(counter,k)
                    jj[counter2 + (1:k)] <- (counter-k+1):counter
                    Fs.d[k] <- Sdiag[k,k] - t(Sdiag[1:(k-1),k])%*%tmp
                    Fsi.d[k] <- 1/Fs.d[k]
                    counter <- counter + 1
                    counter2 <- counter2 + k
                }    
            }
        } else {
            if(!equally_spaced){
                
                #Sigma <- rbind(cbind(matern.k.joint(loc[i-1],loc[i-1],kappa,alpha),
                #                     matern.k.joint(loc[i-1],loc[i],kappa,alpha)),
                #               cbind(matern.k.joint(loc[i],loc[i-1],kappa,alpha),
                #                     matern.k.joint(loc[i],loc[i],kappa,alpha)))
                
                Sigma[1:fa,(fa+1):(2*fa)] <- matern.k.joint(loc[i-1],loc[i],kappa,alpha)
                Sigma[(fa+1):(2*fa),1:fa] <- Stransp*Sigma[1:fa,(fa+1):(2*fa)]
            } 
            for(k in (fa+1):(2*fa)) {
                #tmp <- solve(Sigma[1:(k-1),1:(k-1)], Sigma[1:(k-1),k])
                if(0) { #k > 2
                    ss <- Sigma[1:(k-1),1:(k-1)]
                    prec <- diag(1/sqrt(diag(ss)))
                    tmp <- prec%*%solve(prec%*%ss%*%prec, prec%*%Sigma[1:(k-1),k])  
                } else {
                    tmp <- solve(Sigma[1:(k-1),1:(k-1)], Sigma[1:(k-1),k])   
                }
                val[counter2 + (1:k)] <- c(-t(tmp),1)
                ii[counter2 + (1:k)] <- rep(counter,k)
                jj[counter2 + (1:k)] <- (counter-k+1):counter
                Fs.d[counter] <- Sigma[k,k] - sum(Sigma[1:(k-1),k]*tmp)#t(Sigma[1:(k-1),k])%*%tmp
                Fsi.d[counter] <- 1/Fs.d[counter]
                counter <- counter + 1
                counter2 <- counter2 + k
            }    
        }
    }
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(fa*n, fa*n))
    Fs <-  Matrix::Diagonal(fa*n,Fs.d)
    Fsi <-  Matrix::Diagonal(fa*n,Fsi.d)
    A <-  Matrix::sparseMatrix(i   = 1:n,
                               j    = seq(from=1,to=n*fa,by=fa),
                               x    = rep(1,n),
                               dims = c(n, fa*n))
    return(list(Bs=Bs, Fs = Fs, Fsi = Fsi, A=A))    
}




matern.p.chol.pred <- function(loc,loc.obs,kappa,p,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    
    if(alpha%%1 == 0) {
        fa <- alpha
    } else {
        fa <- floor(alpha) + 1    
    }
    
    #Bs <- Fs <- Fsi <- Diagonal(fa*n)
    if(fa == 1) {
        N <- 2*n
    } else {
        N <- 2*n*fa
    }
    
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    Sigma <- matrix(0,nrow=2*fa, ncol = 2*fa)
    Stransp <- outer((-1)^(0:(fa-1)),(-1)^(0:(fa-1)))
    Sdiag <- matern.p.joint(0,0,kappa,p,alpha)
    
    Sod <- matern.p.joint(loc[1],loc[2],kappa,p,alpha)
    di <- abs(loc[2]-loc[1])
    
    Sigma[1:fa,1:fa] <- Sdiag
    Sigma[(fa+1):(2*fa),(fa+1):(2*fa)] <- Sdiag
    Sigma[1:fa,(fa+1):(2*fa)] <- Sod
    Sigma[(fa+1):(2*fa),1:fa] <- Stransp*Sod
    
    
    Fs.d <- Fsi.d <- numeric(fa*n)
    
    for(i in 1:n){
        n1 <- which(loc.obs <  loc[i])
        n1 <- n1[length(n1)]
        n2 <- which(loc.obs >  loc[i])[1]
        if(is.na(n1) || is.na(n2)) {
            Sigma.nn <- Sdiag
            if(is.na(n1)) {
                Sigma.in <- matern.p.joint(loc[i],loc.obs[n2],kappa,p,alpha) 
                nn <- n2
            } else {
                Sigma.in <- matern.p.joint(loc.obs[n1],loc[i],kappa,p,alpha)    
                nn <- n1
            }
            tmp <- t(solve(Sigma.nn,Sigma.in))
            for(k in 1:fa) {
                val[counter2 + (1:fa)] <- tmp[k,]
                ii[counter2 + (1:fa)] <- rep(counter,fa)
                jj[counter2 + (1:fa)] <- 1+nn*(0:(fa-1))
                counter <- counter + 1
                counter2 <- counter2 + 2*fa
            }
        } else {
            Sod <- matern.p.joint(loc.obs[n1],loc.obs[n2],kappa,p,alpha)
            Sigma[1:fa,(fa+1):(2*fa)] <- Sod
            Sigma[(fa+1):(2*fa),1:fa] <- Stransp*Sod
            Sigma.nn <- Sigma
            Sigma.in <- cbind(matern.p.joint(loc[i],loc.obs[n1],kappa,p,alpha),matern.p.joint(loc[i],loc.obs[n2],kappa,p,alpha))            
            tmp <- t(solve(Sigma.nn,Sigma.in))
            
            for(k in 1:fa) {
                val[counter2 + (1:(2*fa))] <- tmp[k,]
                ii[counter2 + (1:(2*fa))] <- rep(counter,2*fa)
                jj[counter2 + (1:(2*fa))] <- c(1+n1*(0:(fa-1)),1+n2*(0:(fa-1)))
                counter <- counter + 1
                counter2 <- counter2 + 2*fa
            }
        }
    }
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(fa*n, length(loc.obs)*fa))
    
    return(Bs)    
}


matern.k.precision <- function(loc,kappa,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    if(alpha==1) {
        Q <- exp_precision(loc,kappa)/(2*kappa)
        A <- Diagonal(n)
        return(list(Q=Q,A=A))
    } else {
        fa <- floor(alpha)
        if(fa == 0) {
            N <- n 
        } else {
            N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
        }
        da <- max(fa,1)
        ii <- numeric(N)
        jj <- numeric(N)
        val <- numeric(N)
        
        if(equally_spaced){
            Q.1 <- solve(rbind(cbind(matern.k.joint(loc[1],loc[1],kappa,alpha), 
                                     matern.k.joint(loc[1],loc[1+1],kappa,alpha)),
                               cbind(matern.k.joint(loc[1+1],loc[1],kappa,alpha), 
                                     matern.k.joint(loc[1+1],loc[1+1],kappa,alpha))))
            
            Q.i <- solve(rbind(cbind(matern.k.joint(loc[1],loc[1],kappa,alpha), 
                                     matern.k.joint(loc[1],loc[2],kappa,alpha),
                                     matern.k.joint(loc[1],loc[3],kappa,alpha)),
                               cbind(matern.k.joint(loc[2],loc[1],kappa,alpha),
                                     matern.k.joint(loc[2],loc[2],kappa,alpha),
                                     matern.k.joint(loc[2],loc[3],kappa,alpha)),
                               cbind(matern.k.joint(loc[3],loc[1],kappa,alpha),
                                     matern.k.joint(loc[3],loc[2],kappa,alpha),
                                     matern.k.joint(loc[3],loc[3],kappa,alpha))))[-c(1:da),-c(1:da)]
            
        }
    }
   
    for(i in 1:max((n-1),1)){
        if(i==1){
            if(!equally_spaced){
                Q.1 <- solve(rbind(cbind(matern.k.joint(loc[i],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i],loc[i+1],kappa,alpha)),
                                   cbind(matern.k.joint(loc[i+1],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i+1],loc[i+1],kappa,alpha))))
                
            } 
            counter <- 1
            for(ki in 1:da) {
                for(kj in ki:(2*da)) {
                    ii[counter] <- ki
                    jj[counter] <- kj
                    val[counter] <- Q.1[ki,kj]
                    counter <- counter + 1
                }
            }
        } else {
            if(!equally_spaced){
                Q.i <- solve(rbind(cbind(matern.k.joint(loc[i-1],loc[i-1],kappa,alpha),
                                         matern.k.joint(loc[i-1],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i-1],loc[i+1],kappa,alpha)),
                                   cbind(matern.k.joint(loc[i],loc[i-1],kappa,alpha),
                                         matern.k.joint(loc[i],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i],loc[i+1],kappa,alpha)),
                                   cbind(matern.k.joint(loc[i+1],loc[i-1],kappa,alpha),
                                         matern.k.joint(loc[i+1],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i+1],loc[i+1],kappa,alpha))))[-c(1:da),-c(1:da)]
                
            } 
            # Q[(2*n-1):(2*n), (2*n-3):(2*n)] = Q.i[3:4,]
            for(ki in 1:da){
                for(kj in ki:(2*da)){
                    ii[counter] <- da*(i-1) + ki
                    jj[counter] <- da*(i-1) + kj
                    val[counter] <- Q.i[ki,kj]
                    counter <- counter + 1
                }
            }     
        }
    }
    if(n <= 2){
        Q.i <- Q.1
    }
        for(ki in 1:da){
            for(kj in ki:da){
                ii[counter] <- da*(n-1) + ki
                jj[counter] <- da*(n-1) + kj
                val[counter] <- Q.i[ki+da,kj+da]
                counter <- counter + 1
            }
        }

    Q <- Matrix::sparseMatrix(i   = ii,
                              j    = jj,
                              x    = val,
                              dims = c(da*n, da*n), symmetric=TRUE)
    
    A <-  Matrix::sparseMatrix(i   = 1:n,
                               j    = seq(from=1,to=n*da,by=da),
                               x    = rep(1,n),
                               dims = c(n, da*n))
    return(list(Q=Q,A=A))    
}

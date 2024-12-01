#' Summarise rSPDE objects without FEM
#'
#' Summary method for class "rSPDEobj1d"
#'
#' @param object an object of class "rSPDEobj1d", usually, a result of a call
#'   to [matern.rational()].
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary rSPDEobj1d
summary.rSPDEobj1d <- function(object, ...) {
    out <- list()
    class(out) <- "summary.rSPDEobj1d"
    
    out$kappa <- object$kappa
    out$sigma <- object$sigma
    out$nu <- object$nu
    out$range <- object$range
    out$tau <- object$tau
    out$alpha <- object$alpha
    out$order <- object$order
    out$parameterization <- object$parameterization
    return(out)
}

#' @param x an object of class "summary.rSPDEobj1d", usually, a result of a call
#'   to [summary.rSPDEobj1d()].
#' @export
#' @method print summary.rSPDEobj1d
#' @rdname summary.rSPDEobj1d
print.summary.rSPDEobj1d <- function(x, ...) {
    
    cat("Parameterization: ", x$parameterization, "\n")
    
    if(x$parameterization == "matern") {
        cat("Parameters of covariance function: range = ", x$range, ", sigma = ", x$sigma, 
            ", nu = ", x$nu, "\n")    
    } else {
        cat("Parameters of covariance function: kappa = ", x$kappa, ", tau = ", x$tau, 
            ", alpha = ", x$alpha, "\n")    
    }
    
    cat("Order or rational approximation: ", x$order, "\n")
}

#' @export
#' @method print rSPDEobj1d
#' @rdname summary.rSPDEobj1d
print.rSPDEobj1d <- function(x, ...) {
    print.summary.rSPDEobj1d(summary(x))
}




#' Summarise rSPDE objects
#'
#' Summary method for class "rSPDEobj"
#'
#' @param object an object of class "rSPDEobj", usually, a result of a call
#'   to [fractional.operators()], [matern.operators()], or
#'   [spde.matern.operators()].
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary rSPDEobj
summary.rSPDEobj <- function(object, ...) {
    out <- list()
    class(out) <- "summary.rSPDEobj"
    out$type <- object$type
    if (out$type == "Matern approximation") {
        out$kappa <- object$kappa
        out$sigma <- object$sigma
        out$tau <- object$tau
        out[["range"]] <- object[["range"]]
        out$nu <- object$nu
    }
    out$m <- object$m
    out$stationary <- object$stationary
    out$parameterization <- object$parameterization
    out$n <- dim(object$L)[1]
    return(out)
}

#' @param x an object of class "summary.rSPDEobj", usually, a result of a call
#'   to [summary.rSPDEobj()].
#' @export
#' @method print summary.rSPDEobj
#' @rdname summary.rSPDEobj
print.summary.rSPDEobj <- function(x, ...) {
    cat("Type of approximation: ", x$type, "\n")
    cat("Parameterization: ", x$parameterization, "\n")
    if (x$type == "Matern approximation") {
        if (x$parameterization == "spde") {
            cat(
                "Parameters of covariance function: kappa = ",
                x$kappa, ", tau = ", x$tau, ", nu = ", x$nu, "\n"
            )
        } else {
            cat(
                "Parameters of covariance function: range = ",
                x[["range"]], ", sigma = ", x$sigma, ", nu = ", x$nu, "\n"
            )
        }
    }
    cat("Order or rational approximation: ", x$m, "\n")
    cat("Size of discrete operators: ", x$n, " x ", x$n, "\n")
    if (x$stationary) {
        cat("Stationary Model\n")
    } else {
        cat("Non-Stationary Model")
    }
}

#' @export
#' @method print rSPDEobj
#' @rdname summary.rSPDEobj
print.rSPDEobj <- function(x, ...) {
    print.summary.rSPDEobj(summary(x))
}


#' Summarise CBrSPDE objects
#'
#' Summary method for class "CBrSPDEobj"
#'
#' @param object an object of class "CBrSPDEobj", usually, a result of a call
#'   to [matern.operators()].
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary CBrSPDEobj
#' @examples
#' # Compute the covariance-based rational approximation of a
#' # Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#' range <- sqrt(8 * nu) / kappa
#'
#' # create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' # compute rational approximation of covariance function at 0.5
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
#'   (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
#' op_cov <- matern.operators(
#'   loc_mesh = x, nu = nu,
#'   range = range, sigma = sigma, d = 1, m = 2,
#'   parameterization = "matern"
#' )
#'
#' op_cov
summary.CBrSPDEobj <- function(object, ...) {
    out <- list()
    class(out) <- "summary.CBrSPDEobj"
    out$type <- object$type
    out$kappa <- object$kappa
    out$sigma <- object$sigma
    out$theta <- object$theta
    out$tau <- object$tau
    out[["range"]] <- object[["range"]]
    out$nu <- object$nu
    out$m <- object$m
    out$stationary <- object$stationary
    out$parameterization <- object$parameterization
    out$n <- dim(object$C)[1]
    out[["type_rational_approximation"]] <-
        object[["type_rational_approximation"]]
    return(out)
}

#' @param x an object of class "summary.CBrSPDEobj", usually, a result of a call
#'   to [summary.CBrSPDEobj()].
#' @export
#' @method print summary.CBrSPDEobj
#' @rdname summary.CBrSPDEobj
print.summary.CBrSPDEobj <- function(x, ...) {
    cat("Type of approximation: ", x$type, "\n")
    cat("Parameterization: ", x$parameterization, "\n")
    cat(
        "Type of rational approximation: ",
        x[["type_rational_approximation"]], "\n"
    )
    if (x$stationary) {
        if (x$parameterization == "spde") {
            cat(
                "Parameters of covariance function: kappa = ",
                x$kappa, ", tau = ", x$tau, ", nu = ", x$nu, "\n"
            )
        } else {
            cat(
                "Parameters of covariance function: range = ",
                x[["range"]], ", sigma = ", x$sigma, ", nu = ", x$nu, "\n"
            )
        }
    } else if (!is.null(x$theta)) {
        cat(
            "Parameters of covariance function: theta = ",
            x$theta, ", nu = ", x$nu, "\n"
        )
    }
    
    cat("Order or rational approximation: ", x$m, "\n")
    cat("Size of discrete operators: ", x$n, " x ", x$n, "\n")
    if (x$stationary) {
        cat("Stationary Model\n")
    } else {
        cat("Non-Stationary Model")
    }
}

#' @export
#' @method print CBrSPDEobj
#' @rdname summary.CBrSPDEobj
print.CBrSPDEobj <- function(x, ...) {
    print.summary.CBrSPDEobj(summary(x))
}



#' Summarise CBrSPDEobj2d objects
#'
#' Summary method for class "CBrSPDEobj2d"
#'
#' @param object an object of class "CBrSPDEobj2d", usually, a result of a call
#'   to [matern2d.operators()].
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary CBrSPDEobj2d
#' @examples
#' library(fmesher)
#' n_loc <- 2000
#' loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
#' mesh_2d <- fm_mesh_2d(loc = loc_2d_mesh, cutoff = 0.03, max.edge = c(0.1, 0.5))
#' op <- matern2d.operators(mesh = mesh_2d)
#' op
summary.CBrSPDEobj2d <- function(object, ...) {
    out <- list()
    class(out) <- "summary.CBrSPDEobj2d"
    out$type <- object$type
    out$hx <- object$hx
    out$hy <- object$hy
    out$hxy <- object$hxy
    out$sigma <- object$sigma
    out$nu <- object$nu
    out$m <- object$m
    out$stationary <- object$stationary
    out$n <- dim(object$C)[1]
    out[["type_rational_approximation"]] <-
        object[["type_rational_approximation"]]
    return(out)
}

#' @param x an object of class "summary.CBrSPDEobj2d", usually, a result of a call
#'   to [summary.CBrSPDEobj2d()].
#' @export
#' @method print summary.CBrSPDEobj2d
#' @rdname summary.CBrSPDEobj2d
print.summary.CBrSPDEobj2d <- function(x, ...) {
    cat("Type of approximation: ", x$type, "\n")
    cat(
        "Type of rational approximation: ",
        x[["type_rational_approximation"]], "\n"
    )

    cat("Parameters of covariance function: sigma = ",
        x$sigma, ", hx = ", x$hx, ", hy = ", x$hy, 
        ", hxy = ", x$hxy, ", nu = ", x$nu, "\n"
            )
        
    cat("Order or rational approximation: ", x$m, "\n")
    cat("Size of discrete operators: ", x$n, " x ", x$n, "\n")
}

#' @export
#' @method print CBrSPDEobj2d
#' @rdname summary.CBrSPDEobj2d
print.CBrSPDEobj2d <- function(x, ...) {
    print.summary.CBrSPDEobj2d(summary(x))
}


#' @name create_summary_from_density
#' @title Creates a summary from a density data frame
#' @description Auxiliar function to create summaries from density data drames
#' @param density_df A density data frame
#' @param name Name of the parameter
#' @return A data frame containing a basic summary
#' @noRd

create_summary_from_density <- function(density_df, name) {
    min_x <- min(density_df[, "x"])
    max_x <- max(density_df[, "x"])
    denstemp <- function(x) {
        dens <- sapply(x, function(z) {
            if (z < min_x) {
                return(0)
            } else if (z > max_x) {
                return(0)
            } else {
                return(approx(x = density_df[, "x"], y = density_df[, "y"], xout = z)$y)
            }
        })
        return(dens)
    }
    
    ptemp <- function(q) {
        prob_temp <- sapply(q, function(v) {
            if (v <= min_x) {
                return(0)
            } else if (v >= max_x) {
                return(1)
            } else {
                stats::integrate(
                    f = denstemp, lower = min_x, upper = v,
                    stop.on.error = FALSE
                )$value
            }
        })
        return(prob_temp)
    }
    
    mean_temp <- stats::integrate(
        f = function(z) {
            denstemp(z) * z
        }, lower = min_x, upper = max_x,
        subdivisions = nrow(density_df),
        stop.on.error = FALSE
    )$value
    
    sd_temp <- sqrt(stats::integrate(
        f = function(z) {
            denstemp(z) * (z - mean_temp)^2
        }, lower = min_x, upper = max_x,
        stop.on.error = FALSE
    )$value)
    
    mode_temp <- density_df[which.max(density_df[, "y"]), "x"]
    
    qtemp <- function(p) {
        quant_temp <- sapply(p, function(x) {
            if (x < 0 | x > 1) {
                return(NaN)
            } else {
                return(stats::uniroot(function(y) {
                    ptemp(y) - x
                }, lower = min_x, upper = max_x)$root)
            }
        })
        return(quant_temp)
    }
    
    out <- data.frame(
        mean = mean_temp, sd = sd_temp, `0.025quant` = qtemp(0.025),
        `0.5quant` = qtemp(0.5), `0.975quant` = qtemp(0.975), mode = mode_temp
    )
    rownames(out) <- name
    colnames(out) <- c(
        "mean", "sd", "0.025quant",
        "0.5quant", "0.975quant", "mode"
    )
    return(out)
}


#' Summarise spacetime objects
#'
#' Summary method for class "spacetimeobj"
#'
#' @param object an object of class "spacetimeobj", usually, a result of a call
#'   to [spacetime.operators()].
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary spacetimeobj
summary.spacetimeobj <- function(object, ...) {
    out <- list()
    class(out) <- "summary.spacetimeobj"
    
    out$kappa <- object$kappa
    out$sigma <- object$sigma
    out$gamma <- object$gamma
    out$alpha <- object$alpha
    out$beta <- object$beta
    out$rho <- object$rho
    if(object$has_graph) {
        out$domain <- "graph"
    } else {
        if(object$d == 1) {
            out$domain <- "interval"
        } else {
            out$domain <- "2d region"
        }
    }
    return(out)
}

#' @param x an object of class "summary.spacetimeobj", usually, a result of a call
#'   to [summary.spacetimeobj()].
#' @export
#' @method print summary.spacetimeobj
#' @rdname summary.spacetimeobj
print.summary.spacetimeobj <- function(x, ...) {
    
    
    cat("Model parameters: kappa = ", x$kappa, ", sigma = ", x$sigma, 
        ", gamma = ", x$gamma, ", rho = ", x$rho, ", alpha = ", x$alpha, 
        ", beta = ", x$beta,"\n")
        
    cat("Type of spatial domain: ", x$domain, "\n")
}

#' @export
#' @method print spacetimeobj
#' @rdname summary.spacetimeobj
print.spacetimeobj <- function(x, ...) {
    print.summary.spacetimeobj(summary(x))
}

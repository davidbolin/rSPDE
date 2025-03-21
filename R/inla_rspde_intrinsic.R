
#' Rational approximations of fractional intrinsic fields
#'
#' `rspde.intrinsic` computes a Finite Element Method (FEM) approximation of a
#' Gaussian random field defined as the solution to the stochastic partial
#' differential equation (SPDE):
#' \deqn{(-\Delta)^{(\nu+d/2)/2}\tau u = W}.
#'
#' @param mesh Spatial mesh for the FEM approximation.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not
#' be estimated. If nu is `NULL`, it will be estimated.
#' @param nu.upper.bound Upper bound for the smoothness parameter \eqn{\nu}. If `NULL`, it will be set to 2.
#' @param mean.correction Add mean correction for extreme value models?
#' @param rspde.order The order of the covariance-based rational SPDE approach. The default order is 1.
#' @param prior.tau A list specifying the prior for the variance parameter \eqn{\tau}.
#' This list may contain two elements: `mean` and/or `precision`, both of which must
#' be numeric scalars. 
#' @param prior.nu a list containing the elements `mean` and `prec`
#' for beta distribution, or `loglocation` and `logscale` for a
#' truncated lognormal distribution. `loglocation` stands for
#' the location parameter of the truncated lognormal distribution in the log
#' scale. `prec` stands for the precision of a beta distribution.
#' `logscale` stands for the scale of the truncated lognormal
#' distribution on the log scale. Check details below.
#' @param prior.nu.dist The distribution of the smoothness parameter.
#' The current options are "beta" or "lognormal". The default is "lognormal".
#' @param nu.prec.inc Amount to increase the precision in the beta prior
#' distribution. Check details below.
#' @param diagonal Number added to diagonal of Q for increased stability. 
#' @param type.rational.approx Which type of rational approximation
#' should be used? The current types are "brasil", "chebfun" or "chebfunLB".
#' @param shared_lib String specifying which shared library to use for the Cgeneric
#' implementation. Options are "detect", "INLA", or "rSPDE". You may also specify the
#' direct path to a .so (or .dll) file.
#' @param debug Logical value indicating whether to enable INLA debug mode.
#' @param cache Use caching internally in the estimation?
#' @param opts A list of options passed to RSpectra::eigs function. 
#' See RSpectra documentation for available options.
#' @param scaling A positive numeric value of length 1 for scaling the model.
#'   If NULL (default), it will be computed using RSpectra::eigs.
#'   Must be positive if provided.
#' @param ... Additional arguments passed internally for configuration purposes.
#' @return An object of class `inla_rspde_intrinsic` representing the FEM approximation of
#' the intrinsic Gaussian random field.
#' @export
rspde.intrinsic <- function(mesh,
                            nu = NULL,
                            nu.upper.bound = 2,
                            mean.correction = FALSE,
                            rspde.order = 1,
                            prior.tau = NULL,
                            prior.nu = NULL,
                            prior.nu.dist = "lognormal",
                            nu.prec.inc = 0.01,
                            diagonal = 1e-5,
                            type.rational.approx = "brasil",
                            shared_lib = "detect",
                            debug = FALSE,
                            cache = TRUE,
                            scaling = NULL,
                            opts = NULL,
                            ...) {
    # Validate mesh input
    if (inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d"))) {
        d <- fmesher::fm_manifold_dim(mesh)
    } else if (!is.null(mesh$d)) {
        d <- mesh$d
    } else {
        stop("The mesh object should either be an INLA mesh object or contain d, the dimension!")
    }
    fem_mesh <- fm_fem(mesh)
    G <- fem_mesh$g1
    C <- fem_mesh$c0
    Ci <- Diagonal(dim(C)[1],1/diag(C))
    
    if (nu.upper.bound - floor(nu.upper.bound) == 0) {
        nu.upper.bound <- nu.upper.bound - 1e-5
    }
    
    if(!is.null(nu)){
        nu.upper.bound <- nu
    }
    
    prior.tau <- set_prior(prior.tau, 0, 0.1, p = 1)
    
    est_nu <- 0L
    
    if(is.null(nu)){
        est_nu <- 1L
        nu <- -1.0
    }
    
    result_nu <- handle_prior_nu(prior.nu, nu.upper.bound = nu.upper.bound, nu.prec.inc = nu.prec.inc, prior.nu.dist = prior.nu.dist)
    
    prior.nu <- result_nu$prior.nu
    start.nu <- result_nu$start.nu
    
    rational_table <- as.matrix(get_rational_coefficients(rspde.order, type.rational.approx))
    
    rspde_lib <- get_shared_library(shared_lib)
    
    #Q <- G
    #alpha_ub <- nu.upper.bound + d/2
    #if(alpha>1){
    #    for(i in 1:ceiling(alpha)){
    #        Q <- G%*%G
    #    }
    #}
    op <- matern.operators(kappa = 1, tau = 1, alpha = nu.upper.bound + d/2, G = G, C = C, d = d,
                           mesh = mesh, m = rspde.order, type_rational_approximation = type.rational.approx)
    
    to.inla.matrix <- function(M) {
        out  <-  as(as(as(M, "dMatrix"), "generalMatrix"), "TsparseMatrix")
        ii <- out@i
        out@i <- out@j
        out@j <- ii
        idx <- which(out@i <= out@j)
        out@i <- out@i[idx]
        out@j <- out@j[idx]
        out@x <- out@x[idx]  
        return(out)
    }

    # Add scaling parameter validation
    if (!is.null(scaling)) {
        if (!is.numeric(scaling) || length(scaling) != 1 || scaling <= 0) {
            stop("scaling must be a positive numeric value of length 1")
        }
    } else {
        Cmatrix <- to.inla.matrix(op$Q)
        D <- Diagonal(dim(op$Q)[1], diagonal)
        Cd <- Diagonal(dim(C)[1], 1/sqrt(diag(C)))
        Gg <- Cd%*%G%*%Cd        
        # Use opts argument with RSpectra::eigs
        if(is.null(opts)) { 
            opts = list(tol = 1e-10, maxitr = 1e4)
        }
        scaling <- RSpectra::eigs_sym(as(Gg, "CsparseMatrix"), 2, which = "SM", 
                                      opts = opts)$values[1]
        if(is.na(scaling)){
            stop("Computation of scaling failed, provide the scaling manually or change opts to allow for higher maxitr or lower tol")
        }
    }    
    
    list_args <- 
        list(
            model = "inla_cgeneric_rspde_fintrinsic_model",
            shlib = rspde_lib,
            n = as.integer(nrow(op$Q)),
            est_nu = as.integer(est_nu),
            rspde_order = as.integer(rspde.order),
            dim = as.integer(d),
            mean_correction = as.integer(mean.correction),
            use_cache = as.integer(cache),
            prior.tau.mean = prior.tau$mean,
            prior.tau.precision = prior.tau$precision,
            start_nu = start.nu,
            nu = nu,
            prior.nu.loglocation = prior.nu$loglocation,
            prior.nu.mean = prior.nu$mean,
            prior.nu.prec = prior.nu$prec,
            prior.nu.logscale = prior.nu$logscale,
            nu_upper_bound = nu.upper.bound,
            scaling = scaling,
            prior.nu.dist = prior.nu.dist,
            Q = Cmatrix,
            C = C,
            Ci = Ci,
            G = G,
            D = D,
            rational_table = rational_table
        )
    
    model <- do.call(INLA::inla.cgeneric.define, list_args)
    
    rspde_check_cgeneric_symbol(model)
    
    model$prior.tau <- prior.tau
    model$mesh <- mesh
    model$rspde.order <- rspde.order
    model$type_rational_approximation <- type.rational.approx
    model$est_nu <- est_nu
    model$nu <- nu
    model$nu_upper_bound <- nu.upper.bound
    model$rspde_version <- as.character(packageVersion("rSPDE"))
    
    ### The following objects are provided for backward compatibility
    if(!est_nu){
        model$integer.nu <- (nu %% 1) == 0
    } else{
        model$integer.nu <- FALSE
    }
    model$n.spde <- mesh$n
    ### 
    
    class(model) <- c("inla_rspde_fintrinsic", class(model))
    
    return(model)
}


#' @name rspde.matern.intrinsic
#' @title Intrinsic Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary intrinsic Matern model.
#' Currently, alpha is fixed to 2 and beta is fixed to 1.
#' @param mesh The mesh to build the model. It can be an `inla.mesh` or
#' an `inla.mesh.1d` object. Otherwise, should be a list containing elements d, the dimension, C, the mass matrix,
#' and G, the stiffness matrix.
#' @param alpha Smoothness parameter, need to be 1 or 2.
#' @param mean.correction Add mean correction for extreme value models?
#' @param start.lkappa Starting value for log of kappa.
#' @param prior.lkappa.mean Prior on log kappa to be used for the priors and for the starting values.
#' @param prior.ltau.mean Prior on log tau to be used for the priors and for the starting values.
#' @param prior.lkappa.prec Precision to be used on the prior on log kappa to be used for the priors and for the starting values.
#' @param prior.ltau.prec Precision to be used on the prior on log tau to be used for the priors and for the starting values.
#' @param start.ltau Starting value for log of tau.
#' @param true.scaling Compute the true normalizing constant manually? Default `TRUE`.
#' The alternative is to set this to `FALSE` and set the `diagonal` argument to some small
#' positive value. In the latter case, the model is approximated by a non-intrinsic model
#' with a precision matrix that has the `diagonal` value added to the diagonal.
#' @param diagonal Value of diagonal correction for INLA stability. Default 0.
#' @param debug INLA debug argument
#' @param shared_lib Which shared lib to use for the cgeneric implementation?
#' If "detect", it will check if the shared lib exists locally, in which case it will
#' use it. Otherwise it will use INLA's shared library.
#' If "INLA", it will use the shared lib from INLA's installation. If 'rSPDE', then
#' it will use the local installation (does not work if your installation is from CRAN).
#' Otherwise, you can directly supply the path of the .so (or .dll) file.
#' @param ... Only being used internally.
#'
#' @return An INLA model.
#' @export

rspde.intrinsic.matern <- function(mesh,
                                   alpha = 2,
                                   mean.correction = FALSE,
                                   prior.lkappa.mean = NULL,
                                   prior.ltau.mean = 1,
                                   prior.lkappa.prec = 0.1,
                                   prior.ltau.prec = 0.1,
                                   start.ltau = NULL,
                                   start.lkappa = NULL,
                                   true.scaling = TRUE,
                                   diagonal = 0,
                                   debug = FALSE,
                                   shared_lib = "detect",
                                   ...) {

    cache <- TRUE
    eigen.version <- FALSE
    if(mean.correction || true.scaling) {
        eigen.version = TRUE
    }
    if(diagonal <0){
        stop("diagonal correction needs to be non-negative.")
    }
    if(true.scaling && diagonal > 0) {
        warning("If the true intrinsic scaling is used, there is no need to use the diagonal correction. Consider setting diagonal = 0.")
    }
    if (!(alpha %in% c(1,2))){
        stop("Only alpha = 1 or alpha = 2 implemented.")
    }

    if (prior.lkappa.prec < 0 || prior.ltau.prec < 0) {
        stop("Need positive precisions for the priors.")
    }
    ### Location of object files

    rspde_lib <- get_shared_library(shared_lib)

    if (inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d"))) {
        d <- fmesher::fm_manifold_dim(mesh)
    } else if (!is.null(mesh$d)) {
        d <- mesh$d
    } else {
        stop("The mesh object should either be an INLA mesh object or contain d, the dimension!")
    }

    ### Priors and starting values
    if(is.null(prior.lkappa.mean)){
        mesh.range <- ifelse(d == 2, (max(c(diff(range(mesh$loc[,1])),
                                            diff(range(mesh$loc[, 2])),
                                            diff(range(mesh$loc[,3]))))),
                             diff(mesh$interval))
        prior.lkappa.mean <- log(sqrt(8)/(0.2*mesh.range))
    }

    theta.prior.mean <- c(prior.ltau.mean, prior.lkappa.mean)

    theta.prior.prec <- diag(c(prior.ltau.prec, prior.lkappa.prec))
    start.theta <- theta.prior.mean


    if (!is.null(start.lkappa)) {
        start.theta[2] <- start.lkappa
    }
    if (!is.null(start.ltau)) {
        start.theta[1] <- start.ltau
    }


    ### FEM matrices
    if (inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d"))) {

        if (d == 1) {
            fem_mesh <- fem_mesh_order_1d(mesh, m_order = alpha + 1)
        } else {
            fem_mesh <- fm_fem(mesh, order = alpha + 1)
        }

    } else {
        if (is.null(mesh$C) || is.null(mesh$G)) {
            stop("If mesh is not an fm_mesh_1d/2d object, you should manually supply a list with elements c0, g1, g2...")
        }
        fem_mesh <- generic_fem_mesh_order(mesh, m_order = alpha + 1)
    }

    n_cgeneric <- ncol(fem_mesh[["c0"]])

    fem_mesh_orig <- fem_mesh

    if(eigen.version) {
        C <- fem_mesh[["c0"]]
        G <- fem_mesh[["g1"]]
        fem_mesh <- fem_mesh[setdiff(names(fem_mesh), c("ta", "va"))]

        fem_mesh <- lapply(fem_mesh, transpose_cgeneric)

        if(alpha==1) {
            graph_opt <- fem_mesh[["g2"]]
        } else {
            graph_opt <- fem_mesh[["g3"]]
        }
        if(cache) {
            model <- do.call(
                eval(parse(text = "INLA::inla.cgeneric.define")),
                list(
                    model = "inla_cgeneric_rspde_intrinsic_eigen_cache",
                    shlib = rspde_lib,
                    n = as.integer(n_cgeneric),
                    debug = debug,
                    graph_opt_i = graph_opt@i,
                    graph_opt_j = graph_opt@j,
                    C = C,
                    G = G,
                    theta.prior.mean = theta.prior.mean,
                    theta.prior.prec = theta.prior.prec,
                    start.theta = start.theta,
                    alpha = as.integer(alpha),
                    mean_correction = as.integer(mean.correction),
                    true_scaling = as.integer(true.scaling)
                )
            )
        } else {
            model <- do.call(
                eval(parse(text = "INLA::inla.cgeneric.define")),
                list(
                    model = "inla_cgeneric_rspde_intrinsic_eigen",
                    shlib = rspde_lib,
                    n = as.integer(n_cgeneric),
                    debug = debug,
                    graph_opt_i = graph_opt@i,
                    graph_opt_j = graph_opt@j,
                    C = C,
                    G = G,
                    theta.prior.mean = theta.prior.mean,
                    theta.prior.prec = theta.prior.prec,
                    start.theta = start.theta,
                    alpha = as.integer(alpha),
                    mean_correction = as.integer(mean.correction),
                    true_scaling = as.integer(true.scaling)
                )
            )
        }

    } else {
        fem_mesh <- fem_mesh[setdiff(names(fem_mesh), c("ta", "va"))]

        fem_mesh <- lapply(fem_mesh, transpose_cgeneric)

        C_list <- symmetric_part_matrix(fem_mesh$c0)
        G_1_list <- symmetric_part_matrix(fem_mesh$g1)
        G_2_list <- symmetric_part_matrix(fem_mesh$g2)


        idx_matrices <- list()

        if(alpha == 2) {
            G_3_list <- symmetric_part_matrix(fem_mesh$g3)

            idx_matrices[[1]] <- G_1_list$idx
            idx_matrices[[2]] <- G_2_list$idx
            idx_matrices[[3]] <- G_3_list$idx

            positions_matrices_less <- list()
            positions_matrices_less[[1]] <- match(G_1_list$M, G_3_list$M)
            positions_matrices_less[[2]] <- match(G_2_list$M, G_3_list$M)

            n_tmp <- length(fem_mesh[["g3"]]@x[idx_matrices[[3]]])
            tmp <- rep(0, n_tmp)

            tmp[positions_matrices_less[[1]]] <- fem_mesh$g1@x[idx_matrices[[1]]]
            matrices_less <- tmp

            tmp <- rep(0, n_tmp)
            tmp[positions_matrices_less[[2]]] <- fem_mesh$g2@x[idx_matrices[[2]]]
            matrices_less <- c(matrices_less, tmp)

            tmp <- fem_mesh[["g3"]]@x[idx_matrices[[3]]]
            matrices_less <- c(matrices_less, tmp)

            graph_opt <- fem_mesh[["g3"]]
        } else if (alpha == 1) {
            idx_matrices[[1]] <- G_1_list$idx
            idx_matrices[[2]] <- G_2_list$idx

            positions_matrices_less <- list()
            positions_matrices_less[[1]] <- match(G_1_list$M, G_2_list$M)

            n_tmp <- length(fem_mesh[["g2"]]@x[idx_matrices[[2]]])
            tmp <- rep(0, n_tmp)

            tmp[positions_matrices_less[[1]]] <- fem_mesh$g1@x[idx_matrices[[1]]]
            matrices_less <- tmp

            tmp <- fem_mesh[["g2"]]@x[idx_matrices[[2]]]
            matrices_less <- c(matrices_less, tmp)

            graph_opt <- fem_mesh[["g2"]]
        }

        model <- do.call(
            eval(parse(text = "INLA::inla.cgeneric.define")),
            list(
                model = "inla_cgeneric_rspde_intrinsic_int_model",
                shlib = rspde_lib,
                n = as.integer(n_cgeneric),
                debug = debug,
                matrices_less = as.double(matrices_less),
                graph_opt_i = graph_opt@i,
                graph_opt_j = graph_opt@j,
                theta.prior.mean = theta.prior.mean,
                theta.prior.prec = theta.prior.prec,
                start.theta = start.theta,
                d = as.integer(d),
                alpha = as.integer(alpha)
            )
        )
    }

    model$cgeneric_type <- "intrinsic_matern"
    model$f$diagonal <- diagonal
    model$theta.prior.mean <- theta.prior.mean
    model$theta.prior.prec <- theta.prior.prec
    model$start.theta <- start.theta

    rspde_check_cgeneric_symbol(model)

    class(model) <- c("inla_rspde", "intrinsic", class(model))
    model$parameterization <- "spde"
    model$stationary <- TRUE
    model$est_nu <- FALSE
    model$dim <- d
    model$n.spde <- mesh$n
    model$debug <- debug
    model$mesh <- mesh
    model$alpha <- alpha
    model$fem_mesh <- fem_mesh_orig
    model$rspde_version <- as.character(packageVersion("rSPDE"))
    return(model)
}


#' result summary for intrinsic models
#' @noRd
rspde.intrinsic.result <- function(inla, name, rspde, 
                                    compute.summary = TRUE, 
                                    n_samples = 5000, 
                                    n_density = 1024) {

    nu.upper.bound <- rspde$nu_upper_bound
    result <- list()
    
    if (!rspde$est_nu) {
        row_names <- c("tau")
    } else {
        row_names <- c("tau", "nu")
    }
        
    result$summary.values <- inla$summary.random[[name]]
        
    if (!is.null(inla$marginals.random[[name]])) {
        result$marginals.values <- inla$marginals.random[[name]]
    }
        
    name_theta1 <- "tau"
    name_theta2 <- "nu"

    name_theta1_model <- "tau"
    name_theta2_model <- "nu"

    result[[paste0("summary.log.", name_theta1_model)]] <- INLA::inla.extract.el(
            inla$summary.hyperpar,
            paste("Theta1 for ", name, "$", sep = "")
    )
    rownames(result[[paste0("summary.log.", name_theta1_model)]]) <- paste0("log(", name_theta1_model, ")")
        
    if (rspde$est_nu) {
        result$summary.logit.nu <- INLA::inla.extract.el(inla$summary.hyperpar,
                                                         paste("Theta2 for ", name, "$", sep = ""))
        rownames(result$summary.logit.nu) <- "logit(nu)"
    }
        
    if (!is.null(inla$marginals.hyperpar[[paste0("Theta1 for ", name)]])) {
        result[[paste0("marginals.log.", name_theta1_model)]] <- INLA::inla.extract.el(
            inla$marginals.hyperpar,
            paste("Theta1 for ", name, "$", sep = "")
        )
        names(result[[paste0("marginals.log.", name_theta1_model)]]) <- name_theta1_model
        
        if (rspde$est_nu) {
            result$marginals.logit.nu <- INLA::inla.extract.el(
                inla$marginals.hyperpar,
                paste("Theta2 for ", name, "$", sep = "")
            )
            names(result$marginals.logit.nu) <- "nu"
        }
            
            
        result[[paste0("marginals.", name_theta1)]] <- lapply(
            result[[paste0("marginals.log.", name_theta1)]],
            function(x) {
                    INLA::inla.tmarginal(function(y) exp(y), x)
                }
            )
    
        if (rspde$est_nu) {
            result$marginals.nu <- lapply(result$marginals.logit.nu,
                function(x) {
                    INLA::inla.tmarginal(
                        function(y) {
                            nu.upper.bound * exp(y) / (1 + exp(y))
                        }, x)
                })
        }
        
        if (compute.summary) {
            norm_const <- function(density_df) {
                min_x <- min(density_df[, "x"])
                max_x <- max(density_df[, "x"])
                denstemp <- function(x) {
                    dens <- sapply(x, function(z) {
                        if (z < min_x) {
                            return(0)
                        } else if (z > max_x) {
                            return(0)
                        } else {
                            return(approx(
                                x = density_df[, "x"],
                                y = density_df[, "y"], xout = z
                            )$y)
                        }
                    })
                    return(dens)
                }
                norm_const <- stats::integrate(
                    f = function(z) {
                        denstemp(z)
                    }, lower = min_x, upper = max_x,
                    subdivisions = nrow(density_df),
                    stop.on.error = FALSE
                )$value
                return(norm_const)
            }
            
            norm_const_theta1 <- norm_const(result[[paste0("marginals.", name_theta1)]][[name_theta1]])
            result[[paste0("marginals.", name_theta1)]][[name_theta1]][, "y"] <-
                result[[paste0("marginals.", name_theta1)]][[name_theta1]][, "y"] / norm_const_theta1
            
            result[[paste0("summary.", name_theta1)]] <- create_summary_from_density(result[[paste0("marginals.", name_theta1)]][[name_theta1]],
                                                                                     name = name_theta1
            )
            
            
            if (rspde$est_nu) {
                norm_const_nu <- norm_const(result$marginals.nu$nu)
                result$marginals.nu$nu[, "y"] <-
                    result$marginals.nu$nu[, "y"] / norm_const_nu
                
                result$summary.nu <- create_summary_from_density(result$marginals.nu$nu,
                                                                 name = "nu"
                )
            }
        }
    }
    
    result$n_par <- length(rspde$start.theta)
    
    class(result) <- c("rspde_result", "rspde_intrinsic")
    result$stationary <- TRUE
    result$parameterization <- "spde"
    
    result$params <- c(name_theta1)
    if (rspde$est_nu) {
        result$params <- c(result$params, "nu")
    }

    if (!is.null(result$summary.nu)) {
        if(nu.upper.bound - result$summary.nu$mean < 0.1 || nu.upper.bound - result$summary.nu$mode < 0.1){
            warning("the mean or mode of nu is very close to nu.upper.bound, please consider increasing nu.upper.bound, and refitting the model.")
        }
    }
    
    return(result)
}



#'
#' @title rSPDE inlabru mapper
#' @name bru_get_mapper.inla_rspde_fintrinsic
#' @param model An `inla_rspde_fintrinsic` object for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_get_mapper.inla_rspde_fintrinsic
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_get_mapper, inla_rspde_fintrinsic)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_rspde_fintrinsic)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_rspde_fintrinsic)
#'   S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde_fintrinsic)
#' }
bru_get_mapper.inla_rspde_fintrinsic <- function(model, ...) {
    stopifnot(requireNamespace("inlabru"))
    inlabru_version <- as.character(packageVersion("inlabru"))
    if(inlabru_version >= "2.11.1.9022"){
        n_rep <- model[["rspde.order"]] + 1
        if((model[["est_nu"]] == 0L) && (model[["integer.nu"]])){
            n_rep <- 1
        }
        inlabru::bru_mapper_repeat(inlabru::bru_mapper(model[["mesh"]]), n_rep = n_rep)
    } else{
        mapper <- list(model = model)
        inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_inla_rspde_fintrinsic")
    }
}

#' @param mapper A `bru_mapper_inla_rspde` object
#' @rdname bru_get_mapper.inla_rspde
ibm_n.bru_mapper_inla_rspde_fintrinsic <- function(mapper, ...) {
    model <- mapper[["model"]]
    integer_nu <- model$integer.nu
    rspde_order <- model$rspde.order
    if (integer_nu) {
        factor_rspde <- 1
    } else {
        factor_rspde <- rspde_order + 1
    }
    factor_rspde * model$n.spde
}
#' @rdname bru_get_mapper.inla_rspde
ibm_values.bru_mapper_inla_rspde_fintrinsic <- function(mapper, ...) {
    seq_len(inlabru::ibm_n(mapper))
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_get_mapper.inla_rspde
ibm_jacobian.bru_mapper_inla_rspde_fintrinsic <- function(mapper, input, ...) {
    model <- mapper[["model"]]
    if (!is.null(input) && !is.matrix(input) && !inherits(input, "Spatial")) {
        input <- as.matrix(input)
    }
    
    if (model$est_nu) {
        nu <- NULL
    } else {
        nu <- model$nu
    }
    
    rspde_order <- model$rspde.order
    rSPDE::rspde.make.A(
        mesh = model$mesh, loc = input,
        rspde.order = rspde_order,
        nu = nu
    )
}

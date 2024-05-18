#' @name rspde.matern.intrinsic
#' @title Intrinsic Matern rSPDE model object for INLA
#' @description Creates an INLA object for a stationary intrinsic Matern model.
#' Currently, alpha is fixed to 2 and beta is fixed to 1. 
#' @param mesh The mesh to build the model. It can be an `inla.mesh` or
#' an `inla.mesh.1d` object. Otherwise, should be a list containing elements d, the dimension, C, the mass matrix,
#' and G, the stiffness matrix.
#' @param alpha Smoothness parameter, need to be 1 or 2. 
#' @param mean.correction Add mean correction for extreme value models?
#' @param prior.kappa a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale.
#' @param prior.tau a list containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale.
#' @param start.lkappa Starting value for log of kappa.
#' @param prior.kappa.mean Prior kappa to be used for the priors and for the starting values.
#' @param prior.tau.mean Prior tau to be used for the priors and for the starting values.
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
    
    rspde_lib <- shared_lib
    
    if (shared_lib == "INLA") {
        rspde_lib <- INLA::inla.external.lib("rSPDE")
    } else if (shared_lib == "rSPDE") {
        rspde_lib <- system.file("shared", package = "rSPDE")
        if (Sys.info()["sysname"] == "Windows") {
            rspde_lib <- paste0(rspde_lib, "/rspde_cgeneric_models.dll")
        } else {
            rspde_lib <- paste0(rspde_lib, "/rspde_cgeneric_models.so")
        }
    } else if (shared_lib == "detect") {
        rspde_lib_local <- system.file("shared", package = "rSPDE")
        if (Sys.info()["sysname"] == "Windows") {
            rspde_lib_local <- paste0(rspde_lib_local, "/rspde_cgeneric_models.dll")
        } else {
            rspde_lib_local <- paste0(rspde_lib_local, "/rspde_cgeneric_models.so")
        }
        if (file.exists(rspde_lib_local)) {
            rspde_lib <- rspde_lib_local
        } else {
            rspde_lib <- INLA::inla.external.lib("rSPDE")
        }
    }
    
    if (inherits(mesh, c("inla.mesh", "inla.mesh.1d"))) {
        d <- get_inla_mesh_dimension(mesh)
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
    if (inherits(mesh, c("inla.mesh", "inla.mesh.1d"))) {
        
        if (d == 1) {
            fem_mesh <- fem_mesh_order_1d(mesh, m_order = alpha + 1)
        } else {
            fem_mesh <- fm_fem(mesh, order = alpha + 1)
        }
        
    } else {
        if (is.null(mesh$C) || is.null(mesh$G)) {
            stop("If mesh is not an inla.mesh object, you should manually supply a list with elements c0, g1, g2...")
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
    return(model)
}

utils::globalVariables(c("C", "C_inv", "C_inv_G", "G", "d", "loc", "n",
                         "n_m", "nu", "nu_lower_bound", "nu_upper_bound",
                         "do_optimize", "idx_symmetric", "n_Q", "rspde_order",
                         "graph_opt","fem_matrices", "sharp",
                         "prior.kappa", "prior.nu", "prior.tau",
                         "start.lkappa", "start.ltau", "start.lnu"))

#' @importFrom stats dnorm pnorm
#' @importFrom methods as
#' @name 'inla.rgeneric.cov_rspde_general'
#' @title Generic INLA method for the covariance-based rSPDE approach
#' @description Generic INLA method for the covariance-based rSPDE approach
#' @param cmd INLA parameter
#' @param theta Parameters of the model
#' @param ... Additional arguments to the call, or arguments with changed values. Use name = NULL to remove the argument name.
#' @param args Arguments.
#' @return An INLA model
#' @export
'inla.rgeneric.cov_rspde_general' <- function(cmd = c("graph",
                                                      "Q",
                                                      "mu",
                                                      "initial",
                                                      "log.norm.const",
                                                      "log.prior",
                                                      "quit"),
                                              theta = NULL,
                                              args = NULL,
                                              ...) {
  
  initial <- function(n, theta) {
    return(c(start.ltau, start.lkappa, start.lnu))
  }
  
  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      ltau = theta[1L],
      lkappa = theta[2L],
      lnu = theta[3L]
    ))
  }
  
  ######## precision matrix
  Q = function(n, theta) {
    param = interpret.theta(n, theta)
    nu = min(exp(param$lnu) + nu_lower_bound, nu_upper_bound)
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    
    if(do_optimize){
      return(rSPDE::rspde.matern.precision.opt(kappa=kappa, nu=nu, tau=tau, 
                                   rspde_order=rspde_order, 
                                   d=d, fem_matrices = fem_matrices,
                                   sharp=sharp, graph = NULL))
    } else{
      return(rSPDE::rspde.matern.precision(kappa=kappa, nu=nu, tau=tau, 
                               rspde_order=rspde_order, 
                               d=d, fem_mesh_matrices = fem_matrices))
    }
  }
  
  
  ############################# mean
  mu = function(n, theta)
  {
    return(numeric(0))
  }
  ###################### log normal constant
  log.norm.const <- function(n, theta) {
    return(numeric(0))
    
  }
  
  ############################# graph skeleton
  graph <- function(n, theta) {
    return(graph_opt)
  }
  
  
  
  ######################## log prior
  log.prior <- function(n, theta) {
    param = interpret.theta(n, theta)
    
    if (param$lnu > log(nu_upper_bound - nu_lower_bound)) {
      tdnorm_nu = -Inf
    } else{
      tdnorm_nu = dnorm(param$lnu, 0, 1, log = TRUE) -
        pnorm(log(nu_upper_bound - nu_lower_bound), prior.nu$meanlog, 
              prior.nu$sdlog, log.p = TRUE)
    }
    
    res <- tdnorm_nu + dnorm(param$lkappa, prior.kappa$meanlog, 
                             prior.kappa$sdlog, log = TRUE) +
      dnorm(param$ltau, prior.tau$meanlog, prior.tau$sdlog, log = TRUE)
    return(res)
  }
  
  
  quit <- function(n, theta) {
    return(invisible())
  }
  
  if (!length(theta))
    theta = initial(n, theta)
  res <-
    do.call(match.arg(cmd), args = list(n = as.integer(args$n), theta = theta))
  return(res)
  
}

#' @name 'inla.rgeneric.cov_rspde_frac_alpha'
#' @title Non-integer fixed smoothmess generic INLA method 
#' @description Non-integer fixed smoothmess generic INLA method 
#' for the covariance-based rSPDE approach
#' @param cmd INLA parameter
#' @param theta Parameters of the model
#' @param ... Additional arguments to the call, or arguments with changed values. Use name = NULL to remove the argument name.
#' @param args Arguments.
#' @return An INLA model
#' @export
'inla.rgeneric.cov_rspde_frac_alpha' <- function(cmd = c("graph",
                                                         "Q",
                                                         "mu",
                                                         "initial",
                                                         "log.norm.const",
                                                         "log.prior",
                                                         "quit"),
                                                 theta = NULL,
                                                 args = NULL,
                                                 ...) {
  
  initial <- function(n, theta) {
    return(c(start.ltau, start.lkappa))
  }
  
  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      ltau = theta[1L],
      lkappa = theta[2L]
    ))
  }
  
  ######## precision matrix
  Q = function(n, theta) {
    param = interpret.theta(n, theta)
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    
    if(do_optimize){
      return(rSPDE::rspde.matern.precision.opt(kappa=kappa, nu=nu, tau=tau, 
                                   rspde_order=rspde_order, 
                                   d=d, fem_matrices = fem_matrices,
                                   sharp=sharp, graph = NULL))
    } else{
      return(rSPDE::rspde.matern.precision(kappa=kappa, nu=nu, tau=tau, 
                               rspde_order=rspde_order, 
                               d=d, fem_mesh_matrices = fem_matrices))
    }
  }
  
  
  ############################# mean
  mu = function(n, theta)
  {
    return(numeric(0))
  }
  ###################### log normal constant
  log.norm.const <- function(n, theta) {
    return(numeric(0))
    
  }
  
  ############################# graph skeleton
  graph <- function(n, theta) {
    return(graph_opt)
  }
  
  ######################## log prior
  log.prior <- function(n, theta) {
    param = interpret.theta(n, theta)
    
    res <- dnorm(param$lkappa, prior.kappa$meanlog, 
                 prior.kappa$sdlog, log = TRUE) +
      dnorm(param$ltau, prior.tau$meanlog, 
            prior.tau$sdlog, log = TRUE)
    return(res)
  }
  
  
  quit <- function(n, theta) {
    return(invisible())
  }
  
  if (!length(theta))
    theta = initial(n, theta)
  res <-
    do.call(match.arg(cmd), args = list(n = as.integer(args$n), theta = theta))
  return(res)
}


#' @name 'inla.rgeneric.cov_rspde_int_alpha'
#' @title Integer fixed smoothmess generic INLA method 
#' @description Integer fixed smoothmess generic INLA method 
#' for the covariance-based rSPDE approach
#' @param cmd INLA parameter
#' @param theta Parameters of the model
#' @param ... Additional arguments to the call, or arguments with changed values. Use name = NULL to remove the argument name.
#' @param args Arguments.
#' @return An INLA model
#' @export
'inla.rgeneric.cov_rspde_int_alpha' <- function(cmd = c("graph",
                                                        "Q",
                                                        "mu",
                                                        "initial",
                                                        "log.norm.const",
                                                        "log.prior",
                                                        "quit"),
                                                theta = NULL,
                                                args = NULL,
                                                ...) {
  
  initial <- function(n, theta) {
    return(c(start.ltau, start.lkappa))
  }
 
  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      ltau = theta[1L],
      lkappa = theta[2L]
    ))
  }
  
  ######## precision matrix
  Q = function(n, theta) {
    param = interpret.theta(n, theta)
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    

    
    if(do_optimize){
      return(rSPDE::rspde.matern.precision.integer.opt(kappa, nu, tau,
                                       d, fem_matrices, graph = NULL))
    } else{
      return(rSPDE::rspde.matern.precision.integer(kappa=kappa, nu=nu, tau=tau, 
                               d=d, fem_mesh_matrices = fem_matrices))
    }
  }
  
  
  ############################# mean
  mu = function(n, theta)
  {
    return(numeric(0))
  }
  ###################### log normal constant
  log.norm.const <- function(n, theta) {
    return(numeric(0))
    
  }
  
  ############################# graph skeleton
  graph <- function(n, theta) {
    return(graph_opt)
  }
  
  
  
  ######################## log prior
  log.prior <- function(n, theta) {
    param = interpret.theta(n, theta)
    
    res <- dnorm(param$lkappa, prior.kappa$meanlog,
                 prior.kappa$sdlog, log = TRUE) +
      dnorm(param$ltau, prior.tau$meanlog, 
            prior.tau$sdlog, log = TRUE)
    return(res)
  }
  
  
  quit <- function(n, theta) {
    return(invisible())
  }
  
  if (!length(theta))
    theta = initial(n, theta)
  res <-
    do.call(match.arg(cmd), args = list(n = as.integer(args$n), theta = theta))
  return(res)
  
}



#' @name rspde.matern
#' @title Create INLA-based rSPDE models
#' @description Create INLA-based rSPDE models with general smoothness
#' @param mesh An INLA mesh
#' @param nu_lower_bound Lower bound for the smoothness parameter.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not be estimated.
#' If nu is NULL, it will be estimated.
#' @param sharp The sparsity graph should have the correct sparsity (costs more to perform
#' a sparsity analysis) or an upper bound for the sparsity? If TRUE, the graph
#' will have the correct sparsity.
#' @param debug INLA debug argument.
#' @param optimize Should the model be optimized? In this case the sparsities of
#' the matrices will be analyzed.
#' @param prior.kappa a list containing the elements meanlog and sdlog, that is,
#' the mean and standard deviation on the log scale.
#' @param prior.nu a list containing the elements meanlog and sdlog, that is,
#' the mean and standard deviation on the log scale.
#' @param prior.tau a list containing the elements meanlog and sdlog, that is,
#' the mean and standard deviation on the log scale.
#' @param start.lkappa Starting value for log of kappa.
#' @param start.lnu Starting value for log(nu-nu_lower_bound).
#' @param start.ltau Starting value for log of tau.
#' @return An INLA model.
#' @export

rspde.matern <- function(mesh, 
                         nu_lower_bound = get_inla_mesh_dimension(mesh)/4,
                         nu_upper_bound = 4, rspde_order = 2,
                         nu = NULL, sharp = TRUE,
                         debug = FALSE,
                         optimize=TRUE, 
                         prior.kappa = list(meanlog = 0, sdlog = 1),
                         prior.nu = list(meanlog = 0, sdlog = 1),
                         prior.tau = list(meanlog = 0, sdlog = 1),
                         start.lkappa = 0,
                         start.lnu = 0,
                         start.ltau = 0
){
  if(is.null(prior.kappa$meanlog)){
    prior.kappa$meanlog <- 0
  }
  if(is.null(prior.nu$meanlog)){
    prior.nu$meanlog <- 0
  }
  if(is.null(prior.tau$meanlog)){
    prior.tau$meanlog <- 0
  }
  if(is.null(prior.kappa$sdlog)){
    prior.kappa$sdlog <- 1
  }
  if(is.null(prior.nu$sdlog)){
    prior.nu$sdlog <- 1
  }
  if(is.null(prior.tau$sdlog)){
    prior.tau$sdlog <- 1
  }
  
  if(mesh$manifold == "R1"){
    d = 1
  } else if(mesh$manifold == "R2"){
    d = 2
  } else{
    stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
  }
  if(nu_lower_bound<d/4){
    stop("The lower bound for nu cannot be less than d/4")
  }
  
  if(nu_upper_bound-floor(nu_upper_bound)==0){
    nu_upper_bound = nu_upper_bound-1e-5
  }
  
  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  
  if(fixed_nu){
    alpha = nu + d/2
    integer_alpha <- (alpha%%1==0)
  } else{
    integer_alpha = FALSE
  }
  
  if(fixed_nu){
    nu_order <- nu
  } else{
    nu_order <- nu_upper_bound
  }
  
  if(optimize){
    
    beta = nu_order / 2 + d / 4
    
    m_alpha = max(1, floor(2 * beta))
    
    if(integer_alpha){
      fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha)
    } else{
      fem_mesh <- INLA::inla.mesh.fem(mesh, order = m_alpha+1)
    }
  } else{
    beta = nu_order / 2 + d / 4
    
    m_alpha = max(1, floor(2 * beta))
  
    if(integer_alpha){
      fem_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha)
    } else{
      fem_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha+1)
    }
  }
  
  if(optimize){
    if(integer_alpha){
      result_sparsity <- analyze_sparsity_rspde(nu_upper_bound=nu_order, dim=d, 
                                                rspde_order=rspde_order, 
                                                fem_mesh_matrices = fem_mesh,
                                                include_higher_order = FALSE)
    } else{
      if(sharp){
        result_sparsity <- analyze_sparsity_rspde(nu_upper_bound=nu_order, dim=d, 
                                                  rspde_order=rspde_order, 
                                                  fem_mesh_matrices = fem_mesh)
      } else{
        result_sparsity <- analyze_sparsity_rspde(nu_upper_bound=nu_order, dim=d, 
                                                  rspde_order=rspde_order, 
                                                  fem_mesh_matrices = fem_mesh,
                                                  include_lower_order = FALSE)
      }
      positions_matrices <- result_sparsity$positions_matrices
    }

      idx_matrices <- result_sparsity$idx_matrices
      positions_matrices_less <- result_sparsity$positions_matrices_less
  } else{
      positions_matrices <- NULL
      idx_matrices <- NULL
      positions_matrices_less <- NULL
  }
  
  if(optimize){
    
    fem_matrices <- list()
    
    if(sharp || integer_alpha){
      fem_matrices[[paste0("G_",m_alpha,"_less")]] <- fem_mesh[[paste0("g",m_alpha)]]@x[idx_matrices[[m_alpha+1]]]
      
      fem_matrices[["C_less"]] <- rep(0, length(fem_matrices[[paste0("G_",m_alpha,"_less")]]))
      fem_matrices[["C_less"]][positions_matrices_less[[1]]] <- fem_mesh$c0@x[idx_matrices[[1]]]
      
      fem_matrices[["G_less"]] <- rep(0, length(fem_matrices[[paste0("G_",m_alpha,"_less")]]))
      fem_matrices[["G_less"]][positions_matrices_less[[2]]] <- fem_mesh$g1@x[idx_matrices[[2]]]
      
      #The case m_alpha=2 already uses G_2_less defined above
      if(m_alpha > 2){
        for(j in 2:(m_alpha-1)){
          fem_matrices[[paste0("G_",j,"_less")]] <- rep(0, length(fem_matrices[[paste0("G_",m_alpha,"_less")]]))
          fem_matrices[[paste0("G_",j,"_less")]][positions_matrices_less[[j+1]]] <- fem_mesh[[paste0("g",j)]]@x[idx_matrices[[j+1]]]
        }
      }
    }

    
    if(!integer_alpha){
      fem_matrices[[paste0("G_",m_alpha+1)]] <- fem_mesh[[paste0("g",m_alpha+1)]]@x[idx_matrices[[m_alpha+2]]]
      
      fem_matrices[["G"]] <- rep(0, length(fem_matrices[[paste0("G_",m_alpha+1)]]))
      fem_matrices[["G"]][positions_matrices[[2]]] <- fem_mesh$g1@x[idx_matrices[[2]]]
      
      fem_matrices[["C"]] <- rep(0, length(fem_matrices[[paste0("G_",m_alpha+1)]]))
      fem_matrices[["C"]][positions_matrices[[1]]] <- fem_mesh$c0@x[idx_matrices[[1]]]
      if(m_alpha > 1){
        for(j in 2:(m_alpha)){
          fem_matrices[[paste0("G_",j)]] <- rep(0, length(fem_matrices[[paste0("G_",m_alpha+1)]]))
          fem_matrices[[paste0("G_",j)]][positions_matrices[[j+1]]] <- fem_mesh[[paste0("g",j)]]@x[idx_matrices[[j+1]]]
        }
      } 
    }
  }
  
  if(!fixed_nu){
     if(optimize){
       graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_mesh, dim=d,
                                       nu=nu_upper_bound,
                                       rspde_order = rspde_order,
                                       sharp = sharp,
                                       force_non_integer = TRUE)
     } else{
       graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_matrices, dim=d,
                                                    nu=nu_upper_bound,
                                                    rspde_order = rspde_order,
                                                    sharp = TRUE,
                                                    force_non_integer = TRUE)
     }
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_general,
                                        nu_lower_bound=nu_lower_bound,
                                        nu_upper_bound=nu_upper_bound,
                                        fem_matrices=fem_matrices, 
                                        graph_opt=graph_opt,
                                        sharp=sharp,
                                        prior.kappa=prior.kappa,
                                        prior.nu=prior.nu,
                                        prior.tau=prior.tau,
                                        start.lkappa = start.lkappa,
                                        start.lnu = start.lnu,
                                        start.ltau = start.ltau,
                                        d = d, rspde_order = rspde_order, 
                                        n=ncol(C)*(rspde_order+1), 
                                        debug=debug,
                                        do_optimize=optimize, optimize=optimize)
  } else if(!integer_alpha){
    if(optimize){
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_mesh, dim=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   sharp = sharp)
    } else{
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_matrices, dim=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   sharp = TRUE, force_non_integer=TRUE)
    }

    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_frac_alpha,
                                        nu = nu,
                                        fem_matrices=fem_matrices, 
                                        graph_opt=graph_opt,
                                        sharp=sharp,
                                        prior.kappa=prior.kappa,
                                        prior.nu=prior.nu,
                                        prior.tau=prior.tau,
                                        start.lkappa = start.lkappa,
                                        start.ltau = start.ltau,
                                        d = d, rspde_order = rspde_order, 
                                        n=ncol(C)*(rspde_order+1), 
                                        debug=debug,
                                        do_optimize=optimize, optimize=optimize)
  } else{
    if(optimize){
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_mesh, dim=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   sharp = sharp)
    } else{
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_matrices, dim=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   force_non_integer=FALSE)
    }
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_int_alpha,
                                        nu = nu,
                                        fem_matrices=fem_matrices, 
                                        graph_opt=graph_opt,
                                        prior.kappa=prior.kappa,
                                        prior.nu=prior.nu,
                                        prior.tau=prior.tau,
                                        start.lkappa = start.lkappa,
                                        start.ltau = start.ltau,
                                        d = d,
                                        n=ncol(C), 
                                        debug=debug,
                                        do_optimize=optimize, optimize=optimize)
  }
  
  class(model) <- c(class(model), "inla.rspde")
  model$dim = d
  model$est_nu = !fixed_nu
  model$nu_lower_bound = nu_lower_bound
  model$n.spde = mesh$n
  return(model)
}



symmetric_part_matrix <- function(M){
  M  <- as(M, "dgTMatrix")
  idx <- which(M@i <= M@j)
  sM <- cbind(M@i[idx], M@j[idx])
  colnames(sM) <- NULL
  return(list(M = split(sM, seq(nrow(sM))), idx = idx))
}




#' @name rspde.matern.precision.opt
#' @title Optimized precision matrix of the covariance-based rational approximation
#' @description Computes the optimized version of the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param rspde_order The order of the rational approximation
#' @param dim The dimension
#' @param fem_matrices A list containing the FEM-related matrices. The list should contain elements C, G, G_2, G_3, etc.
#' @param graph The sparsity graph of the matrices. If NULL, only a vector
#' of the elements will be returned, if non-NULL, a sparse matrix will be returned.
#' @param sharp The sparsity graph should have the correct sparsity (costs more to perform
#' a sparsity analysis) or an upper bound for the sparsity?
#' @return The precision matrix
#' @export

rspde.matern.precision.opt = function(kappa, nu, tau, rspde_order, dim, fem_matrices, graph=NULL,
                          sharp) {
  
  n_m <- rspde_order
  
  mt <- get(paste0("m", n_m, "t"))
  
  beta = nu / 2 + dim / 4
  
  m_alpha = max(1, floor(2 * beta))
  
  r = sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2*beta))$y
  })
  
  p = sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2*beta))$y
  })
  
  k = approx(mt$alpha, mt$k, cut_decimals(2*beta))$y
  
  if (m_alpha == 1){
    Malpha =  (fem_matrices[["C"]] + fem_matrices[["G"]]/(kappa^2))
  } else if (m_alpha > 1){
    Malpha =  fem_matrices[["C"]] + m_alpha * fem_matrices[["G"]]/(kappa^2)
    for(j in 2:m_alpha){
      Malpha = Malpha +  choose(m_alpha, j) * fem_matrices[[paste0("G_",j)]]/(kappa^(2*j))
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  
  if (m_alpha == 1){
    Malpha2 =  (fem_matrices[["G"]] + fem_matrices[["G_2"]]/(kappa^2))
  } else if (m_alpha > 1){
    Malpha2 =  fem_matrices[["G"]] + m_alpha * fem_matrices[["G_2"]]/(kappa^2)
    for(j in 2:m_alpha){
      Malpha2 = Malpha2 +  choose(m_alpha, j) * fem_matrices[[paste0("G_",j+1)]]/(kappa^(2*j))
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  Q = 1/r[1] * (Malpha + Malpha2/kappa^2 - p[1] * Malpha)
  
  if(length(r)>1){
    for (i in 2:length(r)) {
      Q = c(Q, 1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
    }
  }
  
  # add k_part into Q
  
  if(sharp){
    if (m_alpha == 1){
      Kpart = 1/k * (fem_matrices[["C_less"]] + fem_matrices[["G_less"]]/(kappa^2))
    } else if (m_alpha > 1){
      Kpart = 1/k * fem_matrices[["C_less"]] + 1/k * m_alpha * fem_matrices[["G_less"]]/(kappa^2)
      for(j in 2:m_alpha){
        Kpart = Kpart + 1/k * choose(m_alpha, j) * fem_matrices[[paste0("G_",j,"_less")]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
  } else{
    if (m_alpha == 1){
      Kpart = 1/k * (fem_matrices[["C"]] + fem_matrices[["G"]]/(kappa^2))
    } else if (m_alpha > 1){
      Kpart = 1/k * fem_matrices[["C"]] + 1/k * m_alpha * fem_matrices[["G"]]/(kappa^2)
      for(j in 2:m_alpha){
        Kpart = Kpart + 1/k * choose(m_alpha, j) * fem_matrices[[paste0("G_",j)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
  }
  
  

  
  Q = c(Q, Kpart)
  
  Q = Q * kappa ^ (4 * beta)
  
  Q = tau ^ 2 * Q
  
  if(!is.null(graph)){
    graph = as(graph, "dgTMatrix")
    idx <- which(graph@i <= graph@j)
    Q = Matrix::sparseMatrix(i=graph@i[idx], j = graph@j[idx], x= Q,
                             symmetric=TRUE, index1 = FALSE)
  }
  
  return(Q)
}

#' @name rspde.matern.precision
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param rspde_order The order of the rational approximation
#' @param dim The dimension
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @return The precision matrix
#' @export

rspde.matern.precision = function(kappa, nu, tau, rspde_order, dim, fem_mesh_matrices) {
  
  n_m <- rspde_order
  
  mt <- get(paste0("m", n_m, "t"))
  
  beta = nu / 2 + dim / 4
  
  m_alpha = max(1, floor(2 * beta))
  
  r = sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("r", i)]], cut_decimals(2*beta))$y
  })
  
  p = sapply(1:(n_m), function(i) {
    approx(mt$alpha, mt[[paste0("p", i)]], cut_decimals(2*beta))$y
  })
  
  k = approx(mt$alpha, mt$k, cut_decimals(2*beta))$y
  
  if (m_alpha == 1){
    Malpha =  (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
  } else if (m_alpha > 1){
    Malpha =  fem_mesh_matrices[["c0"]] + m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
    for(j in 2:m_alpha){
      Malpha = Malpha +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  
  if (m_alpha == 1){
    Malpha2 =  (fem_mesh_matrices[["g1"]] + fem_mesh_matrices[["g2"]]/(kappa^2))
  } else if (m_alpha > 1){
    Malpha2 =  fem_mesh_matrices[["g1"]] + m_alpha * fem_mesh_matrices[["g2"]]/(kappa^2)
    for(j in 2:m_alpha){
      Malpha2 = Malpha2 +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j+1)]]/(kappa^(2*j))
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  Q = 1/r[1] * (Malpha + Malpha2/kappa^2 - p[1] * Malpha)
  
  if(length(r)>1){
    for (i in 2:length(r)) {
      Q = bdiag(Q, 1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
    }
  }
  
  # add k_part into Q
  

    if (m_alpha == 1){
      Kpart = 1/k * (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
    } else if (m_alpha > 1){
      Kpart = 1/k * fem_mesh_matrices[["c0"]] + 1/k * m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
      for(j in 2:m_alpha){
        Kpart = Kpart + 1/k * choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
  
  Q = bdiag(Q, Kpart)
  
  Q = Q * kappa ^ (4 * beta)
  
  Q = tau ^ 2 * Q
  
  return(Q)
}


#' @name rspde.matern.precision.integer.opt
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param d The dimension
#' @param fem_matrices A list containing the FEM-related matrices. The list should contain elements C, G, G_2, G_3, etc.
#' @param graph The sparsity graph of the matrices. If NULL, only a vector
#' of the elements will be returned, if non-NULL, a sparse matrix will be returned.
#' @return The precision matrix
#' @export

rspde.matern.precision.integer.opt = function(kappa, nu, tau, d, fem_matrices, graph=NULL) {
  
  beta = nu / 2 + d / 4
  
  n_beta = floor(2*beta)
  
  if (n_beta == 1){
    Q = (kappa^2 * fem_matrices[["C_less"]] + fem_matrices[["G_less"]])
  } else if (n_beta > 1){
    Q = kappa^(2*n_beta) * fem_matrices[["C_less"]] + n_beta * kappa^(2*(n_beta-1))* fem_matrices[["G_less"]]
    for(j in 2:n_beta){
      Q = Q + kappa^(2*(n_beta-j)) * choose(n_beta, j) * fem_matrices[[paste0("G_",j,"_less")]]
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  Q = tau ^ 2 * Q
  
  if(!is.null(graph)){
    graph = as(graph, "dgTMatrix")
    idx <- which(graph@i <= graph@j)
    Q = Matrix::sparseMatrix(i=graph@i[idx], j = graph@j[idx], x= Q,
                             symmetric=TRUE, index1 = FALSE)
  }
  
  return(Q)
}

#' @name rspde.matern.precision.integer
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param dim The dimension
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @return The precision matrix
#' @export

rspde.matern.precision.integer = function(kappa, nu, tau, dim, fem_mesh_matrices) {
  
  beta = nu / 2 + dim / 4
  
  n_beta = floor(2*beta)
  
  if (n_beta == 1){
    Q = (kappa^2 * fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]])
  } else if (n_beta > 1){
    Q = kappa^(2*n_beta) * fem_mesh_matrices[["c0"]] + n_beta * kappa^(2*(n_beta-1)) * fem_mesh_matrices[["g1"]]
    for(j in 2:n_beta){
      Q = Q + kappa^(2*(n_beta-j)) * choose(n_beta, j) * fem_mesh_matrices[[paste0("g",j)]]
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  Q = tau ^ 2 * Q
  
  return(Q)
}

#' @name rspde.make.A
#' @title Create A matrices for rSPDE models
#' @description Create A matrices for INLA-based rSPDE models
#' @param mesh An INLA mesh, optional
#' @param loc Locations, needed if an INLA mesh is provided
#' @param A The A matrix from the standard SPDE approach. Should only be provided if an
#' INLA mesh is not provided.
#' @param dim the dimension. Should only be provided if an
#' INLA mesh is not provided.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If NULL, then the model will assume that nu will be estimated. If
#' nu is fixed, you should provide the value of nu.
#' @param index For each observation/prediction value, an index into loc. Default is seq_len(nrow(A.loc)).
#' @param group For each observation/prediction value, an index into the group model.
#' @param repl For each observation/prediction value, the replicate index.
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @return The A matrix for rSPDE models.
#' @export
rspde.make.A <- function(mesh=NULL,
                           A = NULL,
                           dim = NULL,
                           loc = NULL,
                           rspde_order = 2, nu = NULL,
                           index=NULL,
                           group=NULL,
                           repl=1L,
                           n.group=NULL,
                           n.repl=NULL){
  if(!is.null(mesh)){
    stopifnot(inherits(mesh,"inla.mesh"))
    if(mesh$manifold == "R1"){
      dim = 1
    } else if(mesh$manifold == "R2"){
      dim = 2
    } else{
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
    }
  } else if(is.null(dim)){
    stop("If mesh is not provided, then you should provide the dimension d!")
  }
  
  if(!is.null(mesh)){
    if(is.null(loc)){
      stop("If you provided mesh, you should also provide the locations, loc.")
    }
  }
  
  if(!is.null(mesh)){
    A <- INLA::inla.spde.make.A(mesh=mesh, loc = loc,
                                index=index, group=group,
                                repl=repl, n.group=n.group,
                                n.repl=n.repl)
  } else if(is.null(A)){
    stop("If mesh is not provided, then you should provide the A matrix from
         the standard SPDE approach!")
  }


  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  if(fixed_nu){
    alpha = nu + dim/2
    integer_alpha <- (alpha%%1 == 0)
    if(integer_alpha){
      Abar = A
    } else{
      Abar = kronecker(matrix(1,1,rspde_order+1),A)
    }
    
  } else{
    Abar = kronecker(matrix(1,1,rspde_order+1),A)
  }
  return(Abar)
}


#' @name rspde.make.index
#' @title Create index for INLA-based rSPDE models
#' @description Create index for INLA-based rSPDE models
#' @param name Name
#' @param mesh An INLA mesh
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If NULL, then the model will assume that nu will be estimated. If
#' nu is fixed, you should provide the value of nu.
#' @param n.spde The number of basis functions in the mesh model. 
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @param dim the dimension. Should only be provided if an
#' INLA mesh is not provided.
#' @return index for rSPDE models.
#' @export
rspde.make.index <- function(name, n.spde=NULL, n.group = 1,
                             n.repl = 1, mesh = NULL,
                             rspde_order = 2, nu = NULL, dim = NULL){
  
  if(is.null(n.spde)&&is.null(mesh)){
    stop("You should provide either n.spde or mesh!")
  }
  
  if(!is.null(mesh)){
    n_mesh = mesh$n
    
    if(mesh$manifold == "R1"){
      dim = 1
    } else if(mesh$manifold == "R2"){
      dim = 2
    } else{
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
    }
    
  } else{
    n_mesh = n.spde
    if(is.null(dim)){
      stop("You should provide the dimension d!")
    }
  }
  
  name.group <- paste(name, ".group", sep = "")
  name.repl <- paste(name, ".repl", sep = "")
  
  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  
  if(fixed_nu){
    alpha = nu + dim/2
    integer_alpha <- (alpha%%1 == 0)
    
    if(integer_alpha){
      factor_rspde = 1
    } else{
      factor_rspde = rspde_order+1
    }
    
  } else{
    factor_rspde = rspde_order+1
  }
  
  out <- list()
  out[[name]] <- as.vector(sapply(1:factor_rspde, function(i){rep(rep(((i-1)*n_mesh+1):(i*n_mesh), times = n.group), times = n.repl)}))
  out[[name.group]] <- rep(rep(rep(1:n.group, each = n_mesh), times = n.repl), times = factor_rspde)
  out[[name.repl]] <- rep(rep(1:n.repl, each = n_mesh * n.group), times = factor_rspde)
  return(out)
}

#' @name rspde.precision
#' @title Precision matrices for rSPDE models Calculates the precision matrix 
#' for given parameter values based on an rspde model object.
#' @description Precision matrices for rSPDE models
#' @param rspde An rspde object.
#' @param theta Vector of log(kappa), log(nu)-dim/2, log(tau).
#' @param optimized Logical indicating if only the elements of the precision
#' matrix should be returned.
#' @return Precision matrix for the rspde model.
#' @export
rspde.precision <- function(rspde, theta, optimized = FALSE){
  check_class_inla_rspde(rspde)
  stopifnot(is.logical(optimized))
  if(length(theta)!=length(rspde$f$rgeneric$definition(cmd="initial"))){
      stop("Length of theta is incorrect!")
  }
  if(optimized){
    if(rspde$rspde$f$rgeneric$optimize){
      return(rspde$f$rgeneric$definition(cmd="Q", theta=theta))
    } else{
      return(rspde$f$rgeneric$definition(cmd="Q", theta=theta)@x)
    }
  } else{
    if(rspde$f$rgeneric$optimize){
      graph = rspde$f$rgeneric$definition(cmd="graph")
      entries = rspde$f$rgeneric$definition(cmd="Q", theta=theta)
      return(build_sparse_matrix_rspde(entries=entries, graph=graph))
    } else{
    return(rspde$f$rgeneric$definition(cmd="Q", theta=theta))
  }
  }
}

#' @name rspde.result
#' @title rSPDE result extraction from INLA estimation results
#' @description Extract field and parameter values and distributions for an rspde effect from an inla result object.
#' @param inla An inla object.
#' @param name A character string with the name of the rSPDE effect in the inla formula.
#' @param rspde The inla.rspde object used for the effect in the inla formula. 
#' @return Returns a list containing posterior results.
#' @export
rspde.result <- function(inla, name, rspde)
{
  check_class_inla_rspde(rspde)

result = list()

if(!rspde$est_nu){
  row_names = c("tau", "kappa")
} else{
  row_names <- c("tau","kappa","nu")
}


result$summary.values = inla$summary.random[[name]]

if (!is.null(inla$marginals.random[[name]])) {
  result$marginals.values = inla$marginals.random[[name]]
}


result$summary.log.tau = INLA::inla.extract.el(inla$summary.hyperpar, 
                                           paste("Theta1 for ", name, "$", sep = ""))
rownames(result$summary.log.tau) <- "log(tau)"
result$summary.log.kappa = INLA::inla.extract.el(inla$summary.hyperpar, 
                                             paste("Theta2 for ", name, "$", sep = ""))
rownames(result$summary.log.kappa) <- "log(kappa)"
if(rspde$est_nu){
  result$summary.log.nu = INLA::inla.extract.el(inla$summary.hyperpar, 
                                            paste("Theta3 for ", name, "$", sep = ""))
  rownames(result$summary.log.nu) <- "log(nu-lower_bound)"
}

if (!is.null(inla$marginals.hyperpar[[paste0("Theta1 for ", name)]])) {
  result$marginals.log.tau = INLA::inla.extract.el(inla$marginals.hyperpar, 
                                             paste("Theta1 for ", name, "$", sep = ""))
  names(result$marginals.log.tau) = "tau"
  result$marginals.log.kappa = INLA::inla.extract.el(inla$marginals.hyperpar, 
                                               paste("Theta2 for ", name, "$", sep = ""))
  names(result$marginals.log.kappa) = "kappa"
  
  if(rspde$est_nu){
    result$marginals.log.nu = INLA::inla.extract.el(inla$marginals.hyperpar,
                                            paste("Theta3 for ", name, "$", sep = ""))
    names(result$marginals.log.nu) = "nu"
  }

    result$marginals.tau = lapply(result$marginals.log.tau, 
                                  function(x) INLA::inla.tmarginal(function(y) exp(y), 
                                                             x))
    result$marginals.kappa = lapply(result$marginals.log.kappa, 
                                    function(x) INLA::inla.tmarginal(function(y) exp(y), 
                                                               x))
    if(rspde$est_nu){
      result$marginals.nu = lapply(result$marginals.log.nu, 
                                      function(x) INLA::inla.tmarginal(function(y) exp(y)+rspde$nu_lower_bound, 
                                                                 x))
  }
}

result$summary.tau <- create_summary_from_density(result$marginals.tau$tau, name="tau")
result$summary.kappa <- create_summary_from_density(result$marginals.kappa$kappa, name="kappa")
if(rspde$est_nu){
  result$summary.nu <- create_summary_from_density(result$marginals.nu$nu, name="nu")
}


class(result) <- "rspde.result"
return(result)
}

#' @name plot.rspde.result
#' @title Posterior plots for rSPDE field parameters
#' @description Posterior plots for rSPDE field parameters
#' @param x A rspde.result object.
#' @param which For which parameters the posterior should be plotted?
#' @param caption captions to appear above the plots; character vector or list of valid graphics annotations. Can be set to "" or NA to suppress all captions.
#' @param sub.caption	common title-above the figures if there are more than one. 
#' @param type_plot what type of plot should be drawn. The default is 'l'.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param main character; title to be placed at each plot additionally (and above) all captions.
#' @param cex.caption	controls the size of caption.
#' @param cex.oma.main controls the size of the sub.caption only if that is above the figures when there is more than one.
#' @param ylab Label for y axis.
#' @param xlab Label for x axis.
#' @param ... Additional arguments.
#' @return Called for its side effects.
#' @export
#' @method plot rspde.result

plot.rspde.result <- function(x, which = c("tau","kappa","nu"),
         caption = list("Posterior density for tau",
                        "Posterior density for kappa",
                        "Posterior density for nu"),
         sub.caption = NULL,
         type_plot = "l",
         ask = prod(graphics::par("mfcol")) <
           length(which) && grDevices::dev.interactive(),
         main = "",
         cex.oma.main = 1.25,
         cex.caption = 1,
         ylab = "Density",
         xlab = "x",
         ...){
  result = x
  which = which[which%in%c("tau","kappa","nu")]
  stopifnot(!is.null(which))
  
  one.fig <- prod(graphics::par("mfcol")) == 1
  
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  
  getCaption <- function(k) if (length(caption) < k)
    NA_character_
  else {
      grDevices::as.graphicsAnnot(caption[[k]])
  }

  if("nu"%in%which){
    if(is.null(result$marginals.nu)){
      which <- which[which!="nu"]
    }
  }
  
  for(i in 1:length(which)){
    graph_temp <- result[[paste0("marginals.",which[i])]][[which[i]]]
    graphics::plot(graph_temp, type = type_plot, main = main, ylab = ylab, xlab=xlab,
                   ...)
    graphics::mtext(getCaption(i), side = 3, cex = cex.caption)
    
    if (one.fig)
      graphics::title(sub = sub.caption, ...)
  }
  
  if (!one.fig && graphics::par("oma")[3L] >= 1)
    graphics::mtext(sub.caption, outer = TRUE, cex = 1.25)
  
  grDevices::dev.flush()
  
  invisible()
}


#' @name summary.rspde.result
#' @title Summary for posteriors of rSPDE field parameters
#' @description Summary for posteriors of rSPDE field parameters
#' @param object A rspde.result object.
#' @param digits integer, used for number formatting with signif()
#' @param ... Currently not used.
#' @return Returns a data frame
#' containing the summary.
#' @export
#' @method summary rspde.result
#' 
summary.rspde.result <- function(object,digits=6,...){
  out <- object$summary.tau
  out <- rbind(out, object$summary.kappa)
  if(!is.null(object$summary.nu)){
    out <- rbind(out, object$summary.nu)
  }
  return(signif(out,digits=digits))
}




#' @name rspde.mesh.project
#' @title Calculate a lattice projection to/from an inla.mesh() for rSPDE objects
#' @aliases rspde.mesh.project rspde.mesh.projector rspde.mesh.project.inla.mesh rspde.mesh.project.rspde.mesh.projector rspde.mesh.project.inla.mesh.1d
#' @description Calculate a lattice projection to/from an inla.mesh() for rSPDE objects
#' @param mesh An inla.mesh() or inla.mesh.1d() object.
#' @param nu The smoothness parameter. If NULL, it will be assumed that nu was estimated.
#' @param rspde_order The order of the rational approximation.
#' @param loc	Projection locations. Can be a matrix or a SpatialPoints or a SpatialPointsDataFrame object.
#' @param field Basis function weights, one per mesh basis function, describing the function to be avaluated at the projection locationssFunction values for on the mesh
#' @param projector A rspde.mesh.projector object.
#' @param lattice An inla.mesh.lattice() object.
#' @param xlim X-axis limits for a lattice. For R2 meshes, defaults to covering the domain.
#' @param ylim Y-axis limits for a lattice. For R2 meshes, defaults to covering the domain.
#' @param dims Lattice dimensions.
#' @param projection One of c("default", "longlat", "longsinlat", "mollweide").
#' @param ... Additional parameters.
#' @return A list with projection information for rspde.mesh.project. For rspde.mesh.projector(mesh, ...), 
#' a rspde.mesh.projector object. For rspde.mesh.project(projector, field, ...), a field projected from the mesh onto the locations 
#' given by the projector object.
#' @details This function is built upon the inla.mesh.project and inla.mesh.projector functions from INLA.
#' @rdname rspde.mesh.project
#' @export
#' 
rspde.mesh.project <- function(...) {
  UseMethod("rspde.mesh.project", ...)
}

#' @rdname rspde.mesh.project
#' @export

rspde.mesh.projector <- function(mesh,
                                 nu=NULL,
                                 rspde_order=2,
                                 loc = NULL,
                                 lattice = NULL,
                                 xlim = NULL,
                                 ylim = NULL,
                                 dims = c(100, 100),
                                 projection = NULL,
                                 ...){
  args_list <- list()
  args_list[["mesh"]] = mesh
  if(!is.null(loc)){
    args_list[["loc"]] = loc
  }
  if(!is.null(lattice)){
    args_list[["lattice"]] = lattice
  }
  if(!is.null(xlim)){
    args_list[["xlim"]] = xlim
  }
  if(!is.null(ylim)){
    args_list[["ylim"]] = ylim
  }
  if(!is.null(projection)){
    args_list[["projection"]] = projection
  }
  args_list[["dims"]] = dims
  out <- do.call(INLA::inla.mesh.projector, args_list)
  if(mesh$manifold == "R1"){
    dim = 1
  } else if(mesh$manifold == "R2"){
    dim = 2
  } else{
    stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
  }
  
  out$proj$A <- rspde.make.A(A=out$proj$A,rspde_order = rspde_order, dim = dim,
                             nu=nu)
  
  class(out) <- c(class(out), "rspde.mesh.projector")
  return(out)
}


#' @rdname rspde.mesh.project
#' @method rspde.mesh.project inla.mesh
#' @export

rspde.mesh.project.inla.mesh <- function (mesh, loc = NULL, 
                                          field = NULL, rspde_order = 2,
                                          nu = NULL, ...) {
  stopifnot(inherits(mesh, "inla.mesh"))

    if (!missing(field) && !is.null(field)) {
      proj <- rspde.mesh.projector(mesh, loc = loc, rspde_order = rspde_order, nu=nu,
                                   ...)
      return(INLA::inla.mesh.project(proj, field = field))
    }
    jj <- which(rowSums(matrix(is.na(as.vector(loc)), nrow = nrow(loc), 
                               ncol = ncol(loc))) == 0)
    smorg <- (INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = loc[jj, 
                                                                            , drop = FALSE]))
    ti <- matrix(0L, nrow(loc), 1)
    b <- matrix(0, nrow(loc), 3)
    ti[jj, 1L] <- smorg$p2m.t
    b[jj, ] <- smorg$p2m.b
    ok <- (ti[, 1L] > 0L)
    ti[ti[, 1L] == 0L, 1L] <- NA
    ii <- which(ok)
    A <- (sparseMatrix(dims = c(nrow(loc), mesh$n), i = rep(ii, 
                                                            3), j = as.vector(mesh$graph$tv[ti[ii, 1L], ]), x = as.vector(b[ii, 
                                                            ])))
    
    if(!is.null(nu)){
      if(!is.numeric(nu)){
        stop("nu must be numeric!")
      }
    }
    
    fixed_nu = !is.null(nu)
    if(fixed_nu){
      alpha = nu + 1
      integer_alpha <- (alpha%%1 == 0)
      if(integer_alpha){
        Abar = A
      } else{
        Abar = kronecker(matrix(1,1,rspde_order+1),A)
      }
      
    } else{
      Abar = kronecker(matrix(1,1,rspde_order+1),A)
    }
    
    list(t = ti, bary = b, A = Abar, ok = ok)
}
  

#' @rdname rspde.mesh.project
#' @method rspde.mesh.project rspde.mesh.projector
#' @export
#' 

rspde.mesh.project.rspde.mesh.projector <- function(projector, field, ...){
  return(INLA::inla.mesh.project(projector=projector,field=field,...))
}



#' @rdname rspde.mesh.project
#' @method rspde.mesh.project inla.mesh.1d
#' @export
#' 

rspde.mesh.project.inla.mesh.1d <- function(mesh, loc, field = NULL,
                                            rspde_order=2, nu=NULL,...){
  stopifnot(inherits(mesh, "inla.mesh.1d"))
  if (!missing(field) && !is.null(field)) {
    proj <- rspde.mesh.projector(mesh, loc, rspde_order=rspde_order, nu = nu,...)
    return(INLA::inla.mesh.project(proj, field))
  }
  A <- INLA::inla.mesh.1d.A(mesh, loc = loc)
  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  if(fixed_nu){
    alpha = nu + 1/2
    integer_alpha <- (alpha%%1 == 0)
    if(integer_alpha){
      Abar = A
    } else{
      Abar = kronecker(matrix(1,1,rspde_order+1),A)
    }
    
  } else{
    Abar = kronecker(matrix(1,1,rspde_order+1),A)
  }
  return(list(A = Abar, ok = (loc >= mesh$interval[1]) & (loc <= 
                                                         mesh$interval[2])))
}


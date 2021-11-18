utils::globalVariables(c("C", "C_inv", "C_inv_G", "G", "d", "loc", "n",
                         "n_m", "nu", "nu_lower_bound", "nu_upper_bound",
                         "do_optimize", "idx_symmetric", "n_Q", "rspde_order",
                         "graph_opt","fem_matrices", "sharp",
                         "prior.kappa", "prior.nu", "prior.tau"))

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
    return(c(0, 0, 0))
  }
  
  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      lkappa = theta[1L],
      lnu = theta[2L],
      ltau = theta[3L]
    ))
  }
  
  ######## precision matrix
  Q = function(n, theta) {
    param = interpret.theta(n, theta)
    nu = min(exp(param$lnu) + nu_lower_bound, nu_upper_bound)
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    
    if(do_optimize){
      return(rSPDE::rspde_Prec_opt(kappa=kappa, nu=nu, tau=tau, 
                                   rspde_order=rspde_order, 
                                   d=d, fem_matrices = fem_matrices,
                                   sharp=sharp, graph = NULL))
    } else{
      return(rSPDE::rspde_Prec(kappa=kappa, nu=nu, tau=tau, 
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
    return(c(0, 0))
  }
  
  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      lkappa = theta[1L],
      ltau = theta[2L]
    ))
  }
  
  ######## precision matrix
  Q = function(n, theta) {
    param = interpret.theta(n, theta)
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    
    if(do_optimize){
      return(rSPDE::rspde_Prec_opt(kappa=kappa, nu=nu, tau=tau, 
                                   rspde_order=rspde_order, 
                                   d=d, fem_matrices = fem_matrices,
                                   sharp=sharp, graph = NULL))
    } else{
      return(rSPDE::rspde_Prec(kappa=kappa, nu=nu, tau=tau, 
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
    return(c(0, 0))
  }
 
  ######## parameter
  interpret.theta <- function(n, theta) {
    return(list(
      lkappa = theta[1L],
      ltau = theta[2L]
    ))
  }
  
  ######## precision matrix
  Q = function(n, theta) {
    param = interpret.theta(n, theta)
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    

    
    if(do_optimize){
      return(rSPDE::rspde_Prec_opt_int(kappa, nu, tau,
                                       d, fem_matrices, graph = NULL))
    } else{
      return(rSPDE::rspde_Prec_int(kappa=kappa, nu=nu, tau=tau, 
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



#' @name create_rspde_model
#' @title Create INLA-based rSPDE models
#' @description Create INLA-based rSPDE models with general smoothness
#' @param inla_mesh An INLA mesh
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
#' @return An INLA model.
#' @export

create_rspde_model <- function(inla_mesh, 
                               nu_lower_bound = 
                                 get_inla_mesh_dimension(inla_mesh)/4,
                               nu_upper_bound = 4, rspde_order = 2,
                               nu = NULL, sharp = TRUE,
                               debug = FALSE,
                               optimize=TRUE, 
                               prior.kappa = list(meanlog = 0, sdlog = 1),
                               prior.nu = list(meanlog = 0, sdlog = 1),
                               prior.tau = list(meanlog = 0, sdlog = 1)
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
  
  if(inla_mesh$manifold == "R1"){
    d = 1
  } else if(inla_mesh$manifold == "R2"){
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
      fem_mesh <- INLA::inla.mesh.fem(inla_mesh, order = m_alpha)
    } else{
      fem_mesh <- INLA::inla.mesh.fem(inla_mesh, order = m_alpha+1)
    }
  } else{
    beta = nu_order / 2 + d / 4
    
    m_alpha = max(1, floor(2 * beta))
  
    if(integer_alpha){
      fem_matrices <- INLA::inla.mesh.fem(inla_mesh, order = m_alpha)
    } else{
      fem_matrices <- INLA::inla.mesh.fem(inla_mesh, order = m_alpha+1)
    }
  }
  
  if(optimize){
    if(integer_alpha){
      result_sparsity <- analyze_sparsity_rspde(nu_upper_bound=nu_order, d=d, 
                                                rspde_order=rspde_order, 
                                                fem_mesh_matrices = fem_mesh,
                                                include_higher_order = FALSE)
    } else{
      if(sharp){
        result_sparsity <- analyze_sparsity_rspde(nu_upper_bound=nu_order, d=d, 
                                                  rspde_order=rspde_order, 
                                                  fem_mesh_matrices = fem_mesh)
      } else{
        result_sparsity <- analyze_sparsity_rspde(nu_upper_bound=nu_order, d=d, 
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
       graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_mesh, d=d,
                                       nu=nu_upper_bound,
                                       rspde_order = rspde_order,
                                       sharp = sharp,
                                       force_non_integer = TRUE)
     } else{
       graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_matrices, d=d,
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
                                        m1t=m1t,m2t = m2t, m3t=m3t, 
                                        m4t=m4t, m5t=m5t, m6t=m6t,
                                        m7t=m7t,m8t=m8t,
                                        prior.kappa=prior.kappa,
                                        prior.nu=prior.nu,
                                        prior.tau=prior.tau,
                                        d = d, rspde_order = rspde_order, 
                                        n=ncol(C)*(rspde_order+1), 
                                        debug=debug,
                                        do_optimize=optimize, optimize=optimize)
  } else if(!integer_alpha){
    if(optimize){
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_mesh, d=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   sharp = sharp)
    } else{
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_matrices, d=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   sharp = TRUE, force_non_integer=TRUE)
    }

    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_frac_alpha,
                                        nu = nu,
                                        fem_matrices=fem_matrices, 
                                        graph_opt=graph_opt,
                                        sharp=sharp,
                                        m1t=m1t,m2t = m2t, m3t=m3t, 
                                        m4t=m4t, m5t=m5t, m6t=m6t,
                                        m7t=m7t,m8t=m8t,
                                        prior.kappa=prior.kappa,
                                        prior.nu=prior.nu,
                                        prior.tau=prior.tau,
                                        d = d, rspde_order = rspde_order, 
                                        n=ncol(C)*(rspde_order+1), 
                                        debug=debug,
                                        do_optimize=optimize, optimize=optimize)
  } else{
    if(optimize){
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_mesh, d=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   sharp = sharp)
    } else{
      graph_opt <- rSPDE::get_sparsity_graph_rspde(fem_mesh_matrices = fem_matrices, d=d,
                                                   nu=nu,
                                                   rspde_order = rspde_order,
                                                   force_non_integer=FALSE)
    }
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_int_alpha,
                                        nu = nu,
                                        fem_matrices=fem_matrices, 
                                        graph_opt=graph_opt, 
                                        m1t=m1t,m2t = m2t, m3t=m3t, 
                                        m4t=m4t, m5t=m5t, m6t=m6t,
                                        m7t=m7t,m8t=m8t,
                                        prior.kappa=prior.kappa,
                                        prior.nu=prior.nu,
                                        prior.tau=prior.tau,
                                        d = d,
                                        n=ncol(C), 
                                        debug=debug,
                                        do_optimize=optimize, optimize=optimize)
  }
  return(model)
}



symmetric_part_matrix <- function(M){
  M  <- as(M, "dgTMatrix")
  idx <- which(M@i <= M@j)
  sM <- cbind(M@i[idx], M@j[idx])
  colnames(sM) <- NULL
  return(list(M = split(sM, seq(nrow(sM))), idx = idx))
}



#' @name analyze_sparsity_rspde
#' @title Analyze sparsity of matrices in the rSPDE approach
#' @description Auxiliar function to analyze sparsity of matrices in the rSPDE approach
#' @param nu_upper_bound Upper bound for the smoothness parameter
#' @param d The dimension of the domain
#' @param rspde_order The order of the rational approximation
#' @param fem_mesh_matrices A list containing FEM-related matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @param include_lower_order Logical. Should the lower-order terms be included? They are needed for the cases
#' when alpha = nu + d/2 is integer or for when sharp is set to TRUE.
#' @param include_higher_order Logical. Should be included for when nu is estimated or for when alpha = nu + d/2 is not an integer.
#' @return A list containing informations on sparsity of the precision matrices
#' @export

analyze_sparsity_rspde <- function(nu_upper_bound, d, rspde_order,
                                   fem_mesh_matrices,
                                   include_lower_order = TRUE, 
                                   include_higher_order = TRUE){
  beta = nu_upper_bound / 2 + d / 4
  
  m_alpha = max(1, floor(2 * beta))
  
  positions_matrices <- list()
  
  C_list <- symmetric_part_matrix(fem_mesh_matrices$c0)
  G_1_list <- symmetric_part_matrix(fem_mesh_matrices$g1)
  for(j in 2:(m_alpha)){
      assign(paste0("G_",j,"_list"), symmetric_part_matrix(fem_mesh_matrices[[paste0("g",j)]]))
  }
    
  if(include_higher_order){
    assign(paste0("G_",m_alpha+1,"_list"), symmetric_part_matrix(fem_mesh_matrices[[paste0("g",m_alpha+1)]]))
      
    positions_matrices[[1]] <-  match(C_list$M, get(paste0("G_",m_alpha+1,"_list"))[["M"]])
  }
  

  
  idx_matrices <- list()
  
  idx_matrices[[1]] = C_list$idx
  
  for(i in 1:m_alpha){
    if(include_higher_order){
      positions_matrices[[i+1]] <- match(get(paste0("G_",i,"_list"))[["M"]], get(paste0("G_",m_alpha+1,"_list"))[["M"]])
    }
    idx_matrices[[i+1]] = get(paste0("G_",i,"_list"))[["idx"]]
  }
  
  if(include_higher_order){
    idx_matrices[[m_alpha+2]] = get(paste0("G_",m_alpha+1,"_list"))[["idx"]]
  }
  
  if(include_lower_order){
    positions_matrices_less <- list()
    positions_matrices_less[[1]] <-  match(C_list$M, get(paste0("G_",m_alpha,"_list"))[["M"]])
    if(m_alpha > 1){
      for(i in 1:(m_alpha-1)){
        positions_matrices_less[[i+1]] <- match(get(paste0("G_",i,"_list"))[["M"]], get(paste0("G_",m_alpha,"_list"))[["M"]])
      }
    } else{
      positions_matrices_less[[2]] <-  1:length(get(paste0("G_",m_alpha,"_list"))[["M"]])
    }
  } else{
    positions_matrices_less <- NULL
  }
  
  return(list(positions_matrices = positions_matrices,
              idx_matrices = idx_matrices,
              positions_matrices_less = positions_matrices_less))
}



#' @name rspde_Prec_opt
#' @title Optimized precision matrix of the covariance-based rational approximation
#' @description Computes the optimized version of the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param rspde_order The order of the rational approximation
#' @param d The dimension
#' @param fem_matrices A list containing the FEM-related matrices. The list should contain elements C, G, G_2, G_3, etc.
#' @param graph The sparsity graph of the matrices. If NULL, only a vector
#' of the elements will be returned, if non-NULL, a sparse matrix will be returned.
#' @param sharp The sparsity graph should have the correct sparsity (costs more to perform
#' a sparsity analysis) or an upper bound for the sparsity?
#' @return The precision matrix
#' @export

rspde_Prec_opt = function(kappa, nu, tau, rspde_order, d, fem_matrices, graph=NULL,
                          sharp) {
  
  n_m <- rspde_order
  
  mt <- get(paste0("m", n_m, "t"))
  
  beta = nu / 2 + d / 4
  
  m_alpha = max(1, floor(2 * beta))
  
  r = sapply(1:(n_m), function(i) {
    approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
  })
  
  p = sapply(1:(n_m), function(i) {
    approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
  })
  
  k = approx(mt$nu, mt$k, cut_decimals(nu))$y
  
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

#' @name rspde_Prec
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param rspde_order The order of the rational approximation
#' @param d The dimension
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @return The precision matrix
#' @export

rspde_Prec = function(kappa, nu, tau, rspde_order, d, fem_mesh_matrices) {
  
  n_m <- rspde_order
  
  mt <- get(paste0("m", n_m, "t"))
  
  beta = nu / 2 + d / 4
  
  m_alpha = max(1, floor(2 * beta))
  
  r = sapply(1:(n_m), function(i) {
    approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
  })
  
  p = sapply(1:(n_m), function(i) {
    approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
  })
  
  k = approx(mt$nu, mt$k, cut_decimals(nu))$y
  
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


#' @name rspde_Prec_opt_int
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

rspde_Prec_opt_int = function(kappa, nu, tau, d, fem_matrices, graph=NULL) {
  
  beta = nu / 2 + d / 4
  
  n_beta = floor(2*beta)
  
  if (n_beta == 1){
    Q = (kappa^2 * fem_matrices[["C_less"]] + fem_matrices[["G_less"]])
  } else if (n_beta > 1){
    Q = (kappa^2) * fem_matrices[["C_less"]] + n_beta * fem_matrices[["G_less"]]
    for(j in 2:n_beta){
      Q = Q + kappa^(2*j) * choose(n_beta, j) * fem_matrices[[paste0("G_",j,"_less")]]
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

#' @name rspde_Prec_int
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param kappa Range parameter of the latent process.
#' @param nu The smoothness parameter
#' @param tau The precision parameter
#' @param d The dimension
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @return The precision matrix
#' @export

rspde_Prec_int = function(kappa, nu, tau, d, fem_mesh_matrices) {
  
  beta = nu / 2 + d / 4
  
  n_beta = floor(2*beta)
  
  if (n_beta == 1){
    Q = (kappa^2 * fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]])
  } else if (n_beta > 1){
    Q = (kappa^2) * fem_mesh_matrices[["c0"]] + n_beta * fem_mesh_matrices[["g1"]]
    for(j in 2:n_beta){
      Q = Q + kappa^(2*j) * choose(n_beta, j) * fem_mesh_matrices[[paste0("g",j)]]
    }
  } else {
    stop("Something is wrong with the value of nu!")
  }
  
  Q = tau ^ 2 * Q
  
  return(Q)
}

#' @name get_inla_mesh_dimension
#' @title Get the dimension of an INLA mesh
#' @description Get the dimension of an INLA mesh
#' @param inla_mesh An INLA mesh
#' @return The dimension of an INLA mesh.
#' @export
get_inla_mesh_dimension <- function(inla_mesh){
  if(!(class(inla_mesh) %in% c("inla.mesh", "inla.mesh.1d"))){
    stop("The object should be an INLA mesh")
  } else{
    if(inla_mesh$manifold == "R1"){
      d = 1
    } else if(inla_mesh$manifold == "R2"){
      d = 2
    } else{
      stop("The mesh should be from a flat manifold.")
    }
  }
}



#' @name get_sparsity_graph_rspde
#' @title Sparsity graph for rSPDE models
#' @description Creates the sparsity graph for rSPDE models
#' @param inla_mesh An INLA mesh, optional
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The list should contain elements C, G, G_2, G_3, etc. Optional,
#' should be provided if inla_mesh is not provided.
#' @param d The dimension, optional. Should be provided if inla_mesh is not provided.
#' @param nu The smoothness parameter
#' @param force_non_integer Should nu be treated as non_integer?
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param sharp The graph should have the correct sparsity (costs more to perform
#' a sparsity analysis) or an upper bound for the sparsity?
#' @return The A matrix for rSPDE models.
#' @export
get_sparsity_graph_rspde <- function(inla_mesh=NULL,
                                     fem_mesh_matrices=NULL,
                                     nu,
                                     force_non_integer = FALSE,
                                     rspde_order = 2,
                                     sharp = TRUE,
                                     d=NULL){
  
  if(!is.null(inla_mesh)){
    if(inla_mesh$manifold == "R1"){
      d = 1
    } else if(inla_mesh$manifold == "R2"){
      d = 2
    } else{
      stop("The mesh should be from a flat manifold.")
    }
  } else if(is.null(d)){
    stop("If an INLA mesh is not provided, you should provide the dimension!")
  }
  
  alpha = nu + d/2
  
  m_alpha = max(1, floor(alpha))
  
  integer_alpha <- (alpha%%1 == 0)
  
  if(force_non_integer){
    integer_alpha = FALSE
  }
  
  if(!is.null(fem_mesh_matrices)){
    if(integer_alpha){
      return(fem_mesh_matrices[[paste0("g",m_alpha)]])
    } else{
      if(sharp){
        return(bdiag(kronecker(diag(rep(1,rspde_order)),
                               fem_mesh_matrices[[paste0("g",m_alpha+1)]]), 
                     fem_mesh_matrices[[paste0("g",m_alpha)]]))
      } else{
        return(kronecker(diag(rep(1,rspde_order+1)),
                         fem_mesh_matrices[[paste0("g",m_alpha+1)]]))
      }
    }
  } else if(!is.null(inla_mesh)){
    if(integer_alpha){
      fem_mesh_matrices <- INLA::inla.mesh.fem(inla_mesh, order = m_alpha)
      return(fem_mesh_matrices[[paste0("g",m_alpha)]])
    } else{
      fem_mesh_matrices <- INLA::inla.mesh.fem(inla_mesh, order = m_alpha+1)
      if(sharp){
        return(bdiag(kronecker(diag(rep(1,rspde_order)),
                               fem_mesh_matrices[[paste0("g",m_alpha+1)]]), 
                     fem_mesh_matrices[[paste0("g",m_alpha)]]))
      } else{
        return(kronecker(diag(rep(1,rspde_order+1)),
                         fem_mesh_matrices[[paste0("g",m_alpha+1)]]))
      }
    } 
  } else {
    stop("You should provide either inla_mesh or fem_mesh_matrices!")
  } 
}


#' @name rspde_create_A
#' @title Create A matrices for rSPDE models
#' @description Create A matrices for INLA-based rSPDE models
#' @param inla_mesh An INLA mesh, optional
#' @param loc Locations, needed if an INLA mesh is provided
#' @param A The A matrix from the standard SPDE approach. Should only be provided if an
#' INLA mesh is not provided.
#' @param d the dimension. Should only be provided if an
#' INLA mesh is not provided.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If NULL, then the model will assume that nu will be estimated. If
#' nu is fixed, you should provide the value of nu.
#' @return The A matrix for rSPDE models.
#' @export
rspde_create_A <- function(inla_mesh=NULL,
                           A = NULL,
                           d = NULL,
                           loc = NULL,
                           rspde_order = 2, nu = NULL){
  if(!is.null(inla_mesh)){
    if(inla_mesh$manifold == "R1"){
      d = 1
    } else if(inla_mesh$manifold == "R2"){
      d = 2
    } else{
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
    }
  } else if(is.null(d)){
    stop("If inla_mesh is not provided, then you should provide the dimension d!")
  }
  
  if(!is.null(inla_mesh)){
    if(is.null(loc)){
      stop("If you provided inla_mesh, you should also provide the locations, loc.")
    }
  }
  
  if(!is.null(inla_mesh)){
    A <- INLA::inla.spde.make.A(inla_mesh, loc = loc)
  } else if(is.null(A)){
    stop("If inla_mesh is not provided, then you should provide the A matrix from
         the standard SPDE approach!")
  }


  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  if(fixed_nu){
    alpha = nu + d/2
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


#' @name rspde_make_index
#' @title Create index for INLA-based rSPDE models
#' @description Create index for INLA-based rSPDE models
#' @param name Name
#' @param inla_mesh An INLA mesh
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If NULL, then the model will assume that nu will be estimated. If
#' nu is fixed, you should provide the value of nu.
#' @return index for rSPDE models.
#' @export
rspde_make_index <- function(name, inla_mesh,
                             rspde_order = 2, nu = NULL){
  n_mesh <- inla_mesh$n
  
  if(inla_mesh$manifold == "R1"){
    d = 1
  } else if(inla_mesh$manifold == "R2"){
    d = 2
  } else{
    stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
  }
  
  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  
  if(fixed_nu){
    alpha = nu + d/2
    integer_alpha <- (alpha%%1 == 0)
    
    if(integer_alpha){
      n_rspde = n_mesh
    } else{
      n_rspde = (rspde_order+1)*n_mesh
    }
    
  } else{
    n_rspde = (rspde_order+1)*n_mesh
  }
  return(INLA::inla.spde.make.index(name = name, n.spde = n_rspde))
}

#' @name build_sparse_matrix_rspde
#' @title Create sparse matrix from entries and graph
#' @description Create sparse matrix from entries and graph
#' @param entries The entries of the precision matrix
#' @param graph The sparsity graph of the precision matrix
#' @return index for rSPDE models.
#' @export
build_sparse_matrix_rspde <- function(entries, graph){
  if(!is.null(graph)){
    graph = as(graph, "dgTMatrix")
    idx <- which(graph@i <= graph@j)
    Q = Matrix::sparseMatrix(i=graph@i[idx], j = graph@j[idx], x= entries,
                             symmetric=TRUE, index1 = FALSE)
  }
  return(Q)
}

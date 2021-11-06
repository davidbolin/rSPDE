#Please install pardiso for inla:
#For further details visit: https://www.pardiso-project.org/r-inla/#license
#inla.setOption(pardiso.license = '~/.pardiso.license',
#               smtp = 'pardiso') 


utils::globalVariables(c("C", "C_inv", "C_inv_G", "G", "d", "loc", "n",
                         "n_m", "nu", "nu_lower_bound", "nu_upper_bound",
                         "do_optimize", "idx_symmetric", "n_Q", "positions",
                         "lengths","epsilon_identity"))

#' @importFrom stats dnorm pnorm
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
  
  cut_decimals <- function(nu) {
    temp <- nu - floor(nu)
    if (temp < 10 ^ (-3)) {
      temp <- 10 ^ (-3)
    }
    if (temp > 0.999) {
      temp <- 0.999
    }
    return(temp)
  }
  
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
    nu = exp(param$lnu) + nu_lower_bound
    tau = exp(param$ltau)
    kappa = exp(param$lkappa)
    
    I_C = Diagonal(n = nrow(C), x = 1)
    
    mt <- get(paste0("m", n_m, "t"))
    
    if (nu > nu_upper_bound) {
      nu = nu_upper_bound
    }
    
    beta = nu / 2 + d / 4
    
    m_alpha = max(1, floor(2 * beta))
    
    Id_temp = diag(rep(1, n_m + 1))
    
    r = sapply(1:(n_m), function(i) {
      approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
    })
    
    p = sapply(1:(n_m), function(i) {
      approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
    })
    
    k = approx(mt$nu, mt$k, cut_decimals(nu))$y
    
    # matrix for operator L
    L = C + G / kappa ^ 2
    sizeL = dim(L)[1]
    
    # compute m_alpha part, the interger order of L
    Malpha = I_C
    
    for (i in 1:m_alpha) {
      Malpha = Malpha %*% (I_C + C_inv_G / (kappa ^ 2))
    }
    #construct Q matrix
    
    C_L = (I_C + C_inv_G / (kappa ^ 2))
    
    Q = ((L - p[1] * C) %*% Malpha) / r[1]
    
    for (i in 2:length(r)) {
      Q = bdiag(Q, ((L - p[i] * C) %*% Malpha) / r[i])
    }
    # add k_part into Q
    
    if (m_alpha == 0) {
      Q = bdiag(Q, (1 / k) * C_inv)
      
    } else if (m_alpha == 1) {
      Q = bdiag(Q, (1 / k) * L)
      
    } else if (m_alpha == 2) {
      Q = bdiag(Q, (1 / k) * L %*% C_L)
      
    } else if (m_alpha == 3) {
      Q = bdiag(Q, (1 / k) * L %*% C_L %*% C_L)
    } else {
      dif_m_alpha = m_alpha - 3
      Q_temp = (1 / k) * L %*% C_L %*% C_L
      while (dif_m_alpha > 0) {
        Q_temp = Q_temp %*% C_L
        dif_m_alpha = dif_m_alpha - 1
      }
      Q = bdiag(Q, Q_temp)
    }
    
    Q = Q * kappa ^ (4 * beta)
    
    Q = Q + Diagonal(n = nrow(Q), x = epsilon_identity)
    Q = tau ^ 2 * Q
    
    if(do_optimize){
      #nu_pos <- floor(nu_lower_bound) + floor(nu) + 1
      nu_pos <- which(length(Q@x) == lengths)
      idx_pos <- positions[[nu_pos]]
      Q_x <- rep(0,n_Q)
      Q_x[idx_pos] <- Q@x[idx_symmetric[[nu_pos]]]
      return(Q_x)
    }
  
    return(Q)
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
    return(rSPDE::rspde_Prec(C, G, C_inv_G, C_inv, n_m, d, c(-1, log(nu_upper_bound),-1)))
  }
  
  
  
  ######################## log prior
  log.prior <- function(n, theta) {
    param = interpret.theta(n, theta)
    
    if (param$lnu > log(nu_upper_bound - nu_lower_bound)) {
      tdnorm_nu = -Inf
    } else{
      tdnorm_nu = dnorm(param$lnu, 0, 1, log = TRUE) -
        pnorm(log(nu_upper_bound - nu_lower_bound), 0, 1, log.p = TRUE)
    }
    
    res <- tdnorm_nu + dnorm(param$lkappa, 0, 1, log = TRUE) +
      dnorm(param$ltau, 0, 1, log = TRUE)
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

#' @name analyze_sparsity_rpsde
#' @title Analyze sparsity of matrices in the rSPDE approach
#' @description Auxiliar function to analyze sparsity of matrices in the rSPDE approach
#' @param nu_lower_bound Lower bound for the smoothness parameter
#' @param nu_upper_bound Upper bound for the smoothness parameter
#' @param d The dimension of the domain
#' @param C The mass matrix
#' @param G The stiffness matrix
#' @param C_inv Inverse of the mass matrix
#' @param C_inv_G C^{-1}G
#' @param rspde_order The order of the rational approximation
#' @return A list containing informations on sparsity of the precision matrices
#' @export

analyze_sparsity_rpsde <- function(nu_lower_bound, nu_upper_bound, d, 
                                   C, G, C_inv, C_inv_G, rspde_order){

  
  if(nu_lower_bound==floor(nu_lower_bound)){
    values_theta <- log(2e-3+c(nu_lower_bound:(floor(nu_upper_bound-1e-5))))
  } else{
    values_theta <- log(2e-3+c(nu_lower_bound, ceiling(nu_lower_bound):(floor(nu_upper_bound-1e-5))))
  }
  
  #values_theta <- log(c(nu_lower_bound,nu_lower_bound + 1:(floor(nu_upper_bound)-1) + nu_upper_bound-floor(nu_upper_bound)))
  
  positions = list()
  
  Prec_v1 <- rspde_Prec(C,G,C_inv_G,C_inv, rspde_order, d, c(0, log(nu_upper_bound), 0))
  Prec_v1  <- INLA::inla.as.sparse(Prec_v1)
  length_upper <- length(Prec_v1@x)
  idx_upper_bound <- which(Prec_v1@i <= Prec_v1@j)
  Prec_v1@i <- Prec_v1@i[idx_upper_bound]
  Prec_v1@j <- Prec_v1@j[idx_upper_bound]
  sM_v1 <- cbind(Prec_v1@i, Prec_v1@j)
  colnames(sM_v1) <- NULL
  
  idx_symmetric <- list()
  
  sM_v1 <- split(sM_v1, seq(nrow(sM_v1)))
  
  lengths <- list()
  
  for(i in 1:length(values_theta)){
    Prec_v2 <- rspde_Prec(C,G,C_inv_G,C_inv, rspde_order, d, c(0, values_theta[i], 0))
    Prec_v2  <- INLA::inla.as.sparse(Prec_v2)
    lengths[[i]] <- length(Prec_v2@x)
    idx <- which(Prec_v2@i <= Prec_v2@j)
    Prec_v2@i <- Prec_v2@i[idx]
    Prec_v2@j <- Prec_v2@j[idx]
    sM_v2 <- cbind(Prec_v2@i, Prec_v2@j)
    colnames(sM_v2) <- NULL
    
    sM_v2 <- split(sM_v2, seq(nrow(sM_v2)))
    
    positions[[i]] = match(sM_v2, sM_v1)
    idx_symmetric[[i]] = idx
  }
  n_Q <- length(idx_upper_bound)
  positions[[length(values_theta)+1]] <- sM_v1
  lengths[[length(values_theta)+1]] <- length_upper
  return(list(positions=positions,n_Q=n_Q, idx_symmetric=idx_symmetric,
              idx_upper_bound=idx_upper_bound, lengths=lengths))
}



#' @name rspde_Prec
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param C The mass matrix
#' @param G The stiffness matrix
#' @param C_inv_G C^{-1}G
#' @param C_inv The inverse of the mass matrix
#' @param n_m The order of the rational approximation.
#' @param d The dimension
#' @param theta vector containing the parameters of the model.
#' @param epsilon_identity a multiple of the identity matrix to be summed to 
#' the precision matrix.
#' @return The precision matrix
#' @export

rspde_Prec = function(C, G, C_inv_G, C_inv, n_m, d, theta, epsilon_identity=1e-5) {
  
  cut_decimals <- function(nu) {
    temp <- nu - floor(nu)
    if (temp < 10 ^ (-3)) {
      temp <- 10 ^ (-3)
    }
    if (temp > 0.999) {
      temp <- 0.999
    }
    return(temp)
  }
  
  interpret.theta <- function(theta) {
    return(list(
      lkappa = theta[1L],
      lnu = theta[2L],
      ltau = theta[3L]
    ))
  }
  
  
  
  param = interpret.theta(theta)
  nu = exp(param$lnu)
  tau = exp(param$ltau)
  kappa = exp(param$lkappa)
  
  I_C = Diagonal(n = nrow(C), x = 1)
  
  mt <- get(paste0("m", n_m, "t"))
  
  beta = nu / 2 + d / 4
  
  m_alpha = max(1, floor(2 * beta))

  Id_temp = diag(rep(1, n_m + 1))
  
  r = sapply(1:(n_m), function(i) {
    approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
  })
  
  p = sapply(1:(n_m), function(i) {
    approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
  })
  
  k = approx(mt$nu, mt$k, cut_decimals(nu))$y
  
  # matrix for operator L
  L = C + G / kappa ^ 2
  sizeL = dim(L)[1]
  
  # compute m_alpha part, the interger order of L
  Malpha = I_C
  
  for (i in 1:m_alpha) {
    Malpha = Malpha %*% (I_C + C_inv_G / (kappa ^ 2))
  }
  #construct Q matrix
  
  C_L = (I_C + C_inv_G / (kappa ^ 2))
  
  Q = ((L - p[1] * C) %*% Malpha) / r[1]
  
  for (i in 2:length(r)) {
    Q = bdiag(Q, ((L - p[i] * C) %*% Malpha) / r[i])
  }
  # add k_part into Q
  
  if (m_alpha == 0) {
    Q = bdiag(Q, (1 / k) * C_inv)
    
  } else if (m_alpha == 1) {
    Q = bdiag(Q, (1 / k) * L)
    
  } else if (m_alpha == 2) {
    Q = bdiag(Q, (1 / k) * L %*% C_L)
    
  } else if (m_alpha == 3) {
    Q = bdiag(Q, (1 / k) * L %*% C_L %*% C_L)
  } else {
    dif_m_alpha = m_alpha - 3
    Q_temp = (1 / k) * L %*% C_L %*% C_L
    while (dif_m_alpha > 0) {
      Q_temp = Q_temp %*% C_L
      dif_m_alpha = dif_m_alpha - 1
    }
    Q = bdiag(Q, Q_temp)
  }
  
  Q = Q * kappa ^ (4 * beta)
  
  Q = Q + Diagonal(n = nrow(Q), x = epsilon_identity)
  Q = tau ^ 2 * Q
  
  return(Q)
}



#' @name rspde_Prec_int
#' @title Precision matrix of the covariance-based rational approximation
#' @description Computes the precision matrix for the covariance-based rational SPDE
#' @param C The mass matrix
#' @param G The stiffness matrix
#' @param C_inv_G C^{-1}G
#' @param C_inv The inverse of the mass matrix
#' @param d The dimension of the domain
#' @param theta vector containing the parameters of the model.
#' @param nu the smoothness parameter
#' @param epsilon_identity a multiple of the identity matrix to be summed to 
#' the precision matrix.
#' @return The precision matrix
#' @export

rspde_Prec_int = function(C, G, C_inv_G, C_inv, d, theta, nu, epsilon_identity=1e-5) {
  cut_decimals <- function(nu) {
    temp <- nu - floor(nu)
    if (temp < 10 ^ (-3)) {
      temp <- 10 ^ (-3)
    }
    if (temp > 0.999) {
      temp <- 0.999
    }
    return(temp)
  }
  
  interpret.theta <- function(theta) {
    return(list(
      lkappa = theta[1L],
      ltau = theta[2L]
    ))
  }
  
  param = interpret.theta(theta)
  tau = exp(param$ltau)
  kappa = exp(param$lkappa)
  
  I_C = Diagonal(n = nrow(C), x = 1)
  
  beta = nu / 2 + d / 4
  
  m_alpha = max(1, floor(2 * beta))
  
  n_beta = as.integer(2*beta)-1
  
  K = kappa^2*C + G
  Q = K
  if(n_beta>0){
    CinvK = kappa^2*I_C + C_inv_G
  }
  
  while(n_beta>0){
    Q = Q %*% CinvK
    n_beta = n_beta - 1
  }
  
  Q = Q + Diagonal(n = nrow(Q), x = epsilon_identity)
  Q = tau ^ 2 * Q
  
  return(Q)
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
  
  envir = parent.env(environment())
  if (!exists("cache.done", envir = envir)) {
    Q <- rSPDE::rspde_Prec(C, G, C_inv_G, C_inv, n_m, d, c(-1, log(nu),-1))
    Q <- INLA::inla.as.sparse(Q)
    idx <- which(Q@i <= Q@j)
    assign("idx", idx, envir = envir)
    assign("cache.done", TRUE, envir = envir)
  }
  
  cut_decimals <- function(nu) {
    temp <- nu - floor(nu)
    if (temp < 10 ^ (-3)) {
      temp <- 10 ^ (-3)
    }
    if (temp > 0.999) {
      temp <- 0.999
    }
    return(temp)
  }
  
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
    
    I_C = Diagonal(n = nrow(C), x = 1)
    
    mt <- get(paste0("m", n_m, "t"))
    
    beta = nu / 2 + d / 4
    
    m_alpha = max(1, floor(2 * beta))
    sigma2_spde = gamma(nu) / (gamma(nu + d / 2) * (4 * pi) ^ (0.5) * kappa ^
                                 (2 * nu))
    
    Id_temp = diag(rep(1, n_m + 1))
    
    r = sapply(1:(n_m), function(i) {
      approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
    })
    
    p = sapply(1:(n_m), function(i) {
      approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
    })
    
    k = approx(mt$nu, mt$k, cut_decimals(nu))$y
    
    # matrix for operator L
    L = C + G / kappa ^ 2
    sizeL = dim(L)[1]
    
    # compute m_alpha part, the interger order of L
    Malpha = I_C
    
    for (i in 1:m_alpha) {
      Malpha = Malpha %*% (I_C + C_inv_G / (kappa ^ 2))
    }
    #construct Q matrix
    
    C_L = (I_C + C_inv_G / (kappa ^ 2))
    
    Q = ((L - p[1] * C) %*% Malpha) / r[1]
    
    for (i in 2:length(r)) {
      Q = bdiag(Q, ((L - p[i] * C) %*% Malpha) / r[i])
    }
    # add k_part into Q
    
    if (m_alpha == 0) {
      Q = bdiag(Q, (1 / k) * C_inv)
      
    } else if (m_alpha == 1) {
      Q = bdiag(Q, (1 / k) * L)
      
    } else if (m_alpha == 2) {
      Q = bdiag(Q, (1 / k) * L %*% C_L)
      
    } else if (m_alpha == 3) {
      Q = bdiag(Q, (1 / k) * L %*% C_L %*% C_L)
    } else {
      dif_m_alpha = m_alpha - 3
      Q_temp = (1 / k) * L %*% C_L %*% C_L
      while (dif_m_alpha > 0) {
        Q_temp = Q_temp %*% C_L
        dif_m_alpha = dif_m_alpha - 1
      }
      Q = bdiag(Q, Q_temp)
    }
    
    Q = Q * kappa ^ (4 * beta) * sigma2_spde
    
    Q = Q + Diagonal(n = nrow(Q), x = epsilon_identity)
    Q = tau ^ 2 * Q
    
    if(do_optimize){
      return(Q@x[idx])
    }
    
    return(Q)
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
    return(rSPDE::rspde_Prec(C, G, C_inv_G, C_inv, n_m, d, c(-1, log(nu),-1)))
  }
  
  
  
  ######################## log prior
  log.prior <- function(n, theta) {
    param = interpret.theta(n, theta)
    
    res <- dnorm(param$lkappa, 0, 1, log = TRUE) +
      dnorm(param$ltau, 0, 1, log = TRUE)
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
  
  envir = parent.env(environment())
  if (!exists("cache.done", envir = envir)) {
    Q <- rSPDE::rspde_Prec_int(C, G, C_inv_G, C_inv, d, c(-1, -1),nu)
    Q <- INLA::inla.as.sparse(Q)
    idx <- which(Q@i <= Q@j)
    assign("idx", idx, envir = envir)
    assign("cache.done", TRUE, envir = envir)
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
    
    I_C = Diagonal(n = nrow(C), x = 1)
    
    beta = nu / 2 + d / 4
    
    m_alpha = max(1, floor(2 * beta))
    sigma2_spde = gamma(nu) / (gamma(nu + d / 2) * (4 * pi) ^ (0.5) * kappa ^
                                 (2 * nu))
    
    n_beta = as.integer(2*beta)-1
    
    K = kappa^2*C + G
    Q = K * sigma2_spde
    if(n_beta>0){
      CinvK = kappa^2*I_C + C_inv_G
    }
    
    while(n_beta>0){
      Q = Q %*% CinvK
      n_beta = n_beta - 1
    }
    
    Q = Q + Diagonal(n = nrow(Q), x = epsilon_identity)
    Q = tau ^ 2 * Q
    
    if(do_optimize){
      return(Q@x[idx])
    }
    
    return(Q)
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
    return(rSPDE::rspde_Prec_int(C, G, C_inv_G, C_inv, d, c(-1, -1),nu))
  }
  
  
  
  ######################## log prior
  log.prior <- function(n, theta) {
    param = interpret.theta(n, theta)
    
    res <- dnorm(param$lkappa, 0, 1, log = TRUE) +
      dnorm(param$ltau, 0, 1, log = TRUE)
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

#' @name create_rspde_model
#' @title Create INLA-based rSPDE models
#' @description Create INLA-based rSPDE models with general smoothness
#' @param inla_mesh An INLA mesh
#' @param nu_lower_bound Lower bound for the smoothness parameter.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not be estimated.
#' If nu is NULL, it will be estimated.
#' @param debug INLA debug argument.
#' @param optimize Should the model be optimized? In this case the sparsities of
#' the matrices will be analyzed.
#' @param controls A list of controls. For now, the only control is
#' sum.identity, which sums a multiple of the identity matrix (the value of 
#' the multiple is indicated by the user) to the
#' precision parameter to improve numerical stability.
#' @return An INLA model.
#' @export
create_rspde_model <- function(inla_mesh, 
                               nu_lower_bound = 
                                 get_inla_mesh_dimension(inla_mesh)/4,
                               nu_upper_bound = 7, rspde_order = 2,
                               nu = NULL, debug = FALSE,
                               optimize=FALSE, controls = list()
){
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
  
  if(!is.null(controls[["sum.identity"]])){
    epsilon_identity <- controls[["sum.identity"]]
  } else{
    epsilon_identity <- 1e-5
  }

  
  fem_mesh <- INLA::inla.mesh.fem(inla_mesh)
  C <- fem_mesh$c0
  G <- fem_mesh$g1
  C_inv <- solve(C)
  C_inv_G <- solve(C,G)
  alpha = nu + d/2
  
  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  
  if(!fixed_nu){
    if(optimize){
      result_sparsity <- analyze_sparsity_rpsde(nu_lower_bound, nu_upper_bound, d, 
                                                C, G, C_inv, C_inv_G, rspde_order) 
      positions <- result_sparsity$positions
      n_Q <- result_sparsity$n_Q
      idx_symmetric <- result_sparsity$idx_symmetric
      idx_upper_bound <- result_sparsity$idx_upper_bound
      lengths <- result_sparsity$lengths
    } else{
      positions <- NULL
      n_Q <- NULL
      idx <- NULL
      idx_upper_bound <- NULL
      idx_symmetric <- NULL
    }
  }
  
  
  if(fixed_nu){
  if(alpha - floor(alpha) == 0){
    integer_alpha = TRUE
  } else{
    integer_alpha = FALSE
  }
  }
  if(!fixed_nu){
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_general,
                                  nu_lower_bound=nu_lower_bound,
                                  nu_upper_bound=nu_upper_bound, 
                                  positions=positions, n_Q = n_Q,
                                  idx_symmetric=idx_symmetric,
                                  lengths=lengths,
                                  idx_upper_bound=idx_upper_bound,
                                  epsilon_identity=epsilon_identity,
                                  C = C, C_inv = C_inv,
                                  C_inv_G=C_inv_G, G = G,
                                  m1t=m1t,m2t = m2t, m3t=m3t, 
                                  m4t=m4t, m5t=m5t, m6t=m6t,
                                  m7t=m7t,m8t=m8t,
                                  d = d, n_m = rspde_order, 
                                  n=ncol(C)*(rspde_order+1), 
                                  debug=debug,
                                  do_optimize=optimize, optimize=optimize)
  } else if(!integer_alpha){
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_frac_alpha,
                                  nu = nu,
                                  C = C, C_inv = C_inv,
                                  C_inv_G=C_inv_G, G = G, 
                                  m1t=m1t,m2t = m2t, m3t=m3t, 
                                  m4t=m4t, m5t=m5t, m6t=m6t,
                                  m7t=m7t,m8t=m8t, 
                                  epsilon_identity=epsilon_identity,
                                  d = d, n_m = rspde_order, 
                                  n=ncol(C)*(rspde_order+1), 
                                  debug=debug,
                                  do_optimize=optimize, optimize=optimize)
  } else{
    model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_int_alpha,
                                  nu = nu,
                                  C = C, C_inv = C_inv,
                                  C_inv_G=C_inv_G, G = G, 
                                  m1t=m1t,m2t = m2t, m3t=m3t, 
                                  m4t=m4t, m5t=m5t, m6t=m6t,
                                  m7t=m7t,m8t=m8t,
                                  epsilon_identity=epsilon_identity,
                                  d = d,
                                  n=ncol(C), 
                                  debug=debug,
                                  do_optimize=optimize, optimize=optimize)
  }
  return(model)
}

#' @name rspde_create_A
#' @title Create A matrices for INLA-based rSPDE models
#' @description Create A matrices for INLA-based rSPDE models
#' @param inla_mesh An INLA mesh
#' @param loc Locations
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not be estimated.
#' If nu is NULL, it will be estimated.
#' @return The A matrix for rSPDE models.
#' @export
rspde_create_A <- function(inla_mesh, loc,
                           rspde_order = 2, nu = NULL){
  if(inla_mesh$manifold == "R1"){
    d = 1
  } else if(inla_mesh$manifold == "R2"){
    d = 2
  } else{
    stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
  }
  A <- INLA::inla.spde.make.A(inla_mesh, loc = loc)
  if(!is.null(nu)){
    if(!is.numeric(nu)){
      stop("nu must be numeric!")
    }
  }
  
  fixed_nu = !is.null(nu)
  if(fixed_nu){
    alpha = nu + d/2
    if(alpha - floor(alpha) == 0){
      integer_alpha = TRUE
    } else{
      integer_alpha = FALSE
    }
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
#' @param nu If nu is set to a parameter, nu will be kept fixed and will not be estimated.
#' If nu is NULL, it will be estimated.
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
    if(alpha - floor(alpha) == 0){
      integer_alpha = TRUE
    } else{
      integer_alpha = FALSE
    }
    
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
  


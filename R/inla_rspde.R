#Please install pardiso for inla:
#For further details visit: https://www.pardiso-project.org/r-inla/#license
#inla.setOption(pardiso.license = '~/.pardiso.license',
#               smtp = 'pardiso') 

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
    sigma2_spde = gamma(nu) / (gamma(nu + d / 2) * (4 * pi) ^ (0.5) * kappa ^
                                 (2 * nu))
    
    Id_temp = diag(rep(1, 2 * n_m + 1))
    
    r = sapply(1:(2 * n_m), function(i) {
      approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
    })
    
    p = sapply(1:(2 * n_m), function(i) {
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
    
    Q = Q + Diagonal(n = nrow(Q), x = 1e-5)
    Q = tau ^ 2 * Q
    
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
    return(Q(n, c(-1, log(nu_upper_bound),-1)))
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
    
    Id_temp = diag(rep(1, 2 * n_m + 1))
    
    r = sapply(1:(2 * n_m), function(i) {
      approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
    })
    
    p = sapply(1:(2 * n_m), function(i) {
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
    
    Q = Q + Diagonal(n = nrow(Q), x = 1e-5)
    Q = tau ^ 2 * Q
    
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
    return(Q(n, c(-1, -1)))
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
    
    Q = Q + Diagonal(n = nrow(Q), x = 1e-5)
    Q = tau ^ 2 * Q
    
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
    return(Q(n, c(-1, -1)))
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

#' @name create_rspde_model
#' @title Create INLA-based rSPDE models
#' @description Create INLA-based rSPDE models with general smoothness
#' @param inla_mesh An INLA mesh
#' @param d dimension of the data
#' @param nu_lower_bound Lower bound for the smoothness parameter.
#' @param nu_upper_bound Upper bound for the smoothness parameter.
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param fixed_nu If TRUE, nu will be kept fixed and will not be estimated,
#' in this case the value of nu must be provided in the nu argument. If FALSE
#' nu will be estimated. 
#' @param nu The smoothness parameter to be used if fixed_nu is TRUE.
#' @param debug INLA debug argument.
#' @return An INLA model.
#' @export
create_rspde_model <- function(inla_mesh, d, nu_lower_bound = d/4,
                               nu_upper_bound = 7, rspde_order = 2,
                               fixed_nu = FALSE, 
                               nu = NULL, debug = FALSE
){
  if(nu_lower_bound<d/4){
    stop("The lower bound for nu cannot be less than d/4")
  }
  fem_mesh <- inla.mesh.fem(inla_mesh)
  C <- fem_mesh$c0
  G <- fem_mesh$g1
  C_inv <- solve(C)
  C_inv_G <- solve(C,G)
  alpha = nu + d/2
  if(fixed_nu){
    if(is.null(nu)){
      stop("For fixed nu, nu must be given!")
    }
  if(alpha - floor(alpha) == 0){
    integer_alpha = TRUE
  } else{
    integer_alpha = FALSE
  }
  }
  if(!fixed_nu){
    model <- inla.rgeneric.define(inla.rgeneric.cov_rspde_general,
                                  nu_lower_bound=nu_lower_bound,
                                  nu_upper_bound=nu_upper_bound, 
                                  C = C, C_inv = C_inv,
                                  C_inv_G=C_inv_G, G = G, 
                                  m2t = m2t, m3t=m3t, m4t=m4t, 
                                  d = d, n_m = rspde_order, 
                                  n=ncol(C)*(2*rspde_order+1), 
                                  debug=debug)
  } else if(!integer_alpha){
    model <- inla.rgeneric.define(inla.rgeneric.cov_rspde_frac_alpha,
                                  nu = nu,
                                  C = C, C_inv = C_inv,
                                  C_inv_G=C_inv_G, G = G, 
                                  m2t = m2t, m3t=m3t, m4t=m4t, 
                                  d = d, n_m = rspde_order, 
                                  n=ncol(C)*(2*rspde_order+1), 
                                  debug=debug)
  } else{
    model <- inla.rgeneric.define(inla.rgeneric.cov_rspde_int_alpha,
                                  nu = nu,
                                  C = C, C_inv = C_inv,
                                  C_inv_G=C_inv_G, G = G, 
                                  d = d,
                                  n=ncol(C), 
                                  debug=debug)
  }
  return(model)
}

#' @name rspde_create_A
#' @title Create A matrices for INLA-based rSPDE models
#' @description Create A matrices for INLA-based rSPDE models
#' @param inla_mesh An INLA mesh
#' @param loc Locations
#' @param d dimension of the data
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param fixed_nu If TRUE, nu will be kept fixed and will not be estimated,
#' in this case the value of nu must be provided in the nu argument. If FALSE
#' nu will be estimated. 
#' @param nu The smoothness parameter to be used if fixed_nu is TRUE.
#' @return The A matrix for rSPDE models.
#' @export
rspde_create_A <- function(inla_mesh, loc, fixed_nu = FALSE,
                           d = NULL,
                           rspde_order = 2, nu = NULL){
  A <- inla.spde.make.A(inla_mesh, loc = loc)
  if(fixed_nu){
    if(is.null(nu)){
      stop("For fixed nu, nu must be given!")
    }
    if(is.null(d)){
      stop("For fixed nu, d must be given!")
    }
    alpha = nu + d/2
    if(alpha - floor(alpha) == 0){
      integer_alpha = TRUE
    } else{
      integer_alpha = FALSE
    }
    
    if(integer_alpha){
      Abar = A
    } else{
      Abar = kronecker(matrix(1,1,2*rspde_order+1),A)
    }
    
  } else{
    Abar = kronecker(matrix(1,1,2*rspde_order+1),A)
  }
  return(Abar)
}


#' @name rspde_make_index
#' @title Create index for INLA-based rSPDE models
#' @description Create index for INLA-based rSPDE models
#' @param name Name
#' @param inla_mesh An INLA mesh
#' @param d dimension of the data
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param fixed_nu If TRUE, nu will be kept fixed and will not be estimated,
#' in this case the value of nu must be provided in the nu argument. If FALSE
#' nu will be estimated. 
#' @param nu The smoothness parameter to be used if fixed_nu is TRUE.
#' @return index for rSPDE models.
#' @export
rspde_make_index <- function(name, inla_mesh, fixed_nu = FALSE,
                             rspde_order = 2, nu = NULL, d = NULL){
  n_mesh <- inla_mesh$n
  
  if(fixed_nu){
    if(is.null(nu)){
      stop("For fixed nu, nu must be given!")
    }
    if(is.null(d)){
      stop("For fixed nu, d must be given!")
    }
    alpha = nu + d/2
    if(alpha - floor(alpha) == 0){
      integer_alpha = TRUE
    } else{
      integer_alpha = FALSE
    }
    
    if(integer_alpha){
      n_rspde = n_mesh
    } else{
      n_rspde = (2*rspde_order+1)*n_mesh
    }
    
  } else{
    n_rspde = (2*rspde_order+1)*n_mesh
  }
  return(inla.spde.make.index(name = name, n.spde = n_rspde))
}
  






#' @noRd 
kappa_integral <- function(n,beta,kappa){
  y <- 0
  for(k in 0:n){
    y = y + (-1)^(k+1)*choose(n,k)/(n-k-beta+1)  
  }
  return(kappa^(2*(n-beta+1))*y)
}

#' @noRd 

markov_approx <- function(beta,kappa,d){
  nu <- 2*beta - d/2
  alpha <- nu + d/2
  L <- alpha - floor(alpha)
  p <- ceiling(alpha)
  B <- matrix(0,nrow=p+1,ncol=p+1)
  c <- rep(p+1,1)
  for (i in 0:p){
    c[i+1] <- 2*kappa_integral(i,-alpha+2*p+1+L,kappa)
    for (j in 0:p){
      B[j+1,i+1] = 2*kappa_integral(i+j,2*p+1+L,kappa)
    }
  }
  b <- solve(solve(diag(sqrt(diag(B))),B),solve(diag(sqrt(diag(B))),c))
  return(b*gamma(nu)/(gamma(alpha)*(4*pi)^(d/2)*kappa^(2*nu)))
}

#' @noRd 

markov.Q <- function(beta,kappa,d,fem){
  b <- markov_approx(beta,kappa,d)
  Q <- b[1]*fem$C
  M <- fem$C
  for(i in 2:length(b)){
    M <- M%*%solve(fem$C,fem$G)
    Q <- Q + b[i]*M
  }
  return(Q)
}

# kappa <- 10
# sigma <- 1/sqrt(2)
# nu <- 0.8
# beta <- (nu + 1/2)/2

# nobs <- 101
# x <- seq(from = 0, to = 1, length.out = 101)
# fem <- rSPDE.fem1d(x)

# Q <- 2*markov.Q(beta,kappa,d=1,fem)
# v <- t(rSPDE.A1d(x, 0.5))
# A <- Diagonal(nobs)
# c_cov.approx <- A %*% solve(Q, v)
# c.true <- folded.matern.covariance.1d(rep(0.5, length(x)),
#                                       abs(x), kappa, nu, sigma)

# # plot the result and compare with the true Matern covariance
# plot(x, c.true,
#      type = "l", ylab = "C(h)",
#      xlab = "h", main = "Matern covariance and rational approximations"
# )
# lines(x, c_cov.approx, col = 2)



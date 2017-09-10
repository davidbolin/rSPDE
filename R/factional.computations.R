rSPDE.sample <- function(obj,n=1)
  {
  if(class(obj) != "rSPDEobj")
    stop("Input op is not of class rSPDEobj")
  m = dim(obj$Q)[1]
  z = rnorm(n*m)
  dim(z) <- c(m,n)
  if(obj$commutative){
    R = chol(obj$Q)
    x <- obj$L2%*%solve(R,z)
  } else {
    R = chol(obj$Q)
    x <- solve(obj$L1,t(R)%*%z)
  }
  return(x)
}





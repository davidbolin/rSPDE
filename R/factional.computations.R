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

rSPDE.krig <- function(obj,A,Aprd,Y,nugget)
  {
  if(class(obj) != "rSPDEobj")
    stop("Input op is not of class rSPDEobj")

  if(obj$commutative){
    A = A%*%obj$L2
    Qhat = obj$Q + t(A)%*%A/nugget
    Yhat = as.vector(Aprd%*%(obj$L2%*%solve(Qhat,t(A)%*%Y/nugget)))
  } else {
    M1 = cBind(t(A)%*%A/nugget,t(obj$L1))
    M2 = cBind(obj$L1, -obj$Q)
    M = rBind(M1,M2)
    v = c(as.vector(t(A)%*%Y/nugget),rep(0,dim(obj$Q)[1]))
    tmp = solve(M,v)
    Yhat = tmp[1:dim(obj$Q)[1]]
  }
  return(Yhat)
}

rSPDE.loglike <- function(obj,Y,A,nugget)
  {
  if(length(dim(Y))==2){
    n.rep = dim(Y)[2]
  } else {
    n.rep = 1
  }
  if(obj$commutative) {
    A = A%*%obj$L2
    Q.post = obj$Q + t(A)%*%A/nugget
    R = Matrix::Cholesky(obj$Q)
    prior.ld = 2*c(determinant(R,logarithm = TRUE)$modulus)
    R.post = Matrix::Cholesky(Q.post)
    posterior.ld = 2*c(determinant(R.post,logarithm = TRUE)$modulus)

    AtY = t(A)%*%Y/nugget
    mu.post = solve(Q.post,AtY)

    lik = n.rep*(prior.ld - posterior.ld - dim(A)[1]*(log(nugget)+log(2*pi)))/2
    if(n.rep>1){
      lik = lik - 0.5*sum(colSums(mu.post*(obj$Q%*%mu.post)))
      lik = lik - 0.5*sum(colSums((Y-A%*%mu.post)^2))/nugget
    } else {
      lik = lik - 0.5*(t(mu.post)%*%obj$Q%*%mu.post + t(Y-A%*%mu.post)%*%(Y-A%*%mu.post)/nugget)
    }
  } else {

    L1.lu = lu(op$L1)
    prior.ld = 2*sum(log(abs(diag(L1.lu@U))))

    M1 = cBind(t(A)%*%A/nugget,t(obj$L1))
    M2 = cBind(obj$L1, -obj$Q)
    M = rBind(M1,M2)
    M.lu = lu(M)
    posterior.ld = sum(log(abs(diag(M.lu@U))))

    #compute posterior mean
    v = c(as.vector(t(A)%*%Y/nugget),rep(0,dim(obj$Q)[1]))
    tmp = solve(M,v)
    mu.post = tmp[1:dim(obj$Q)[1]]

    lik = n.rep*(prior.ld - posterior.ld - dim(A)[1]*(log(nugget)+log(2*pi)))/2
    L1.post = op$L1%*%mu.post
    QL1.post = solve(op$Q,L1.post)

    if(n.rep>1){
      lik = lik - 0.5*sum(colSums(L1.post*QL1.post))
      lik = lik - 0.5*sum(colSums((Y-A%*%mu.post)^2))/nugget
    } else {
      lik = lik - 0.5*(t(L1.post)%*%QL1.post + t(Y-A%*%mu.post)%*%(Y-A%*%mu.post)/nugget)
    }
  }
  return(as.double(lik))
}


matern.loglike <- function(kappa,sigma,nu,nugget,Y,G,C,A,d=2,m=1)
{
  op <- matern.operators(kappa,sigma,nu,G,C,d=d,m=m)
  return(rSPDE.loglike(op,Y,A,nugget))
}


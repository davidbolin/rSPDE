fractional.operators <- function(L,beta,C,scale.factor,m=1,commutative=TRUE,tau=1)
  {
  roots <- get.roots(m,beta)
  bc = max(1,floor(beta))
  C = Diagonal(dim(C)[1],rowSums(C))
  Ci = Diagonal(dim(C)[1],1/rowSums(C))
  I = Diagonal(dim(C)[1])
  CiL = Ci%*%L
  L2 = I-CiL*roots$rc[1]
  if(length(roots$rc)>1){
    for(i in 2:length(roots$rc)){
      if(i==length(roots$rc)){
        L2C <- L2%*%(Diagonal(dim(C)[1],sqrt(diag(C)))-L%*%Diagonal(dim(C)[1],1/sqrt(diag(C)))%*%roots$rc[i])
      }
      L2 = L2%*%(I-CiL*roots$rc[i])
    }
  } else {
    L2C = Diagonal(dim(C)[1],sqrt(diag(C)))-L%*%Diagonal(dim(C)[1],1/sqrt(diag(C)))*roots$rc[1]
  }


  if(bc>1){
    L1 = I-CiL*roots$rb[1]
  } else {
    L1 = C-L*roots$rb[1]
  }

  CiL1 = I-CiL*roots$rb[1]
  if(length(roots$rb)>1){
    for(i in 2:length(roots$rb)){
      L1 = L1%*%(I-CiL*roots$rb[i])
      CiL1 = CiL1%*%(I-CiL*roots$rb[i])
    }
  }

  if(bc>1){
    Lp = L
    CiLp = CiL
    if(bc>2){
      for(i in 2:(bc-1)){
        Lp = CiL%*%Lp
        CiLp = CiL%*%CiLp
      }
    }
    L1 = Lp%*%L1
    CiL1 = CiLp%*%CiL1
  }
  L1 <- L1*tau*scale.factor^beta/roots$factor
  CiL1 = CiL1*tau*scale.factor^beta/roots$factor

  if(commutative){
    Q = t(L1)%*%CiL1
  } else {

    Q = L2C%*%t(L2C)
  }

  output <- list(Q = Q,
                 L1 = L1,
                 L2 = L2,
                 CiL1 = CiL1,
                 C = C,
                 m = m,
                 beta = beta,
                 type = "fractional approximation",
                 commutative = commutative)
    class(output) <- "rSPDEobj"
  return(output)
}

matern.operators <- function(kappa,sigma,nu,G,C,d=2,m=1)
{
  s = colSums(C)
  tau = sqrt(gamma(nu)/(sigma^2*kappa^(2*nu)*(4*pi)^(d/2)*gamma(nu+d/2)))
  beta = (nu + d/2)/2
  operators <- fractional.operators(L = G/kappa^2 + C,
                                    beta=beta,
                                    C=C,
                                    scale.factor = kappa^2,
                                    m=m,
                                    tau = tau)
  output <- operators
  output$kappa = kappa
  output$sigma = sigma
  output$nu = nu
  output$type = "Matern approximation"
  return(output)
}




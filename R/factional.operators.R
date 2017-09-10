fractional.operators <- function(L,beta,C,Ci,scale.factor,m=1,commutative=TRUE)
  {
  roots <- get.roots(m,beta)

  I = Diagonal(dim(C)[1])
  CiL = Ci%*%L

  L2 = C-L*roots$rc[1]
  CiL2 = I-CiL*roots$rc[1]
  if(length(roots$rc)>1){
    for(i in 2:length(roots$rc)){
      L2 <- L2%*%(I-CiL*roots$rc[i])
      CiL2 = CiL2%*%(I-CiL*roots$rc[i])
    }
  }


  if(max(1,floor(beta))>1){
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

  if(max(1,floor(beta))>1){
    Lp = L
    CiLp = CiL
    if(max(1,floor(beta))>2){
      for(i in 2:(max(1,floor(beta))-1)){
        Lp = CiL%*%Lp
        CiLp = CiL%*%CiLp
      }
    }
    L1 = Lp%*%L1
    CiL1 = CiLp%*%CiL1
  }
  L1 <- L1*scale.factor^beta/roots$factor
  CiL1 = CiL1*scale.factor^beta/roots$factor
  if(commutative){
    Q = t(L1)%*%CiL1
  } else {
    Q = L2%*%C%*%t(L2)
  }

  output <- list(Q = Q,
                 L1 = L1,
                 L2 = CiL2,
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
  s = rowSums(C)
  Ci = Diagonal(length(s),1/s)
  C = Diagonal(length(s),s)
  tau = sqrt(gamma(nu)/(sigma^2*kappa^(2*nu)*(4*pi)^(d/2)*gamma(nu+d/2)))
  beta = (nu + d/2)/2
  operators <- fractional.operators(L = G/kappa^2 + C,
                                    beta=beta,
                                    C=C,
                                    Ci=Ci,
                                    scale.factor = kappa^(2),
                                    m=m)
  output <- operators
  output$Q=(tau^2)*operators$Q
  output$kappa = kappa
  output$sigma = sigma
  output$nu = nu
  output$type = "Matern approximation"
  output$commutative = TRUE
  return(output)
}




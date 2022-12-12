system("make")

n = 10

y = rnorm(n, sd = 0.1)

i <- rep(1:10,2)
j <- c(1:10,c(1,1:9))
sp_mat <- Matrix::sparseMatrix(i = i, j= j, x =1)


cmodel <- inla.cgeneric.define(model = "inla_cgeneric_cppmodel",
                               shlib = "cgeneric_test.so", n = n,
                               sp_mat = sp_mat)


rc <- inla(
  y ~ -1 + f(idx, model = cmodel), 
  data = data.frame(y, idx = 1:n),
  verbose = TRUE,
  control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))

cat(rc$summary.fitted.values$mean)  

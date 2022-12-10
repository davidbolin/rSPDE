system("make")

n = 10

y = rnorm(n, sd = 0.1)

cmodel <- inla.cgeneric.define(model = "inla_cgeneric_cppmodel",
                               shlib = "cgeneric_test.so", n = n)


rc <- inla(
  y ~ -1 + f(idx, model = cmodel), 
  data = data.frame(y, idx = 1:n),
  control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))

cat(rc$summary.fitted.values$mean)  
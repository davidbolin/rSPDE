library(INLA)
library(rSPDE)

inla.setOption(pardiso.license = "/home/debusta/pardiso.lic")


data(PRprec)
data(PRborder)

Y <- rowMeans(PRprec[, 6 + 1:31])
ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])
alt <- PRprec$Altitude[ind]
seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)


prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 0.25), cutoff = 0.2)

# prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
# prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.35, 0.1), cutoff = 0.01)

Abar.int <- rspde.make.A(
  mesh = prmesh, loc = coords,
  nu = 1
)
mesh.index.int <- rspde.make.index(
  name = "field", mesh = prmesh,
  nu = 1
)


start <- Sys.time()
rspde_model_fix_int1 <- rspde.matern_cgeneric(mesh = prmesh,
      nu = 1
    )
end <- Sys.time()
print('Execution Time')
print(end-start)
stk.dat.int <- inla.stack(
  data = list(y = Y), A = list(Abar.int, 1), tag = "est",
  effects = list(
    c(
      mesh.index.int,
      list(Intercept = 1)
    ),
    list(
      long = inla.group(coords[, 1]),
      lat = inla.group(coords[, 2]),
      seaDist = inla.group(seaDist)
    )
  )
)

f.s.fix.int.1 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_fix_int1)

total_time = c()
for(i in 1:5){
rspde_fix_int_1 <- inla(f.s.fix.int.1,
      family = "Gamma",
      data = inla.stack.data(stk.dat.int), verbose = FALSE, debug=FALSE,
      control.inla = list(int.strategy = "eb"),
      control.predictor = list(
        A = inla.stack.A(stk.dat.int),
        compute = FALSE
      ),
                inla.mode = "experimental"
    )
total_time = c(total_time, rspde_fix_int_1$cpu.used[4][[1]])
}

cat("Mean Total time cgeneric:\n")
cat(mean(total_time))
cat("\n")

cat("Std Dev Total time cgeneric:\n")
cat(sd(total_time))
cat("\n")


# result <- inla.cgeneric.q(rspde_model_fix_int1)
# head(result$Q)


rspde_model_fix_base <- rspde.matern(mesh = prmesh,
      nu = 1
    )

start <- Sys.time()
spde_model_int <- inla.spde2.matern(mesh = prmesh, alpha=2)
end <- Sys.time()
print('Execution Time')
print(end-start)

stk.dat.int <- inla.stack(
  data = list(y = Y), A = list(Abar.int, 1), tag = "est",
  effects = list(
    c(
      mesh.index.int,
      list(Intercept = 1)
    ),
    list(
      long = inla.group(coords[, 1]),
      lat = inla.group(coords[, 2]),
      seaDist = inla.group(seaDist)
    )
  )
)

f.spde.s.fix.int.1 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = spde_model_int)


total_time = c()
for(i in 1:5){
spde_fix_int <- inla(f.spde.s.fix.int.1,
      family = "Gamma",
      data = inla.stack.data(stk.dat.int), verbose = FALSE, debug=FALSE,
      control.inla = list(int.strategy = "eb"),
      control.predictor = list(
        A = inla.stack.A(stk.dat.int),
        compute = FALSE
      ),
                inla.mode = "experimental"
    )
total_time = c(total_time, spde_fix_int$cpu.used[4][[1]])
}

cat("Mean Total time INLA:\n")
cat(mean(total_time))
cat("\n")

cat("Std Dev Total time INLA:\n")
cat(sd(total_time))
cat("\n")

# Rerun cgeneric: 


total_time = c()
for(i in 1:5){
rspde_fix_int_1 <- inla(f.s.fix.int.1,
      family = "Gamma",
      data = inla.stack.data(stk.dat.int), verbose = FALSE, debug=FALSE,
      control.inla = list(int.strategy = "eb"),
      control.predictor = list(
        A = inla.stack.A(stk.dat.int),
        compute = FALSE
      ),
                inla.mode = "experimental"
    )
total_time = c(total_time, rspde_fix_int_1$cpu.used[4][[1]])
}

cat("Mean Total time cgeneric:\n")
cat(mean(total_time))
cat("\n")

cat("Std Dev Total time cgeneric:\n")
cat(sd(total_time))
cat("\n")

# df_results <- data.frame(Processor = c("M1 Pro", "M1 Pro", "Intel i9-12900K", "Intel i9-12900K"), Precision = c("1874x1874(18602)", "1874x1874(18602)"), Model = c("INLA", "rSPDE-cgeneric"), `Mean time` = c(4.399, 3.948), `Std Dev` = c(0.188, 0.188))
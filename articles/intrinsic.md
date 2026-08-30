# Intrinsic models in the rSPDE package

### Introduction

In this vignette we provide a brief introduction to the intrinsic models
implemented in the `rSPDE` package.

### A fractional intrinsic model

A basic intrinsic model which is implemented in `rSPDE` is defined as  
``` math
(-\Delta)^{\beta/2}(\tau u) = \mathcal{W},
```
where $`\beta > d/2`$ and $`d`$ is the dimension of the spatial domain.

To illustrate these models, we begin by defining a mesh over
$`[0,2]\times [0, 2]`$:

``` r

library(fmesher)
bnd <- fm_segm(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2)), is.bnd = TRUE)
mesh_2d <- fm_mesh_2d(
    boundary = bnd, 
    cutoff = 0.02,
    max.edge = c(0.1)
)
plot(mesh_2d, main = "")
```

![](intrinsic_files/figure-html/unnamed-chunk-1-1.png)

We now use the
[`intrinsic.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.operators.md)
function to construct the `rSPDE` representation of the general model.

``` r

library(rSPDE)
tau <- 0.2
beta <- 1.8
fem <- fm_fem(mesh_2d)
op <- intrinsic.operators(tau = tau, beta = beta, mesh = mesh_2d, m = 2)
```

To see that the `rSPDE` model is approximating the true model, we can
compare the variogram of the approximation (implemented in the function
`variogram` in the model object) with the true variogram (implemented in
[`variogram.intrinsic.spde()`](https://davidbolin.github.io/rSPDE/reference/variogram.intrinsic.spde.md))
as follows.

``` r

point <- matrix(c(1,1),1,2)
Gamma <- op$variogram(point)
vario <- variogram.intrinsic.spde(point, mesh_2d$loc[,1:2], tau = tau,
                                  beta = beta, L = 2, d = 2)
d = sqrt((mesh_2d$loc[,1]-point[1])^2 +  (mesh_2d$loc[,2]-point[2])^2)
plot(d, Gamma, xlim = c(0,0.7), ylim = c(0,3),
     ylab = "variogram(h)", xlab = "h")
points(d,vario,col=2)
```

![](intrinsic_files/figure-html/unnamed-chunk-3-1.png)

If we want to increase the accuracy, we can either use a finer mesh or
increase the order of the rational approximation through the argument
`m` in `intrinsic.operators`. The default value of `m` is 1. We can now
use the `simulate` function to simulate a realization of the field
$`u`$:

``` r

u <- simulate(op,nsim = 1)

proj <- fm_evaluator(mesh_2d, dims = c(100, 100))
field <- fm_evaluate(proj, field = as.vector(u))
field.df <- data.frame(x1 = proj$lattice$loc[,1],
                       x2 = proj$lattice$loc[,2], 
                       y = as.vector(field))

library(ggplot2)
library(viridis)
#> Loading required package: viridisLite
ggplot(field.df, aes(x = x1, y = x2, fill = y)) +
    geom_raster() +
    scale_fill_viridis()
```

![](intrinsic_files/figure-html/unnamed-chunk-4-1.png)

By default, the field is simulated with a zero-integral constraint.

### Fitting the model with `R-INLA`

Let us now consider a simple Gaussian linear model where the spatial
field $`u(\mathbf{s})`$ is observed at $`m`$ locations,
$`\{\mathbf{s}_1 , \ldots , \mathbf{s}_m \}`$ under Gaussian measurement
noise. For each $`i = 1,\ldots,m,`$ we have
``` math
\begin{align} 
y_i &= u(\mathbf{s}_i)+\varepsilon_i\\
\end{align},
```
where $`\varepsilon_1,\ldots,\varepsilon_{m}`$ are iid normally
distributed with mean 0 and standard deviation 0.1.

To generate a data set `y` from this model, we first draw some
observation locations at random in the domain and then use the
[`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
functions (that wraps the functions
[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html),
[`fm_block()`](https://inlabru-org.github.io/fmesher/reference/fm_block.html)
and
[`fm_row_kron()`](https://inlabru-org.github.io/fmesher/reference/fm_row_kron.html)
of the `fmesher` package) to construct the observation matrix which can
be used to evaluate the simulated field $`u`$ at the observation
locations. After this we simply add the measurment noise.

``` r

n_loc <- 1000
loc_2d_mesh <- matrix(2*runif(n_loc * 2), n_loc, 2)

A <- spde.make.A(
  mesh = mesh_2d,
  loc = loc_2d_mesh
)
sigma.e <- 0.1
y <- A %*% u + rnorm(n_loc) * sigma.e
```

The generated data can be seen in the following image.

``` r

df <- data.frame(x1 = as.double(loc_2d_mesh[, 1]),
  x2 = as.double(loc_2d_mesh[, 2]), y = as.double(y))
ggplot(df, aes(x = x1, y = x2, col = y)) +
  geom_point() +
  scale_color_viridis()
```

![](intrinsic_files/figure-html/unnamed-chunk-6-1.png)

We will now fit the model using our [`R-INLA`](https://www.r-inla.org)
implementation of the rational SPDE approach. Further details on this
implementation can be found in [R-INLA implementation of the rational
SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md).

``` r

library(INLA)
#> 
rspde.order <- 2
mesh.index <- rspde.make.index(name = "field", mesh = mesh_2d, rspde.order = rspde.order)
Abar <- rspde.make.A(mesh = mesh_2d, loc = loc_2d_mesh, rspde.order = rspde.order)
st.dat <- inla.stack(data = list(y = as.vector(y)), A = Abar, effects = mesh.index)
```

We now create the model object.

``` r

rspde_model <- rspde.intrinsic(mesh = mesh_2d, rspde.order = rspde.order)
```

Finally, we create the formula and fit the model to the data:

``` r

f <- y ~ -1 + f(field, model = rspde_model)
rspde_fit <- inla(f,
                  data = inla.stack.data(st.dat),
                  family = "gaussian",
                  control.predictor = list(A = inla.stack.A(st.dat)))
```

To compare the estimated parameters to the true parameters, we can do
the following:

``` r

result_fit <- rspde.result(rspde_fit, "field", rspde_model)
summary(result_fit)
#>         mean        sd 0.025quant 0.5quant 0.975quant     mode
#> tau 0.124832 0.0257716  0.0823184 0.122031   0.183194 0.116496
#> nu  0.969270 0.0714821  0.8280500 0.969787   1.108280 0.971626
tau <- op$tau
nu <- op$beta - 1 #beta = nu + d/2 
result_df <- data.frame(
    parameter = c("tau", "nu", "sigma.e"),
    true = c(tau, nu, sigma.e), 
    mean = c(result_fit$summary.tau$mean,result_fit$summary.nu$mean,
             sqrt(1/rspde_fit$summary.hyperpar[1,1])),
    mode = c(result_fit$summary.tau$mode, result_fit$summary.nu$mode,
             sqrt(1/rspde_fit$summary.hyperpar[1,6]))
)
print(result_df)
#>   parameter true       mean       mode
#> 1       tau  0.2 0.12483234 0.11649557
#> 2        nu  0.8 0.96927017 0.97162580
#> 3   sigma.e  0.1 0.09790113 0.09826668
```

### Extreme value models

When used for extreme value statistics, one might want to use a
particular form of the mean value of the latent field $`u`$, which is
zero at one location $`k`$ and is given by the diagonal of
$`Q_{-k,-k}^{-1}`$ for the remaining locations. This option can be
specified via the `mean.correction` argument of `rspde.intrinsic`:

``` r

rspde_model2 <- rspde.intrinsic(mesh = mesh_2d, rspde.order = rspde.order,
                                mean.correction = TRUE)
```

We can then fit this model as before:

``` r

f <- y ~ -1 + f(field, model = rspde_model2)
rspde_fit <- inla(f,
                  data = inla.stack.data(st.dat),
                  family = "gaussian",
                  control.predictor = list(A = inla.stack.A(st.dat)))
```

To see the posterior distributions of the parameters we can do:

``` r

result_fit <- rspde.result(rspde_fit, "field", rspde_model2)
posterior_df_fit <- gg_df(result_fit)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](intrinsic_files/figure-html/unnamed-chunk-13-1.png)

### An example with replicates

Let us redo the previous example with replicated data to illustrate that
replicates are handled in the same way as any other `rSPDE` model. We
start by generating some data with 200 observations per replicate

``` r

set.seed(1)
tau <- 0.2
beta <- 1.9
op <- intrinsic.operators(tau = tau, beta = beta, mesh = mesh_2d)
n.rep <- 5
m <- 1000
loc_2d_mesh <- matrix(2*runif(m * 2), m, 2)

A <- spde.make.A(
  mesh = mesh_2d,
  loc = loc_2d_mesh,
  index = rep(1:m, times = n.rep),
  repl = rep(1:n.rep, each = m)
)

u <- simulate(op, nsim = n.rep)
y <- as.vector(A %*% as.vector(u)) +
  rnorm(m * n.rep) * 0.1
```

We now create the stack, A matrix and index and fit the model:

``` r

Abar.rep <- rspde.make.A(
  mesh = mesh_2d, loc = loc_2d_mesh, index = rep(1:m, times = n.rep),
  repl = rep(1:n.rep, each = m)
)
mesh.index.rep <- rspde.make.index(
  name = "field", mesh = mesh_2d,
  n.repl = n.rep
)

st.dat.rep <- inla.stack(
  data = list(y = y),
  A = Abar.rep,
  effects = mesh.index.rep
)

rspde_model.rep <- rspde.intrinsic(mesh = mesh_2d, prior.nu.dist = "beta")

f.rep <-
  y ~ -1 + f(field,
    model = rspde_model.rep,
    replicate = field.repl
  )
rspde_fit.rep <-
  inla(f.rep,
    data = inla.stack.data(st.dat.rep),
    family = "gaussian",
    control.predictor =
      list(A = inla.stack.A(st.dat.rep))
  )
```

We then compare with the true parameter estimates as before

``` r

result_fit <- rspde.result(rspde_fit.rep, "field", rspde_model.rep)
summary(result_fit)
#>         mean         sd 0.025quant 0.5quant 0.975quant     mode
#> tau 0.177477 0.00858324   0.161060 0.177339   0.194753 0.177147
#> nu  0.924320 0.01193950   0.901128 0.924213   0.948022 0.923887
tau <- op$tau
nu <- op$beta - 1 #beta = nu + d/2 
result_df <- data.frame(
    parameter = c("tau", "nu", "sigma.e"),
    true = c(tau, nu, sigma.e), 
    mean = c(result_fit$summary.tau$mean,result_fit$summary.nu$mean,
             sqrt(1/rspde_fit.rep$summary.hyperpar[1,1])),
    mode = c(result_fit$summary.tau$mode, result_fit$summary.nu$mode,
             sqrt(1/rspde_fit.rep$summary.hyperpar[1,6]))
)
print(result_df)
#>   parameter true       mean       mode
#> 1       tau  0.2 0.17747662 0.17714713
#> 2        nu  0.9 0.92431997 0.92388676
#> 3   sigma.e  0.1 0.09994918 0.09997177
```

To see the posterior distributions of the parameters we can do:

``` r

result_fit <- rspde.result(rspde_fit.rep, "field", rspde_model.rep)
posterior_df_fit <- gg_df(result_fit)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](intrinsic_files/figure-html/unnamed-chunk-17-1.png)

## A more general model

The `rSPDE` package also contains a partial implementation of a more
general intrinsic model, which we refer to as an intrinsic Matérn model.
The model is defined as  
``` math
(-\Delta)^{\beta/2}(\kappa^2-\Delta)^{\alpha/2}(\tau u) = \mathcal{W},
```
where $`\alpha + \beta > d/2`$ and $`d`$ is the dimension of the spatial
domain. These models are handled by performing two rational
approximations, one for each fractional operator.

To illustrate this model, we consider the same mesh as before and use
the
[`intrinsic.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.matern.operators.md)
function to construct the `rSPDE` representation of the general model.

``` r

bnd <- fm_segm(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2)), is.bnd = TRUE)
mesh_2d <- fm_mesh_2d(
    boundary = bnd, 
    cutoff = 0.01,
    max.edge = c(0.05)
)

kappa <- 10
tau <- 0.0025
alpha <- 2
beta <- 1
op <- intrinsic.matern.operators(kappa = kappa, tau = tau, alpha = alpha, 
                                 beta = beta, mesh = mesh_2d)
```

To see that the `rSPDE` model is approximating the true model, we can
compare the variogram of the approximation with the true variogram
(implemented in
[`variogram.intrinsic.spde()`](https://davidbolin.github.io/rSPDE/reference/variogram.intrinsic.spde.md))
as follows.

``` r

point <- matrix(c(1,1),1,2)
Gamma <- op$variogram(point)
vario <- variogram.intrinsic.spde(point, mesh_2d$loc[,1:2], kappa = kappa, 
                                  alpha = alpha, tau = tau,
                                  beta = beta, L = 2, d = 2)

d = sqrt((mesh_2d$loc[,1]-point[1])^2 +  (mesh_2d$loc[,2]-point[2])^2)
plot(d, Gamma, xlim = c(0,0.5), ylim = c(0,4),
     ylab = "variogram(h)", xlab = "h")
lines(sort(d),sort(vario),col=2, lwd = 2)
```

![](intrinsic_files/figure-html/unnamed-chunk-19-1.png)

We can now use the `simulate` function to simulate a realization of the
field $`u`$:

``` r

u <- simulate(op,nsim = 1, use_kl = FALSE)

proj <- fm_evaluator(mesh_2d, dims = c(100, 100))
field <- fm_evaluate(proj, field = as.vector(u))
field.df <- data.frame(x1 = proj$lattice$loc[,1],
                       x2 = proj$lattice$loc[,2], 
                       y = as.vector(field))

library(ggplot2)
library(viridis)
ggplot(field.df, aes(x = x1, y = x2, fill = y)) +
    geom_raster() +
    scale_fill_viridis()
```

![](intrinsic_files/figure-html/unnamed-chunk-20-1.png)

By default, the field is simulated with a zero-integral constraint.

### Fitting the model with `R-INLA`

We will now fit the model using our [`R-INLA`](https://www.r-inla.org)
implementation of the rational SPDE approach. Further details on this
implementation can be found in [R-INLA implementation of the rational
SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md).

We begin by simulating some data as before.

``` r

n_loc <- 2000
loc_2d_mesh <- matrix(2*runif(n_loc * 2), n_loc, 2)

A <- spde.make.A(
  mesh = mesh_2d,
  loc = loc_2d_mesh
)
sigma.e <- 0.1
y <- A %*% u + rnorm(n_loc) * sigma.e
```

The generated data can be seen in the following image.

``` r

df <- data.frame(x1 = as.double(loc_2d_mesh[, 1]),
  x2 = as.double(loc_2d_mesh[, 2]), y = as.double(y))
ggplot(df, aes(x = x1, y = x2, col = y)) +
  geom_point() +
  scale_color_viridis()
```

![](intrinsic_files/figure-html/unnamed-chunk-22-1.png)

To fit the model, we create the $`A`$ matrix, the index, and the
`inla.stack` object. For now, these more general models can only be
estimated with $`\beta = 1`$ and $`\alpha = 1`$ or $`\alpha = 2`$. For
these non-fractional models, we can use the standard INLA functions to
make the required elements.

``` r

mesh.index <- inla.spde.make.index(name = "field", n.spde = mesh_2d$n)

st.dat <- inla.stack(data = list(y = as.vector(y)), A = A, effects = mesh.index)
```

We now create the model object.

``` r

rspde_model <- rspde.intrinsic.matern(mesh = mesh_2d, alpha = alpha)
```

Finally, we create the formula and fit the model to the data:

``` r

f <- y ~ -1 + f(field, model = rspde_model)
rspde_fit <- inla(f,
                  data = inla.stack.data(st.dat),
                  family = "gaussian",
                  control.predictor = list(A = inla.stack.A(st.dat)))
```

We can get a summary of the fit:

``` r

summary(rspde_fit)
#> Time used:
#>     Pre = 0.15, Running = 21.8, Post = 0.0575, Total = 22 
#> Random effects:
#>   Name     Model
#>     field CGeneric
#> 
#> Model hyperparameters:
#>                                           mean    sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 100.68 4.491      92.10   100.59
#> Theta1 for field                         -5.98 0.048      -6.07    -5.98
#> Theta2 for field                          2.35 0.086       2.17     2.35
#>                                         0.975quant   mode
#> Precision for the Gaussian observations     109.78 100.44
#> Theta1 for field                             -5.88  -5.98
#> Theta2 for field                              2.51   2.35
#> 
#> Marginal log-Likelihood:  727.86 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```

To get a summary of the fit of the random field only, we can do the
following:

``` r

result_fit <- rspde.result(rspde_fit, "field", rspde_model)
summary(result_fit)
#>              mean          sd 0.025quant    0.5quant  0.975quant        mode
#> tau    0.00253234 0.000121303 0.00230553  0.00252725  0.00278188  0.00251618
#> kappa 10.49090000 0.897678000 8.79855000 10.47020000 12.32000000 10.44900000
tau <- op$tau
result_df <- data.frame(
  parameter = c("tau", "kappa"),
  true = c(tau, kappa), mean = c(result_fit$summary.tau$mean,
                                     result_fit$summary.kappa$mean),
  mode = c(result_fit$summary.tau$mode, result_fit$summary.kappa$mode)
)
print(result_df)
#>   parameter    true         mean         mode
#> 1       tau  0.0025  0.002532338  0.002516184
#> 2     kappa 10.0000 10.490890222 10.448969604
```

### Kriging with `R-INLA` implementation

Let us now obtain predictions (i.e., do kriging) of the latent field on
a dense grid in the region.

We begin by creating the grid of locations where we want to compute the
predictions. To this end, we can use the
[`rspde.mesh.projector()`](https://davidbolin.github.io/rSPDE/reference/rspde.mesh.project.md)
function. This function has the same arguments as the function
[`inla.mesh.projector()`](https://rdrr.io/pkg/INLA/man/inla.mesh.project.html)
the only difference being that the rSPDE version also has an argument
`nu` and an argument `rspde.order`. Thus, we proceed in the same fashion
as we would in [`R-INLA`](https://www.r-inla.org)’s standard SPDE
implementation:

``` r

projgrid <- inla.mesh.projector(mesh_2d,
  xlim = c(0, 2),
  ylim = c(0, 2)
)
#> Warning: `inla.mesh.projector()` was deprecated in INLA 23.06.07.
#> ℹ Please use `fmesher::fm_evaluator()` instead.
#> ℹ For more information, see
#>   https://inlabru-org.github.io/fmesher/articles/inla_conversion.html
#> ℹ To silence these deprecation messages in old legacy code, set
#>   `inla.setOption(fmesher.evolution.warn = FALSE)`.
#> ℹ To ensure visibility of these messages in package tests, also set
#>   `inla.setOption(fmesher.evolution.verbosity = 'warn')`.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

This lattice contains 100 × 100 locations (the default). Let us now
calculate the predictions jointly with the estimation. To this end,
first, we begin by linking the prediction coordinates to the mesh nodes
through an $`A`$ matrix

``` r

A.prd <- projgrid$proj$A
```

We now make a stack for the prediction locations. We have no data at the
prediction locations, so we set `y= NA`. We then join this stack with
the estimation stack.

``` r

ef.prd <- list(c(mesh.index))
st.prd <- inla.stack(
  data = list(y = NA),
  A = list(A.prd), tag = "prd",
  effects = ef.prd
)
st.all <- inla.stack(st.dat, st.prd)
```

Doing the joint estimation takes a while, and we therefore turn off the
computation of certain things that we are not interested in, such as the
marginals for the random effect. We will also use a simplified
integration strategy (actually only using the posterior mode of the
hyper-parameters) through the command
`control.inla = list(int.strategy = "eb")`, i.e. empirical Bayes:

``` r

rspde_fitprd <- inla(f,
  family = "Gaussian",
  data = inla.stack.data(st.all),
  control.predictor = list(
    A = inla.stack.A(st.all),
    compute = TRUE, link = 1
  ),
  control.compute = list(
    return.marginals = FALSE,
    return.marginals.predictor = FALSE
  ),
  control.inla = list(int.strategy = "eb")
)
```

We then extract the indices to the prediction nodes and then extract the
mean and the standard deviation of the response:

``` r

id.prd <- inla.stack.index(st.all, "prd")$data
m.prd <- matrix(rspde_fitprd$summary.fitted.values$mean[id.prd], 100, 100)
sd.prd <- matrix(rspde_fitprd$summary.fitted.values$sd[id.prd], 100, 100)
```

Finally, we plot the results. First the mean:

``` r

field.pred.df <- data.frame(x1 = projgrid$lattice$loc[,1],
                        x2 = projgrid$lattice$loc[,2], 
                        y = as.vector(m.prd))
ggplot(field.pred.df, aes(x = x1, y = x2, fill = y)) +
  geom_raster()  + scale_fill_viridis()
```

![](intrinsic_files/figure-html/plot_pred-1.png)

Then, the marginal standard deviations:

``` r

field.pred.sd.df <- data.frame(x1 = proj$lattice$loc[,1],
                        x2 = proj$lattice$loc[,2], 
                        sd = as.vector(sd.prd))
ggplot(field.pred.sd.df, aes(x = x1, y = x2, fill = sd)) +
  geom_raster() + scale_fill_viridis()
```

![](intrinsic_files/figure-html/plot_pred_sd-1.png)

## Using intrinsic models without `R-INLA`

Currently, the more general model is only implemented in `R-INLA` using
fixed integer values of the smoothness parameters. However, all
intrinsic models are implemented in `rSPDE` in full generality. In this
section, we illustrate the `rSPDE` interface. Let us test a model in one
dimension.

Let us start with generating the model

``` r

L = 20
x <- seq(from = 0, to = L, length.out = 101)
mesh <- fm_mesh_1d(x)
beta <- 1.1
alpha <- 0
kappa <- 10
tau <- 10
op <- intrinsic.matern.operators(kappa = kappa, tau = tau, alpha = alpha,
                                 beta = beta, mesh = mesh, d = 1)

vario <- variogram.intrinsic.spde(c(L/2), mesh$loc, tau = tau,
                                  beta = beta, alpha = alpha, kappa = kappa, L = L, d = 1)
plot(x, vario, type = "l", col = 2, lwd = 2)
points(x,op$variogram(L/2),col=1)
```

![](intrinsic_files/figure-html/unnamed-chunk-28-1.png)

We now generate some data. The option to use a mean value correction for
extremes models is also implemented, so we generate some data using
this.

``` r

n.rep <- 100
u <- simulate(op,nsim = n.rep, integral.constraint = FALSE, use_kl = TRUE)

drift <- op$mean_correction()
u <- u + matrix(rep(drift, times = n.rep), nrow = op$n, ncol= n.rep)

sigma.e <- 0.01
n.obs <- 300
obs.loc <- runif(n = n.obs, min = 0, max = L)
A <- rSPDE.A1d(x, obs.loc)
Y <- as.matrix(A %*% u + sigma.e * matrix(rnorm(n.obs*n.rep),n.obs,n.rep))
```

Let us now show how to do kriging prediction for this model.

``` r

A <- make_A(op, loc = obs.loc)
A.krig <- make_A(op, loc = x)
u.krig <- predict(op,
  A = A, Aprd = A.krig, Y = Y[,1], sigma.e = sigma.e,
  compute.variances = TRUE
)


plot(obs.loc, Y[,1],
  ylab = "u(x)", xlab = "x", main = "Data and prediction",
  ylim = c(
    min(c(min(u.krig$mean - 2 * sqrt(u.krig$variance)),min(u[,1]))),
    max(c(max(u.krig$mean + 2 * sqrt(u.krig$variance)), max(u[,1])))
  )
)
lines(x,u[,1],col=3)
lines(x, u.krig$mean)
lines(x, u.krig$mean + 2 * sqrt(u.krig$variance), col = 2)
lines(x, u.krig$mean - 2 * sqrt(u.krig$variance), col = 2)
```

![](intrinsic_files/figure-html/unnamed-chunk-30-1.png)

We now use `rspde_lme` to fit the parameters based on this data. Since
we generated data with `alpha=0`, we specify this in the function to
indicate that this parameter should not be fitted but kept fixed at
`alpha=0` by setting `fix_alpha=0` in `model_options`. We also specify
`mean_correction=TRUE` to indicate that we should use the mean value
correction when fitting.

``` r

data = data.frame(y = c(Y), loc = rep(obs.loc, n.rep), rep  = rep(1:n.rep, each = n.obs))

fit <- rspde_lme(y ~ -1, loc = "loc", repl  = "rep", data = data,
                 model = op, mean_correction = TRUE, parallel = TRUE,
                 model_options = list(fix_alpha = 0))

rbind(c(fit$coeff$random_effects[c("beta", "tau")], fit$coeff$measurement_error), 
      c(beta, tau, sigma.e))
#>          beta       tau    std. dev
#> [1,] 1.124194  9.776033 0.009983577
#> [2,] 1.100000 10.000000 0.010000000
```

### An example with estimated alpha and beta parameters

In the previous example, we fixed the alpha parameter and only estimated
beta. Now, let us demonstrate how to estimate both alpha and beta
simultaneously. We will set up a new model with different parameter
values:

``` r

L = 20
x <- seq(from = 0, to = L, length.out = 101)
mesh <- fm_mesh_1d(x)
beta <- 1.2
alpha <- 0.3
kappa <- 15
tau <- 7
op <- intrinsic.matern.operators(kappa = kappa, tau = tau, alpha = alpha,
                                 beta = beta, mesh = mesh, d = 1)

vario <- variogram.intrinsic.spde(c(L/2), mesh$loc, tau = tau,
                                  beta = beta, alpha = alpha, kappa = kappa, L = L, d = 1)
plot(x, vario, type = "l", col = 2, lwd = 2)
points(x, op$variogram(L/2), col = 1)
```

![](intrinsic_files/figure-html/unnamed-chunk-32-1.png)

We can note here that the variogram of the approximate model is not
particularly close to the variogram of the true continuous model. The
reason for this is that the value of `alpha` is very small, and we
therefore need a larger order of the rational approximation than the
default value of 2. We can adjust the orders of the rational
approximations through the `m_alpha` and `m_beta` values in
`intrinsic.matern.operators`. Let us increase the value of `m_alpha` and
decrease the value of `m_beta`.

``` r

op <- intrinsic.matern.operators(kappa = kappa, tau = tau, alpha = alpha,
                                 beta = beta, mesh = mesh, d = 1, m_alpha = 6, 
                                 m_beta = 1)

vario <- variogram.intrinsic.spde(c(L/2), mesh$loc, tau = tau,
                                  beta = beta, alpha = alpha, kappa = kappa, L = L, d = 1)
plot(x, vario, type = "l", col = 2, lwd = 2)
points(x, op$variogram(L/2), col = 1)
```

![](intrinsic_files/figure-html/unnamed-chunk-33-1.png)

We now have a better approximation. Similar to the previous example, we
will generate data with the mean value correction for extremes models:

``` r

n.rep <- 100
u <- simulate(op, nsim = n.rep, integral.constraint = FALSE, use_kl = TRUE)

drift <- op$mean_correction()
u <- u + matrix(rep(drift, times = n.rep), nrow = op$n, ncol = n.rep)

sigma.e <- 0.015
n.obs <- 300
obs.loc <- runif(n = n.obs, min = 0, max = L)
A <- rSPDE.A1d(x, obs.loc)
Y <- as.matrix(A %*% u + sigma.e * matrix(rnorm(n.obs*n.rep), n.obs, n.rep))
```

Let’s visualize the data and predictions for this model:

``` r

A <- make_A(op, loc = obs.loc)
A.krig <- make_A(op, loc = x)
u.krig <- predict(op,
  A = A, Aprd = A.krig, Y = Y[,1], sigma.e = sigma.e,
  compute.variances = TRUE
)

plot(obs.loc, Y[,1],
  ylab = "u(x)", xlab = "x", main = "Data and prediction with fractional alpha and beta",
  ylim = c(
    min(c(min(u.krig$mean - 2 * sqrt(u.krig$variance)), min(u[,1]))),
    max(c(max(u.krig$mean + 2 * sqrt(u.krig$variance)), max(u[,1])))
  )
)
lines(x, u[,1], col = 3)
lines(x, u.krig$mean)
lines(x, u.krig$mean + 2 * sqrt(u.krig$variance), col = 2)
lines(x, u.krig$mean - 2 * sqrt(u.krig$variance), col = 2)
```

![](intrinsic_files/figure-html/unnamed-chunk-35-1.png)

Now, we will use `rspde_lme` to fit the parameters but this time we will
not fix alpha, allowing both alpha and beta to be estimated. Unlike the
previous example where we set `fix_alpha=0`, we do not include this
constraint:

``` r

data = data.frame(y = c(Y), loc = rep(obs.loc, n.rep), rep = rep(1:n.rep, each = n.obs))
op <- intrinsic.matern.operators(kappa = kappa, tau = tau, alpha = 1.3, beta = 1.05, mesh = mesh, d = 1, m_alpha = 3, m_beta = 1)
fit <- rspde_lme(y ~ -1, loc = "loc", repl = "rep", data = data,
                 model = op, mean_correction = TRUE, parallel = FALSE)
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01594345 , lik =  -19712.32 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01594345 , lik =  -19712.32 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  10.58863 , beta = 1.05 , sigma_e =  0.01594345 , lik =  -87316.1 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01594345 , lik =  -114617 nz =  8 , nz.p =  7 
#> alpha =  1.96646 , tau =  7 , beta = 1.05 , sigma_e =  0.01594345 , lik =  -673790.5 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.588295 , sigma_e =  0.01594345 , lik =  11672.99 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.02411705 , lik =  3191.202 nz =  8 , nz.p =  7 
#> alpha =  0.8594124 , tau =  8.26028 , beta = 1.239042 , sigma_e =  0.01881391 , lik =  52555.98 nz =  8 , nz.p =  7 
#> alpha =  0.5681458 , tau =  8.973112 , beta = 1.345967 , sigma_e =  0.02043749 , lik =  66531.23 nz =  8 , nz.p =  7 
#> alpha =  0.9335766 , tau =  9.122897 , beta = 1.368435 , sigma_e =  0.02077864 , lik =  57317.98 nz =  8 , nz.p =  7 
#> alpha =  1.014141 , tau =  8.538353 , beta = 1.280753 , sigma_e =  0.01944726 , lik =  47675.72 nz =  8 , nz.p =  7 
#> alpha =  1.157752 , tau =  9.062685 , beta = 1.151997 , sigma_e =  0.01749221 , lik =  8955.246 nz =  8 , nz.p =  7 
#> alpha =  0.7807311 , tau =  9.530452 , beta = 1.578857 , sigma_e =  0.02397375 , lik =  60314.43 nz =  8 , nz.p =  7 
#> alpha =  0.8868741 , tau =  8.822865 , beta = 1.425786 , sigma_e =  0.02164948 , lik =  55448.65 nz =  8 , nz.p =  7 
#> alpha =  0.6366852 , tau =  10.78251 , beta = 1.858678 , sigma_e =  0.01581098 , lik =  61375.88 nz =  8 , nz.p =  7 
#> alpha =  0.7610784 , tau =  9.678638 , beta = 1.611389 , sigma_e =  0.01757115 , lik =  59135.78 nz =  8 , nz.p =  7 
#> alpha =  0.5628285 , tau =  8.92783 , beta = 2.051366 , sigma_e =  0.02093554 , lik =  65711.28 nz =  8 , nz.p =  7 
#> alpha =  0.6740405 , tau =  8.961355 , beta = 1.775805 , sigma_e =  0.02001586 , lik =  65066.58 nz =  8 , nz.p =  7 
#> alpha =  0.3586093 , tau =  12.73984 , beta = 1.789585 , sigma_e =  0.02561338 , lik =  67465.53 nz =  8 , nz.p =  7 
#> alpha =  0.1883478 , tau =  17.18687 , beta = 1.990638 , sigma_e =  0.03246454 , lik =  63971.26 nz =  8 , nz.p =  7 
#> alpha =  0.3405678 , tau =  11.17228 , beta = 2.2202 , sigma_e =  0.02136841 , lik =  66395.13 nz =  8 , nz.p =  7 
#> alpha =  0.4382179 , tau =  10.62037 , beta = 1.92206 , sigma_e =  0.02121941 , lik =  67617.75 nz =  8 , nz.p =  7 
#> alpha =  0.3232417 , tau =  11.1679 , beta = 2.083994 , sigma_e =  0.01763813 , lik =  67027.04 nz =  8 , nz.p =  7 
#> alpha =  0.4029684 , tau =  10.73388 , beta = 1.916267 , sigma_e =  0.01904466 , lik =  68707.49 nz =  8 , nz.p =  7 
#> alpha =  0.3300978 , tau =  9.853262 , beta = 1.776166 , sigma_e =  0.0288108 , lik =  66673.59 nz =  8 , nz.p =  7 
#> alpha =  0.3890122 , tau =  10.07778 , beta = 1.776944 , sigma_e =  0.0247974 , lik =  68172.67 nz =  8 , nz.p =  7 
#> alpha =  0.3221263 , tau =  12.49112 , beta = 1.517014 , sigma_e =  0.02328299 , lik =  68669.69 nz =  8 , nz.p =  7 
#> alpha =  0.3703511 , tau =  11.48517 , beta = 1.619457 , sigma_e =  0.02267254 , lik =  68708.92 nz =  8 , nz.p =  7 
#> alpha =  0.2689081 , tau =  13.71772 , beta = 2.356651 , sigma_e =  0.02486134 , lik =  66018.96 nz =  8 , nz.p =  7 
#> alpha =  0.471244 , tau =  9.97764 , beta = 1.537599 , sigma_e =  0.02146356 , lik =  68171.25 nz =  8 , nz.p =  7 
#> alpha =  0.475217 , tau =  8.76207 , beta = 1.700318 , sigma_e =  0.01848327 , lik =  68719.72 nz =  8 , nz.p =  7 
#> alpha =  0.5470503 , tau =  7.266547 , beta = 1.689376 , sigma_e =  0.01570127 , lik =  66199.06 nz =  8 , nz.p =  7 
#> alpha =  0.4017021 , tau =  9.732359 , beta = 1.521985 , sigma_e =  0.02111206 , lik =  69715.16 nz =  8 , nz.p =  7 
#> alpha =  0.3846016 , tau =  9.3166 , beta = 1.360867 , sigma_e =  0.02105859 , lik =  70386.14 nz =  8 , nz.p =  7 
#> alpha =  0.3443931 , tau =  10.0791 , beta = 1.781244 , sigma_e =  0.0207135 , lik =  70240.17 nz =  8 , nz.p =  7 
#> alpha =  0.3724795 , tau =  10.05364 , beta = 1.723138 , sigma_e =  0.02089852 , lik =  69910.29 nz =  8 , nz.p =  7 
#> alpha =  0.3973504 , tau =  9.979462 , beta = 1.560618 , sigma_e =  0.01668353 , lik =  69971.73 nz =  8 , nz.p =  7 
#> alpha =  0.3952493 , tau =  10.00395 , beta = 1.612145 , sigma_e =  0.01842118 , lik =  70092.89 nz =  8 , nz.p =  7 
#> alpha =  0.3806326 , tau =  9.109217 , beta = 1.357759 , sigma_e =  0.02143574 , lik =  70016.66 nz =  8 , nz.p =  7 
#> alpha =  0.3860977 , tau =  9.490739 , beta = 1.477279 , sigma_e =  0.0208112 , lik =  69965.72 nz =  8 , nz.p =  7 
#> alpha =  0.4187166 , tau =  7.759573 , beta = 1.485738 , sigma_e =  0.01760625 , lik =  71145.59 nz =  8 , nz.p =  7 
#> alpha =  0.4452188 , tau =  6.37805 , beta = 1.418411 , sigma_e =  0.01551494 , lik =  71617.25 nz =  8 , nz.p =  7 
#> alpha =  0.3179175 , tau =  8.957795 , beta = 1.332092 , sigma_e =  0.02013616 , lik =  71681.14 nz =  8 , nz.p =  7 
#> alpha =  0.2600313 , tau =  9.057291 , beta = 1.192711 , sigma_e =  0.02101723 , lik =  71759.55 nz =  8 , nz.p =  7 
#> alpha =  0.3408069 , tau =  8.596738 , beta = 1.560212 , sigma_e =  0.0172263 , lik =  72328.65 nz =  8 , nz.p =  7 
#> alpha =  0.322485 , tau =  8.351414 , beta = 1.667689 , sigma_e =  0.01544254 , lik =  72763.62 nz =  8 , nz.p =  7 
#> alpha =  0.302551 , tau =  7.282467 , beta = 1.343171 , sigma_e =  0.01867985 , lik =  72851.33 nz =  8 , nz.p =  7 
#> alpha =  0.2647052 , tau =  6.213439 , beta = 1.234156 , sigma_e =  0.01881055 , lik =  72744.7 nz =  8 , nz.p =  7 
#> alpha =  0.3296954 , tau =  6.347041 , beta = 1.093949 , sigma_e =  0.01593989 , lik =  73739.58 nz =  8 , nz.p =  7 
#> alpha =  0.3225834 , tau =  5.036701 , beta = 0.8736344 , sigma_e =  0.01398304 , lik =  73086.99 nz =  0 , nz.p =  0 
#> alpha =  0.2775856 , tau =  5.889264 , beta = 1.29212 , sigma_e =  0.01402594 , lik =  73706.24 nz =  8 , nz.p =  7 
#> alpha =  0.3011622 , tau =  6.604807 , beta = 1.31001 , sigma_e =  0.01552588 , lik =  73756.12 nz =  8 , nz.p =  7 
#> alpha =  0.205081 , tau =  8.723673 , beta = 1.19511 , sigma_e =  0.01904292 , lik =  72726.83 nz =  8 , nz.p =  7 
#> alpha =  0.2489358 , tau =  8.066711 , beta = 1.249741 , sigma_e =  0.01809205 , lik =  73179.41 nz =  8 , nz.p =  7 
#> alpha =  0.2790896 , tau =  8.124859 , beta = 1.253174 , sigma_e =  0.01872405 , lik =  72856.18 nz =  8 , nz.p =  7 
#> alpha =  0.2626079 , tau =  6.291134 , beta = 0.9621539 , sigma_e =  0.01946267 , lik =  71909.36 nz =  0 , nz.p =  0 
#> alpha =  0.3063438 , tau =  7.780406 , beta = 1.438104 , sigma_e =  0.01636212 , lik =  73320.33 nz =  8 , nz.p =  7 
#> alpha =  0.2813085 , tau =  7.407981 , beta = 1.193035 , sigma_e =  0.01525928 , lik =  73822.31 nz =  8 , nz.p =  7 
#> alpha =  0.2712532 , tau =  7.471546 , beta = 1.127269 , sigma_e =  0.0137916 , lik =  73667.04 nz =  8 , nz.p =  7 
#> alpha =  0.3059239 , tau =  6.39906 , beta = 1.25082 , sigma_e =  0.01402689 , lik =  73714.26 nz =  8 , nz.p =  7 
#> alpha =  0.2989826 , tau =  6.792679 , beta = 1.251587 , sigma_e =  0.01507721 , lik =  73815.56 nz =  8 , nz.p =  7 
#> alpha =  0.3690567 , tau =  6.016907 , beta = 1.245584 , sigma_e =  0.01349606 , lik =  73240.41 nz =  8 , nz.p =  7 
#> alpha =  0.3344581 , tau =  6.474469 , beta = 1.249986 , sigma_e =  0.01452202 , lik =  73667.68 nz =  8 , nz.p =  7 
#> alpha =  0.3106164 , tau =  5.796267 , beta = 1.035978 , sigma_e =  0.01422755 , lik =  73602.25 nz =  8 , nz.p =  7 
#> alpha =  0.3095427 , tau =  6.238955 , beta = 1.122639 , sigma_e =  0.01473355 , lik =  73814.57 nz =  8 , nz.p =  7 
#> alpha =  0.2758312 , tau =  6.86293 , beta = 1.138624 , sigma_e =  0.01612332 , lik =  73696.8 nz =  8 , nz.p =  7 
#> alpha =  0.2894463 , tau =  6.763682 , beta = 1.164829 , sigma_e =  0.01570716 , lik =  73794.28 nz =  8 , nz.p =  7 
#> alpha =  0.265616 , tau =  7.181037 , beta = 1.322142 , sigma_e =  0.01460297 , lik =  73775.07 nz =  8 , nz.p =  7 
#> alpha =  0.2803619 , tau =  6.962788 , beta = 1.263818 , sigma_e =  0.0149263 , lik =  73819.73 nz =  8 , nz.p =  7 
#> alpha =  0.2825752 , tau =  7.048027 , beta = 1.099394 , sigma_e =  0.01475804 , lik =  73795.05 nz =  8 , nz.p =  7 
#> alpha =  0.2871115 , tau =  6.934509 , beta = 1.147294 , sigma_e =  0.01494637 , lik =  73833.16 nz =  8 , nz.p =  7 
#> alpha =  0.293062 , tau =  6.951483 , beta = 1.225309 , sigma_e =  0.01430087 , lik =  73756.74 nz =  8 , nz.p =  7 
#> alpha =  0.290346 , tau =  6.810151 , beta = 1.179564 , sigma_e =  0.01534313 , lik =  73830.99 nz =  8 , nz.p =  7 
#> alpha =  0.2671069 , tau =  7.80489 , beta = 1.292076 , sigma_e =  0.01549514 , lik =  73765.94 nz =  8 , nz.p =  7 
#> alpha =  0.2983402 , tau =  6.598202 , beta = 1.164088 , sigma_e =  0.01492036 , lik =  73837.26 nz =  8 , nz.p =  7 
#> alpha =  0.2763033 , tau =  7.08584 , beta = 1.13115 , sigma_e =  0.01507874 , lik =  73825.72 nz =  8 , nz.p =  7 
#> alpha =  0.2818065 , tau =  7.011384 , beta = 1.159486 , sigma_e =  0.01507836 , lik =  73837.11 nz =  8 , nz.p =  7 
#> alpha =  0.2952609 , tau =  6.931983 , beta = 1.080637 , sigma_e =  0.01529307 , lik =  73790.16 nz =  8 , nz.p =  7 
#> alpha =  0.2840146 , tau =  6.955074 , beta = 1.215283 , sigma_e =  0.01501716 , lik =  73836.12 nz =  8 , nz.p =  7 
#> alpha =  0.2953964 , tau =  6.35304 , beta = 1.152825 , sigma_e =  0.01486395 , lik =  73820.08 nz =  8 , nz.p =  7 
#> alpha =  0.2847662 , tau =  7.128863 , beta = 1.18298 , sigma_e =  0.01515947 , lik =  73835.38 nz =  8 , nz.p =  7 
#> alpha =  0.2839887 , tau =  7.038337 , beta = 1.167716 , sigma_e =  0.01471169 , lik =  73828.52 nz =  8 , nz.p =  7 
#> alpha =  0.2887435 , tau =  6.866495 , beta = 1.176583 , sigma_e =  0.01518277 , lik =  73837.56 nz =  8 , nz.p =  7 
#> alpha =  0.2878399 , tau =  6.88492 , beta = 1.212976 , sigma_e =  0.01519732 , lik =  73831.29 nz =  8 , nz.p =  7 
#> alpha =  0.2872934 , tau =  6.922078 , beta = 1.163289 , sigma_e =  0.01500871 , lik =  73838.16 nz =  8 , nz.p =  7 
#> alpha =  0.2912381 , tau =  6.618823 , beta = 1.168205 , sigma_e =  0.01492389 , lik =  73835.87 nz =  8 , nz.p =  7 
#> alpha =  0.2896064 , tau =  6.742806 , beta = 1.171906 , sigma_e =  0.01498244 , lik =  73838.75 nz =  8 , nz.p =  7 
#> alpha =  0.294295 , tau =  6.700618 , beta = 1.120639 , sigma_e =  0.01505138 , lik =  73820.22 nz =  8 , nz.p =  7 
#> alpha =  0.2865505 , tau =  6.890568 , beta = 1.190954 , sigma_e =  0.01502571 , lik =  73839.83 nz =  8 , nz.p =  7 
#> alpha =  0.2985882 , tau =  6.600751 , beta = 1.187505 , sigma_e =  0.01496934 , lik =  73839.8 nz =  8 , nz.p =  7 
#> alpha =  0.2943013 , tau =  6.701098 , beta = 1.180396 , sigma_e =  0.01499652 , lik =  73839.9 nz =  8 , nz.p =  7 
#> alpha =  0.2805072 , tau =  7.057644 , beta = 1.188829 , sigma_e =  0.01515869 , lik =  73838.92 nz =  8 , nz.p =  7 
#> alpha =  0.284863 , tau =  6.939868 , beta = 1.182742 , sigma_e =  0.01509875 , lik =  73839.76 nz =  8 , nz.p =  7 
#> alpha =  0.2882658 , tau =  6.810777 , beta = 1.179072 , sigma_e =  0.01486367 , lik =  73837.27 nz =  8 , nz.p =  7 
#> alpha =  0.288624 , tau =  6.852523 , beta = 1.177205 , sigma_e =  0.01510236 , lik =  73839.32 nz =  8 , nz.p =  7 
#> alpha =  0.2902567 , tau =  6.72885 , beta = 1.198357 , sigma_e =  0.0150735 , lik =  73840.13 nz =  8 , nz.p =  7 
#> alpha =  0.2917498 , tau =  6.634268 , beta = 1.216446 , sigma_e =  0.015106 , lik =  73838.49 nz =  8 , nz.p =  7 
#> alpha =  0.288197 , tau =  6.902028 , beta = 1.200107 , sigma_e =  0.01513657 , lik =  73837.76 nz =  8 , nz.p =  7 
#> alpha =  0.2892534 , tau =  6.782264 , beta = 1.178888 , sigma_e =  0.01502083 , lik =  73839.93 nz =  8 , nz.p =  7 
#> alpha =  0.2894301 , tau =  6.76358 , beta = 1.195398 , sigma_e =  0.01498391 , lik =  73840.02 nz =  8 , nz.p =  7 
#> alpha =  0.2892283 , tau =  6.785707 , beta = 1.190815 , sigma_e =  0.01501343 , lik =  73840.28 nz =  8 , nz.p =  7 
#> alpha =  0.2950408 , tau =  6.618712 , beta = 1.193015 , sigma_e =  0.01495355 , lik =  73839.52 nz =  8 , nz.p =  7 
#> alpha =  0.287374 , tau =  6.858146 , beta = 1.185304 , sigma_e =  0.01506232 , lik =  73840.17 nz =  8 , nz.p =  7 
#> alpha =  0.2936401 , tau =  6.653502 , beta = 1.182475 , sigma_e =  0.01504088 , lik =  73840.15 nz =  8 , nz.p =  7 
#> alpha =  0.2918514 , tau =  6.711993 , beta = 1.184611 , sigma_e =  0.01503708 , lik =  73840.31 nz =  8 , nz.p =  7 
#> alpha =  0.2849523 , tau =  6.846074 , beta = 1.194683 , sigma_e =  0.01508645 , lik =  73840.29 nz =  8 , nz.p =  7 
#> alpha =  0.2872613 , tau =  6.809538 , beta = 1.191137 , sigma_e =  0.01506391 , lik =  73840.4 nz =  8 , nz.p =  7 
#> alpha =  0.2891248 , tau =  6.775012 , beta = 1.20131 , sigma_e =  0.0150793 , lik =  73839.84 nz =  8 , nz.p =  7 
#> alpha =  0.2892213 , tau =  6.78045 , beta = 1.184446 , sigma_e =  0.01503542 , lik =  73840.29 nz =  8 , nz.p =  7 
#> alpha =  0.2877138 , tau =  6.849691 , beta = 1.176319 , sigma_e =  0.01501141 , lik =  73839.71 nz =  8 , nz.p =  7 
#> alpha =  0.2896189 , tau =  6.758859 , beta = 1.192791 , sigma_e =  0.01505795 , lik =  73840.36 nz =  8 , nz.p =  7 
#> alpha =  0.2915059 , tau =  6.681465 , beta = 1.192231 , sigma_e =  0.01502081 , lik =  73840.52 nz =  8 , nz.p =  7 
#> alpha =  0.293594 , tau =  6.594839 , beta = 1.19572 , sigma_e =  0.0150001 , lik =  73840.49 nz =  8 , nz.p =  7 
#> alpha =  0.2905471 , tau =  6.711105 , beta = 1.187266 , sigma_e =  0.01507268 , lik =  73840.29 nz =  8 , nz.p =  7 
#> alpha =  0.2902168 , tau =  6.729678 , beta = 1.188154 , sigma_e =  0.01505785 , lik =  73840.41 nz =  8 , nz.p =  7 
#> alpha =  0.2909538 , tau =  6.696145 , beta = 1.195162 , sigma_e =  0.01505961 , lik =  73840.49 nz =  8 , nz.p =  7 
#> alpha =  0.2905197 , tau =  6.717123 , beta = 1.192469 , sigma_e =  0.01505356 , lik =  73840.49 nz =  8 , nz.p =  7 
#> alpha =  0.2878048 , tau =  6.766509 , beta = 1.198112 , sigma_e =  0.01506455 , lik =  73840.34 nz =  8 , nz.p =  7 
#> alpha =  0.2888112 , tau =  6.752838 , beta = 1.194734 , sigma_e =  0.01505768 , lik =  73840.44 nz =  8 , nz.p =  7 
#> alpha =  0.2896995 , tau =  6.717195 , beta = 1.190701 , sigma_e =  0.01504356 , lik =  73840.55 nz =  8 , nz.p =  7 
#> alpha =  0.2897398 , tau =  6.696459 , beta = 1.189657 , sigma_e =  0.01503637 , lik =  73840.54 nz =  8 , nz.p =  7 
#> alpha =  0.2906013 , tau =  6.699306 , beta = 1.191466 , sigma_e =  0.01503218 , lik =  73840.56 nz =  8 , nz.p =  7 
#> alpha =  0.2899581 , tau =  6.723434 , beta = 1.189427 , sigma_e =  0.0150507 , lik =  73840.49 nz =  8 , nz.p =  7 
#> alpha =  0.289255 , tau =  6.734993 , beta = 1.192716 , sigma_e =  0.01505062 , lik =  73840.53 nz =  8 , nz.p =  7 
#> alpha =  0.2901093 , tau =  6.717159 , beta = 1.191584 , sigma_e =  0.01504856 , lik =  73840.53 nz =  8 , nz.p =  7 
#> alpha =  0.291378 , tau =  6.673902 , beta = 1.191432 , sigma_e =  0.01503651 , lik =  73840.57 nz =  8 , nz.p =  7 
#> alpha =  0.292839 , tau =  6.629692 , beta = 1.191681 , sigma_e =  0.01502791 , lik =  73840.58 nz =  8 , nz.p =  7 
#> alpha =  0.2910392 , tau =  6.675787 , beta = 1.193843 , sigma_e =  0.01503043 , lik =  73840.59 nz =  8 , nz.p =  7 
#> alpha =  0.2915813 , tau =  6.65209 , beta = 1.196061 , sigma_e =  0.0150203 , lik =  73840.58 nz =  8 , nz.p =  7 
#> alpha =  0.2924651 , tau =  6.640833 , beta = 1.190987 , sigma_e =  0.01502245 , lik =  73840.56 nz =  8 , nz.p =  7 
#> alpha =  0.2916593 , tau =  6.664249 , beta = 1.191422 , sigma_e =  0.01502949 , lik =  73840.58 nz =  8 , nz.p =  7 
#> alpha =  0.292226 , tau =  6.637435 , beta = 1.19206 , sigma_e =  0.01501688 , lik =  73840.57 nz =  8 , nz.p =  7 
#> alpha =  0.2916954 , tau =  6.657277 , beta = 1.191942 , sigma_e =  0.01502479 , lik =  73840.58 nz =  8 , nz.p =  7 
#> alpha =  0.2934442 , tau =  6.613653 , beta = 1.19344 , sigma_e =  0.01501438 , lik =  73840.56 nz =  8 , nz.p =  7 
#> alpha =  0.2925035 , tau =  6.639388 , beta = 1.192756 , sigma_e =  0.01502167 , lik =  73840.58 nz =  8 , nz.p =  7 
#> alpha =  0.2932981 , tau =  6.607526 , beta = 1.193191 , sigma_e =  0.01502154 , lik =  73840.61 nz =  8 , nz.p =  7 
#> alpha =  0.2946558 , tau =  6.562108 , beta = 1.194052 , sigma_e =  0.01501622 , lik =  73840.62 nz =  8 , nz.p =  7 
#> alpha =  0.2934315 , tau =  6.601374 , beta = 1.194292 , sigma_e =  0.01501892 , lik =  73840.6 nz =  8 , nz.p =  7 
#> alpha =  0.2929874 , tau =  6.617037 , beta = 1.193574 , sigma_e =  0.01502156 , lik =  73840.6 nz =  8 , nz.p =  7 
#> alpha =  0.2929554 , tau =  6.610891 , beta = 1.19357 , sigma_e =  0.01502564 , lik =  73840.63 nz =  8 , nz.p =  7 
#> alpha =  0.2931817 , tau =  6.596689 , beta = 1.193977 , sigma_e =  0.01502763 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2943649 , tau =  6.569059 , beta = 1.195201 , sigma_e =  0.01502365 , lik =  73840.64 nz =  8 , nz.p =  7 
#> alpha =  0.2936952 , tau =  6.591004 , beta = 1.194385 , sigma_e =  0.01502394 , lik =  73840.63 nz =  8 , nz.p =  7 
#> alpha =  0.2938255 , tau =  6.572194 , beta = 1.196877 , sigma_e =  0.01501883 , lik =  73840.64 nz =  8 , nz.p =  7 
#> alpha =  0.2935786 , tau =  6.586522 , beta = 1.195575 , sigma_e =  0.0150211 , lik =  73840.64 nz =  8 , nz.p =  7 
#> alpha =  0.2967714 , tau =  6.486112 , beta = 1.195904 , sigma_e =  0.01501168 , lik =  73840.64 nz =  8 , nz.p =  7 
#> alpha =  0.2953279 , tau =  6.53302 , beta = 1.195393 , sigma_e =  0.01501637 , lik =  73840.64 nz =  8 , nz.p =  7 
#> alpha =  0.2951114 , tau =  6.531973 , beta = 1.195908 , sigma_e =  0.01502215 , lik =  73840.67 nz =  8 , nz.p =  7 
#> alpha =  0.295955 , tau =  6.497547 , beta = 1.196717 , sigma_e =  0.01502377 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2944028 , tau =  6.545123 , beta = 1.197216 , sigma_e =  0.01502788 , lik =  73840.69 nz =  8 , nz.p =  7 
#> alpha =  0.2944661 , tau =  6.549365 , beta = 1.196424 , sigma_e =  0.01502497 , lik =  73840.68 nz =  8 , nz.p =  7 
#> alpha =  0.2954666 , tau =  6.524298 , beta = 1.194522 , sigma_e =  0.01502889 , lik =  73840.67 nz =  8 , nz.p =  7 
#> alpha =  0.2950555 , tau =  6.536239 , beta = 1.195112 , sigma_e =  0.01502638 , lik =  73840.67 nz =  8 , nz.p =  7 
#> alpha =  0.2938551 , tau =  6.564713 , beta = 1.195894 , sigma_e =  0.01503537 , lik =  73840.69 nz =  8 , nz.p =  7 
#> alpha =  0.2942226 , tau =  6.556775 , beta = 1.195769 , sigma_e =  0.01503061 , lik =  73840.68 nz =  8 , nz.p =  7 
#> alpha =  0.2946121 , tau =  6.526969 , beta = 1.196365 , sigma_e =  0.01503276 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2947358 , tau =  6.506025 , beta = 1.196948 , sigma_e =  0.01503732 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2963776 , tau =  6.472067 , beta = 1.198552 , sigma_e =  0.01503083 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2950224 , tau =  6.506195 , beta = 1.198789 , sigma_e =  0.01503387 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2950059 , tau =  6.491225 , beta = 1.200633 , sigma_e =  0.01503761 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2966911 , tau =  6.44887 , beta = 1.199903 , sigma_e =  0.01502577 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2959795 , tau =  6.477637 , beta = 1.198899 , sigma_e =  0.01502817 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2970577 , tau =  6.429954 , beta = 1.19965 , sigma_e =  0.01503241 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2963917 , tau =  6.458555 , beta = 1.199041 , sigma_e =  0.01503128 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.295674 , tau =  6.46146 , beta = 1.201083 , sigma_e =  0.01503953 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2974778 , tau =  6.41 , beta = 1.202898 , sigma_e =  0.01503167 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2967588 , tau =  6.439044 , beta = 1.201259 , sigma_e =  0.01503195 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73545.41 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73547.61 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73543.2 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73545.89 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73544.92 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73546.28 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73544.16 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73546.96 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73543.85 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73544.91 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73545.87 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73801.02 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73801.88 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73800.16 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73801.25 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73800.79 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73800.47 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73801.62 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73800.42 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73800.77 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73801.23 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73836.15 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73836.45 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73835.85 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73836.23 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73836.07 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73836.35 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73835.94 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73836.36 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73835.94 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73836.05 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73836.21 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.23 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.33 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.13 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.26 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.2 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.3 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.16 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.3 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.15 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.18 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.24 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.68 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.64 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.69 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.67 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.66 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.67 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.483416 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.476936 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.476936 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.457535 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.463995 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.463995 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.476936 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.463995 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2963363 , tau =  6.470463 , beta = 1.199398 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.20069 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.198699 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.476936 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.463995 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2951533 , tau =  6.470463 , beta = 1.200581 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.201282 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.20069 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.201282 , sigma_e =  0.01503559 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.201984 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01505063 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.198699 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198001 , sigma_e =  0.01503559 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01505063 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2954486 , tau =  6.470463 , beta = 1.200286 , sigma_e =  0.01505063 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01505063 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01505063 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01506569 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.476936 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.7 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.463995 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2960401 , tau =  6.470463 , beta = 1.199694 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.200986 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.198995 , sigma_e =  0.01502056 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01500555 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  9.845171 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  72902.32 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.01503559 , lik =  73768.55 nz =  8 , nz.p =  7 
#> alpha =  0.4499914 , tau =  6.470463 , beta = 1.045743 , sigma_e =  0.01503559 , lik =  73056.28 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.719321 , sigma_e =  0.01503559 , lik =  73436.95 nz =  8 , nz.p =  7 
#> alpha =  0.2957442 , tau =  6.470463 , beta = 1.19999 , sigma_e =  0.0228775 , lik =  70745.28 nz =  8 , nz.p =  7 
#> alpha =  0.3498094 , tau =  7.653331 , beta = 1.327956 , sigma_e =  0.009881721 , lik =  66223.73 nz =  8 , nz.p =  7 
#> alpha =  0.3084218 , tau =  6.74783 , beta = 1.229996 , sigma_e =  0.01854661 , lik =  72988.59 nz =  8 , nz.p =  7 
#> alpha =  0.355732 , tau =  4.32453 , beta = 1.341974 , sigma_e =  0.01635227 , lik =  73638.11 nz =  8 , nz.p =  7 
#> alpha =  0.339681 , tau =  5.312024 , beta = 1.303983 , sigma_e =  0.01601267 , lik =  73745.98 nz =  8 , nz.p =  7 
#> alpha =  0.3545395 , tau =  5.733719 , beta = 1.339151 , sigma_e =  0.01250011 , lik =  72578.68 nz =  8 , nz.p =  7 
#> alpha =  0.3193559 , tau =  6.478615 , beta = 1.255876 , sigma_e =  0.01680455 , lik =  73546.05 nz =  8 , nz.p =  7 
#> alpha =  0.2118526 , tau =  5.982517 , beta = 1.571849 , sigma_e =  0.01612057 , lik =  73554.35 nz =  8 , nz.p =  7 
#> alpha =  0.2557564 , tau =  6.10094 , beta = 1.448958 , sigma_e =  0.0158422 , lik =  73730.79 nz =  8 , nz.p =  7 
#> alpha =  0.3041486 , tau =  5.843439 , beta = 0.9656147 , sigma_e =  0.01646108 , lik =  73326.91 nz =  0 , nz.p =  0 
#> alpha =  0.2978233 , tau =  6.307665 , beta = 1.481301 , sigma_e =  0.01537995 , lik =  73676.91 nz =  8 , nz.p =  7 
#> alpha =  0.2739054 , tau =  5.773993 , beta = 1.388749 , sigma_e =  0.01421554 , lik =  73728.24 nz =  8 , nz.p =  7 
#> alpha =  0.2846225 , tau =  5.942617 , beta = 1.355532 , sigma_e =  0.01482276 , lik =  73802.79 nz =  8 , nz.p =  7 
#> alpha =  0.2884305 , tau =  5.790732 , beta = 1.144953 , sigma_e =  0.0153047 , lik =  73782.17 nz =  8 , nz.p =  7 
#> alpha =  0.2907505 , tau =  5.915851 , beta = 1.219138 , sigma_e =  0.01532348 , lik =  73828.63 nz =  8 , nz.p =  7 
#> alpha =  0.3535501 , tau =  5.913621 , beta = 1.069967 , sigma_e =  0.01466159 , lik =  73765.23 nz =  8 , nz.p =  7 
#> alpha =  0.326058 , tau =  5.959904 , beta = 1.160914 , sigma_e =  0.01494823 , lik =  73824.59 nz =  8 , nz.p =  7 
#> alpha =  0.2618759 , tau =  7.11184 , beta = 1.155259 , sigma_e =  0.01411182 , lik =  73736.35 nz =  8 , nz.p =  7 
#> alpha =  0.3182932 , tau =  5.714003 , beta = 1.263968 , sigma_e =  0.01551471 , lik =  73811.23 nz =  8 , nz.p =  7 
#> alpha =  0.3097526 , tau =  5.555315 , beta = 1.279337 , sigma_e =  0.01521863 , lik =  73796.01 nz =  8 , nz.p =  7 
#> alpha =  0.3061895 , tau =  5.77119 , beta = 1.258771 , sigma_e =  0.01517266 , lik =  73829.47 nz =  8 , nz.p =  7 
#> alpha =  0.3313985 , tau =  5.97835 , beta = 1.094133 , sigma_e =  0.01558189 , lik =  73778.07 nz =  8 , nz.p =  7 
#> alpha =  0.2956579 , tau =  5.95153 , beta = 1.286574 , sigma_e =  0.015009 , lik =  73829.02 nz =  8 , nz.p =  7 
#> alpha =  0.2877248 , tau =  6.319662 , beta = 1.187306 , sigma_e =  0.01469092 , lik =  73825.9 nz =  8 , nz.p =  7 
#> alpha =  0.29508 , tau =  6.16248 , beta = 1.205719 , sigma_e =  0.01489267 , lik =  73838.01 nz =  8 , nz.p =  7 
#> alpha =  0.2698781 , tau =  6.140469 , beta = 1.305255 , sigma_e =  0.01522495 , lik =  73805.69 nz =  8 , nz.p =  7 
#> alpha =  0.3110019 , tau =  6.004542 , beta = 1.197308 , sigma_e =  0.01501693 , lik =  73838.85 nz =  8 , nz.p =  7 
#> alpha =  0.3109131 , tau =  6.223127 , beta = 1.239166 , sigma_e =  0.01473255 , lik =  73825.4 nz =  8 , nz.p =  7 
#> alpha =  0.2956652 , tau =  5.991218 , beta = 1.224124 , sigma_e =  0.01517356 , lik =  73837.23 nz =  8 , nz.p =  7 
#> alpha =  0.3057555 , tau =  6.202286 , beta = 1.15111 , sigma_e =  0.01510698 , lik =  73826.89 nz =  8 , nz.p =  7 
#> alpha =  0.2981506 , tau =  6.013252 , beta = 1.251278 , sigma_e =  0.01503344 , lik =  73837.01 nz =  8 , nz.p =  7 
#> alpha =  0.3026082 , tau =  5.945826 , beta = 1.236885 , sigma_e =  0.01510125 , lik =  73837.27 nz =  8 , nz.p =  7 
#> alpha =  0.3017757 , tau =  6.212262 , beta = 1.175492 , sigma_e =  0.01505399 , lik =  73837.79 nz =  8 , nz.p =  7 
#> alpha =  0.3008653 , tau =  6.1619 , beta = 1.193964 , sigma_e =  0.01504885 , lik =  73839.9 nz =  8 , nz.p =  7 
#> alpha =  0.3064422 , tau =  6.305564 , beta = 1.18931 , sigma_e =  0.01486582 , lik =  73837.7 nz =  8 , nz.p =  7 
#> alpha =  0.3037117 , tau =  6.225464 , beta = 1.197996 , sigma_e =  0.01494216 , lik =  73839.88 nz =  8 , nz.p =  7 
#> alpha =  0.2998475 , tau =  6.471585 , beta = 1.162668 , sigma_e =  0.01487385 , lik =  73834.93 nz =  8 , nz.p =  7 
#> alpha =  0.3019156 , tau =  6.073119 , beta = 1.217774 , sigma_e =  0.01504408 , lik =  73839.93 nz =  8 , nz.p =  7 
#> alpha =  0.3103271 , tau =  6.207688 , beta = 1.196909 , sigma_e =  0.01514332 , lik =  73839.03 nz =  8 , nz.p =  7 
#> alpha =  0.306443 , tau =  6.196355 , beta = 1.19918 , sigma_e =  0.01508026 , lik =  73840.39 nz =  8 , nz.p =  7 
#> alpha =  0.2927056 , tau =  6.451609 , beta = 1.205961 , sigma_e =  0.01504331 , lik =  73838.28 nz =  8 , nz.p =  7 
#> alpha =  0.3063233 , tau =  6.113317 , beta = 1.199567 , sigma_e =  0.01502352 , lik =  73840.4 nz =  8 , nz.p =  7 
#> alpha =  0.30076 , tau =  6.177578 , beta = 1.206166 , sigma_e =  0.01515146 , lik =  73839.35 nz =  8 , nz.p =  7 
#> alpha =  0.302971 , tau =  6.213458 , beta = 1.200039 , sigma_e =  0.01499421 , lik =  73840.55 nz =  8 , nz.p =  7 
#> alpha =  0.3044535 , tau =  6.26214 , beta = 1.212761 , sigma_e =  0.01502218 , lik =  73839.74 nz =  8 , nz.p =  7 
#> alpha =  0.3017584 , tau =  6.186809 , beta = 1.19862 , sigma_e =  0.01504218 , lik =  73840.53 nz =  8 , nz.p =  7 
#> alpha =  0.3022689 , tau =  6.15348 , beta = 1.208598 , sigma_e =  0.0150396 , lik =  73840.56 nz =  8 , nz.p =  7 
#> alpha =  0.2972148 , tau =  6.256291 , beta = 1.203494 , sigma_e =  0.01497395 , lik =  73840.14 nz =  8 , nz.p =  7 
#> alpha =  0.3010461 , tau =  6.331926 , beta = 1.19962 , sigma_e =  0.01505791 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2989888 , tau =  6.309981 , beta = 1.204283 , sigma_e =  0.0150376 , lik =  73840.69 nz =  8 , nz.p =  7 
#> alpha =  0.3009873 , tau =  6.289355 , beta = 1.199812 , sigma_e =  0.01502956 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2993358 , tau =  6.340658 , beta = 1.200029 , sigma_e =  0.01501489 , lik =  73840.68 nz =  8 , nz.p =  7 
#> alpha =  0.2987362 , tau =  6.327046 , beta = 1.199318 , sigma_e =  0.01503888 , lik =  73840.67 nz =  8 , nz.p =  7 
#> alpha =  0.2964786 , tau =  6.362476 , beta = 1.20174 , sigma_e =  0.01500474 , lik =  73840.58 nz =  8 , nz.p =  7 
#> alpha =  0.2998977 , tau =  6.339549 , beta = 1.200156 , sigma_e =  0.0150446 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2961066 , tau =  6.425941 , beta = 1.201682 , sigma_e =  0.01503906 , lik =  73840.61 nz =  8 , nz.p =  7 
#> alpha =  0.2997596 , tau =  6.323227 , beta = 1.200286 , sigma_e =  0.01503193 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2987464 , tau =  6.386118 , beta = 1.202585 , sigma_e =  0.01502696 , lik =  73840.69 nz =  8 , nz.p =  7 
#> alpha =  0.2987438 , tau =  6.371299 , beta = 1.201767 , sigma_e =  0.01502994 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2979119 , tau =  6.384708 , beta = 1.202564 , sigma_e =  0.015057 , lik =  73840.68 nz =  8 , nz.p =  7 
#> alpha =  0.2989792 , tau =  6.351642 , beta = 1.200664 , sigma_e =  0.01502541 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2982539 , tau =  6.432659 , beta = 1.196883 , sigma_e =  0.01502939 , lik =  73840.65 nz =  8 , nz.p =  7 
#> alpha =  0.2988049 , tau =  6.340429 , beta = 1.202427 , sigma_e =  0.01503554 , lik =  73840.72 nz =  8 , nz.p =  7 
#> alpha =  0.2971075 , tau =  6.426164 , beta = 1.201713 , sigma_e =  0.0150365 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2977683 , tau =  6.400273 , beta = 1.201358 , sigma_e =  0.01503536 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2988837 , tau =  6.365627 , beta = 1.200736 , sigma_e =  0.0150386 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2967335 , tau =  6.437896 , beta = 1.201987 , sigma_e =  0.01504507 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2972933 , tau =  6.416223 , beta = 1.201658 , sigma_e =  0.01504015 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2963897 , tau =  6.436097 , beta = 1.200842 , sigma_e =  0.01504462 , lik =  73840.71 nz =  8 , nz.p =  7 
#> alpha =  0.2981535 , tau =  6.387437 , beta = 1.201537 , sigma_e =  0.01503361 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2960706 , tau =  6.486569 , beta = 1.199829 , sigma_e =  0.01503823 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2967518 , tau =  6.449722 , beta = 1.200478 , sigma_e =  0.01503756 , lik =  73840.73 nz =  8 , nz.p =  7 
#> alpha =  0.2951461 , tau =  6.4949 , beta = 1.201405 , sigma_e =  0.01503476 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2932949 , tau =  6.560518 , beta = 1.201725 , sigma_e =  0.01503285 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2939309 , tau =  6.542318 , beta = 1.200685 , sigma_e =  0.01503945 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2949809 , tau =  6.503249 , beta = 1.200901 , sigma_e =  0.01503799 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2934432 , tau =  6.563866 , beta = 1.200179 , sigma_e =  0.01503263 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2915368 , tau =  6.638956 , beta = 1.199431 , sigma_e =  0.01502888 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.291412 , tau =  6.639672 , beta = 1.199207 , sigma_e =  0.01503323 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2886053 , tau =  6.749071 , beta = 1.197938 , sigma_e =  0.0150316 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2910323 , tau =  6.662658 , beta = 1.200612 , sigma_e =  0.0150332 , lik =  73840.75 nz =  8 , nz.p =  7 
#> alpha =  0.2887046 , tau =  6.760887 , beta = 1.200901 , sigma_e =  0.015032 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2868721 , tau =  6.811252 , beta = 1.20025 , sigma_e =  0.01502901 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2893108 , tau =  6.719013 , beta = 1.200329 , sigma_e =  0.01503114 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2868229 , tau =  6.823928 , beta = 1.199904 , sigma_e =  0.01502294 , lik =  73840.74 nz =  8 , nz.p =  7 
#> alpha =  0.2885836 , tau =  6.752409 , beta = 1.20011 , sigma_e =  0.01502706 , lik =  73840.75 nz =  8 , nz.p =  7 
#> alpha =  0.285589 , tau =  6.883925 , beta = 1.198224 , sigma_e =  0.01502722 , lik =  73840.75 nz =  8 , nz.p =  7 
#> alpha =  0.2874962 , tau =  6.801609 , beta = 1.199107 , sigma_e =  0.01502863 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2857116 , tau =  6.868774 , beta = 1.200382 , sigma_e =  0.0150311 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2871569 , tau =  6.810584 , beta = 1.200153 , sigma_e =  0.01503054 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2841572 , tau =  6.938212 , beta = 1.200964 , sigma_e =  0.01502566 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2859538 , tau =  6.862341 , beta = 1.200539 , sigma_e =  0.01502755 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2858934 , tau =  6.866586 , beta = 1.200267 , sigma_e =  0.01503203 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2845577 , tau =  6.924397 , beta = 1.200339 , sigma_e =  0.01503451 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2857978 , tau =  6.865877 , beta = 1.201768 , sigma_e =  0.01503282 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2862215 , tau =  6.849754 , beta = 1.201103 , sigma_e =  0.01503177 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2861595 , tau =  6.871632 , beta = 1.200968 , sigma_e =  0.01503354 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2863375 , tau =  6.856487 , beta = 1.200789 , sigma_e =  0.01503241 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2854781 , tau =  6.896881 , beta = 1.201389 , sigma_e =  0.01503321 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2858969 , tau =  6.875205 , beta = 1.201081 , sigma_e =  0.01503254 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2860873 , tau =  6.861412 , beta = 1.20074 , sigma_e =  0.01503028 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2827269 , tau =  7.002765 , beta = 1.200893 , sigma_e =  0.01503332 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2842096 , tau =  6.941496 , beta = 1.200903 , sigma_e =  0.01503299 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2837819 , tau =  6.973264 , beta = 1.200629 , sigma_e =  0.01503418 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2825699 , tau =  7.035851 , beta = 1.200387 , sigma_e =  0.01503538 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2825135 , tau =  7.031676 , beta = 1.200846 , sigma_e =  0.01503771 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2834027 , tau =  6.988718 , beta = 1.200822 , sigma_e =  0.01503585 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2813506 , tau =  7.068828 , beta = 1.200553 , sigma_e =  0.01503537 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2789766 , tau =  7.169538 , beta = 1.200325 , sigma_e =  0.01503628 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2794353 , tau =  7.153029 , beta = 1.199709 , sigma_e =  0.01503693 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2809339 , tau =  7.088114 , beta = 1.200135 , sigma_e =  0.015036 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2814374 , tau =  7.079093 , beta = 1.19992 , sigma_e =  0.01503788 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2807948 , tau =  7.117569 , beta = 1.199433 , sigma_e =  0.01504017 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2781416 , tau =  7.238443 , beta = 1.200082 , sigma_e =  0.01503896 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2829399 , tau =  7.001607 , beta = 1.200284 , sigma_e =  0.01503562 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2790927 , tau =  7.177111 , beta = 1.1994 , sigma_e =  0.01503753 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.280164 , tau =  7.129542 , beta = 1.199759 , sigma_e =  0.01503711 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2785323 , tau =  7.185942 , beta = 1.199441 , sigma_e =  0.01503886 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2765351 , tau =  7.262184 , beta = 1.198956 , sigma_e =  0.01504059 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2783914 , tau =  7.202518 , beta = 1.199231 , sigma_e =  0.01504008 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2771287 , tau =  7.26041 , beta = 1.198775 , sigma_e =  0.01504211 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741325 , tau =  7.398142 , beta = 1.198436 , sigma_e =  0.01504305 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2763083 , tau =  7.296951 , beta = 1.198913 , sigma_e =  0.01504119 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2769603 , tau =  7.275943 , beta = 1.197873 , sigma_e =  0.01504436 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.277463 , tau =  7.249195 , beta = 1.198486 , sigma_e =  0.01504234 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2761974 , tau =  7.297378 , beta = 1.198427 , sigma_e =  0.01504503 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2769184 , tau =  7.267124 , beta = 1.198672 , sigma_e =  0.01504316 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2727164 , tau =  7.432226 , beta = 1.197955 , sigma_e =  0.01504434 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2747139 , tau =  7.352281 , beta = 1.198338 , sigma_e =  0.0150433 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2765037 , tau =  7.271457 , beta = 1.198282 , sigma_e =  0.01504416 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2764548 , tau =  7.277822 , beta = 1.19844 , sigma_e =  0.01504342 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2762448 , tau =  7.312557 , beta = 1.198032 , sigma_e =  0.01504589 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2763174 , tau =  7.299931 , beta = 1.198263 , sigma_e =  0.01504456 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2748657 , tau =  7.346126 , beta = 1.198408 , sigma_e =  0.01504503 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2735762 , tau =  7.395076 , beta = 1.198363 , sigma_e =  0.01504638 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2771558 , tau =  7.25996 , beta = 1.19857 , sigma_e =  0.0150453 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2765433 , tau =  7.282931 , beta = 1.198513 , sigma_e =  0.0150448 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2758053 , tau =  7.308786 , beta = 1.198518 , sigma_e =  0.01504348 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2756094 , tau =  7.314496 , beta = 1.198563 , sigma_e =  0.0150427 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2752111 , tau =  7.343174 , beta = 1.198554 , sigma_e =  0.0150448 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2755215 , tau =  7.326782 , beta = 1.198526 , sigma_e =  0.01504446 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2750305 , tau =  7.331694 , beta = 1.198835 , sigma_e =  0.01504362 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2753517 , tau =  7.32374 , beta = 1.198692 , sigma_e =  0.01504385 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2743315 , tau =  7.365251 , beta = 1.198654 , sigma_e =  0.015043 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2726433 , tau =  7.430588 , beta = 1.198334 , sigma_e =  0.01504604 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2737578 , tau =  7.387673 , beta = 1.198449 , sigma_e =  0.01504506 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2762572 , tau =  7.292353 , beta = 1.198788 , sigma_e =  0.01504125 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2755844 , tau =  7.317899 , beta = 1.198683 , sigma_e =  0.01504253 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2745686 , tau =  7.361034 , beta = 1.198459 , sigma_e =  0.01504324 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2763366 , tau =  7.294201 , beta = 1.198655 , sigma_e =  0.01504119 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2748923 , tau =  7.341722 , beta = 1.198635 , sigma_e =  0.01504048 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2750495 , tau =  7.337984 , beta = 1.198607 , sigma_e =  0.01504148 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.394526 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.387135 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.387135 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.365007 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.372376 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.372376 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.75 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.387135 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.372376 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2747268 , tau =  7.379752 , beta = 1.197793 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.19904 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.197095 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.387135 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.372376 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2736301 , tau =  7.379752 , beta = 1.198889 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.199589 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.197644 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.01505799 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.19904 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.199589 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.200289 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.197095 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.197644 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.196398 , sigma_e =  0.01504294 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.01505799 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.75 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01505799 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.01505799 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.01505799 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.01505799 , lik =  73840.76 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01507306 , lik =  73840.69 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.387135 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.77 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.372376 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2744523 , tau =  7.379752 , beta = 1.198067 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2739039 , tau =  7.379752 , beta = 1.198616 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.199315 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.19737 , sigma_e =  0.0150279 , lik =  73840.78 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01504294 , lik =  73840.79 nz =  8 , nz.p =  7 
#> alpha =  0.2741779 , tau =  7.379752 , beta = 1.198342 , sigma_e =  0.01501288 , lik =  73840.72 nz =  8 , nz.p =  7
#> Warning in rspde_lme(y ~ -1, loc = "loc", repl = "rep", data = data, model =
#> op, : All optimization methods failed to provide a numerically
#> positive-definite Hessian. The optimization method with largest likelihood was
#> chosen. You can try to obtain a positive-definite Hessian by setting
#> 'improve_hessian' to TRUE.

# Compare estimated and true parameter values
rbind(c(fit$coeff$random_effects[c("alpha", "beta", "tau", "kappa")], fit$coeff$measurement_error), 
      c(alpha, beta, tau, kappa, sigma.e))
#>          alpha     beta      tau    kappa   std. dev
#> [1,] 0.2741779 1.198342 7.379752 16.42387 0.01504294
#> [2,] 0.3000000 1.200000 7.000000 15.00000 0.01500000
```

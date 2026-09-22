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
#> tau 0.124812 0.0271957  0.0792896 0.122180   0.185688 0.116990
#> nu  0.968878 0.0769668  0.8203080 0.967942   1.121990 0.964971
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
#> 1       tau  0.2 0.12481241 0.11699049
#> 2        nu  0.8 0.96887843 0.96497053
#> 3   sigma.e  0.1 0.09786177 0.09800655
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
#> tau 0.176702 0.00914724   0.158071 0.177065   0.193787 0.178719
#> nu  0.925094 0.01327720   0.901391 0.924163   0.953232 0.921062
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
#> 1       tau  0.2 0.17670226 0.17871882
#> 2        nu  0.9 0.92509398 0.92106179
#> 3   sigma.e  0.1 0.09994821 0.09960707
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
#>     Pre = 0.149, Running = 20.8, Post = 0.0565, Total = 21 
#> Random effects:
#>   Name     Model
#>     field CGeneric
#> 
#> Model hyperparameters:
#>                                           mean    sd 0.025quant 0.5quant
#> Precision for the Gaussian observations 100.64 4.494      92.09   100.53
#> Theta1 for field                         -5.98 0.048      -6.07    -5.98
#> Theta2 for field                          2.35 0.087       2.17     2.35
#>                                         0.975quant   mode
#> Precision for the Gaussian observations     109.78 100.31
#> Theta1 for field                             -5.88  -5.98
#> Theta2 for field                              2.52   2.35
#> 
#> Marginal log-Likelihood:  727.87 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```

To get a summary of the fit of the random field only, we can do the
following:

``` r

result_fit <- rspde.result(rspde_fit, "field", rspde_model)
summary(result_fit)
#>              mean          sd 0.025quant    0.5quant  0.975quant       mode
#> tau    0.00253262 0.000121607 0.00230336  0.00252838  0.00278045  0.0025204
#> kappa 10.49390000 0.907260000 8.80740000 10.46260000 12.36790000 10.4089000
tau <- op$tau
result_df <- data.frame(
  parameter = c("tau", "kappa"),
  true = c(tau, kappa), mean = c(result_fit$summary.tau$mean,
                                     result_fit$summary.kappa$mean),
  mode = c(result_fit$summary.tau$mode, result_fit$summary.kappa$mode)
)
print(result_df)
#>   parameter    true         mean         mode
#> 1       tau  0.0025  0.002532622  0.002520404
#> 2     kappa 10.0000 10.493912989 10.408923875
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
                 model = op, mean_correction = TRUE, parallel = FALSE,
                 model_options = list(fix_alpha = 0))
#> alpha =  0 , tau =  10 , beta = 1.1 , sigma_e =  0.0192686 , lik =  73414.15 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10 , beta = 1.1 , sigma_e =  0.0192686 , lik =  73414.15 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  14.84277 , beta = 1.1 , sigma_e =  0.0192686 , lik =  72966.64 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10 , beta = 1.390566 , sigma_e =  0.0192686 , lik =  73038.77 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10 , beta = 1.1 , sigma_e =  0.02859994 , lik =  66278.7 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  13.01198 , beta = 1.280719 , sigma_e =  0.0129818 , lik =  76747.43 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  14.84277 , beta = 1.390566 , sigma_e =  0.008746214 , lik =  71913.45 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  8.029969 , beta = 1.430515 , sigma_e =  0.01480835 , lik =  76684.05 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.362982 , beta = 1.333836 , sigma_e =  0.01581585 , lik =  76036.16 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.29686 , beta = 1.145525 , sigma_e =  0.01242447 , lik =  78750.16 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.44858 , beta = 1.049587 , sigma_e =  0.009976811 , lik =  79587.57 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.60253 , beta = 1.403688 , sigma_e =  0.008011349 , lik =  75427.89 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.44858 , beta = 1.31574 , sigma_e =  0.009976811 , lik =  78464.31 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  15.7371 , beta = 1.033743 , sigma_e =  0.008011349 , lik =  77179.49 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  13.30062 , beta = 1.11331 , sigma_e =  0.00934127 , lik =  78780.34 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.854795 , beta = 1.041607 , sigma_e =  0.007338228 , lik =  77163.66 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.56383 , beta = 1.093454 , sigma_e =  0.008463055 , lik =  78980.75 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  12.3626 , beta = 0.9192927 , sigma_e =  0.008556408 , lik =  78907.79 nz =  0 , nz.p =  0 
#> alpha =  0 , tau =  11.8535 , beta = 0.9951946 , sigma_e =  0.008891326 , lik =  79053.11 nz =  0 , nz.p =  0 
#> alpha =  0 , tau =  8.993838 , beta = 0.9835632 , sigma_e =  0.00884269 , lik =  78594.49 nz =  0 , nz.p =  0 
#> alpha =  0 , tau =  12.06119 , beta = 1.077927 , sigma_e =  0.009214049 , lik =  79305.07 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  12.37014 , beta = 0.9909874 , sigma_e =  0.01032949 , lik =  79328.19 nz =  0 , nz.p =  0 
#> alpha =  0 , tau =  11.89149 , beta = 1.014813 , sigma_e =  0.00982745 , lik =  79559.02 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  11.0478 , beta = 1.103852 , sigma_e =  0.01051047 , lik =  79520.04 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  11.24394 , beta = 1.074635 , sigma_e =  0.01007995 , lik =  79585.22 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.36139 , beta = 1.015441 , sigma_e =  0.01076821 , lik =  79370.01 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.76244 , beta = 1.030399 , sigma_e =  0.01035667 , lik =  79543.04 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  11.61174 , beta = 1.061629 , sigma_e =  0.009580182 , lik =  79525.31 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.96876 , beta = 1.03804 , sigma_e =  0.01015684 , lik =  79586.64 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.958314 , beta = 1.095908 , sigma_e =  0.01032044 , lik =  79589.22 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.112983 , beta = 1.141127 , sigma_e =  0.01057613 , lik =  79523.11 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.83992 , beta = 1.067589 , sigma_e =  0.01011511 , lik =  79607.81 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.878444 , beta = 1.105368 , sigma_e =  0.01011615 , lik =  79611.59 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.374629 , beta = 1.142128 , sigma_e =  0.01009586 , lik =  79605.36 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.989483 , beta = 1.132099 , sigma_e =  0.01039435 , lik =  79581.14 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.33186 , beta = 1.069146 , sigma_e =  0.0100796 , lik =  79604.01 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.74178 , beta = 1.065378 , sigma_e =  0.009891323 , lik =  79610.67 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.54032 , beta = 1.072861 , sigma_e =  0.0099969 , lik =  79614.77 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.49221 , beta = 1.094536 , sigma_e =  0.0100722 , lik =  79615.15 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.57331 , beta = 1.107652 , sigma_e =  0.0100685 , lik =  79602.33 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.785434 , beta = 1.114889 , sigma_e =  0.01000843 , lik =  79616.05 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.297306 , beta = 1.139997 , sigma_e =  0.00995551 , lik =  79603.9 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.67035 , beta = 1.082545 , sigma_e =  0.009936233 , lik =  79614.11 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.46661 , beta = 1.088169 , sigma_e =  0.009980909 , lik =  79618.27 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.953632 , beta = 1.12652 , sigma_e =  0.01004403 , lik =  79615 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.09717 , beta = 1.112652 , sigma_e =  0.01003223 , lik =  79618.57 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.746686 , beta = 1.115881 , sigma_e =  0.009942557 , lik =  79615.05 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.30064 , beta = 1.099802 , sigma_e =  0.01003963 , lik =  79618.72 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.81434 , beta = 1.085714 , sigma_e =  0.01002669 , lik =  79607.13 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.0331 , beta = 1.107462 , sigma_e =  0.01001299 , lik =  79619.43 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.829392 , beta = 1.12564 , sigma_e =  0.01007587 , lik =  79615.63 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.30353 , beta = 1.097321 , sigma_e =  0.01000456 , lik =  79619.48 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.32739 , beta = 1.090576 , sigma_e =  0.01000589 , lik =  79619.1 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.26935 , beta = 1.096019 , sigma_e =  0.01001247 , lik =  79619.59 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.10288 , beta = 1.10069 , sigma_e =  0.009980472 , lik =  79619.46 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15196 , beta = 1.100468 , sigma_e =  0.009995228 , lik =  79619.83 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.45405 , beta = 1.088554 , sigma_e =  0.009995186 , lik =  79618.4 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.13672 , beta = 1.102678 , sigma_e =  0.01000854 , lik =  79619.83 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.06949 , beta = 1.102119 , sigma_e =  0.01000625 , lik =  79619.16 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.24452 , beta = 1.098517 , sigma_e =  0.01000499 , lik =  79619.79 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.08671 , beta = 1.105119 , sigma_e =  0.009993373 , lik =  79619.83 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  9.996611 , beta = 1.109721 , sigma_e =  0.009983839 , lik =  79619.43 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.08671 , beta = 1.105119 , sigma_e =  0.009993373 , lik =  79619.83 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.0968 , beta = 1.105119 , sigma_e =  0.009993373 , lik =  79619.88 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.07663 , beta = 1.105119 , sigma_e =  0.009993373 , lik =  79619.77 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.08671 , beta = 1.105725 , sigma_e =  0.009993373 , lik =  79619.86 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.08671 , beta = 1.104514 , sigma_e =  0.009993373 , lik =  79619.8 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.08671 , beta = 1.105119 , sigma_e =  0.01000337 , lik =  79619.77 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.08671 , beta = 1.105119 , sigma_e =  0.009983384 , lik =  79619.85 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.07923 , beta = 1.410647 , sigma_e =  0.00583406 , lik =  45987.34 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.10032 , beta = 1.410647 , sigma_e =  0.00583406 , lik =  45938.07 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.05816 , beta = 1.410647 , sigma_e =  0.00583406 , lik =  46036.55 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.07923 , beta = 1.411558 , sigma_e =  0.00583406 , lik =  45922.26 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.07923 , beta = 1.409736 , sigma_e =  0.00583406 , lik =  46052.3 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.07923 , beta = 1.410647 , sigma_e =  0.005839897 , lik =  46046.82 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  21.07923 , beta = 1.410647 , sigma_e =  0.005828229 , lik =  45927.75 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11328 , beta = 1.106003 , sigma_e =  0.009974193 , lik =  79619.9 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.1234 , beta = 1.106003 , sigma_e =  0.009974193 , lik =  79619.89 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.10317 , beta = 1.106003 , sigma_e =  0.009974193 , lik =  79619.9 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11328 , beta = 1.106609 , sigma_e =  0.009974193 , lik =  79619.86 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11328 , beta = 1.105397 , sigma_e =  0.009974193 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11328 , beta = 1.106003 , sigma_e =  0.009984172 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11328 , beta = 1.106003 , sigma_e =  0.009964224 , lik =  79619.83 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11014 , beta = 1.105606 , sigma_e =  0.009984055 , lik =  79619.94 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.12025 , beta = 1.105606 , sigma_e =  0.009984055 , lik =  79619.94 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.10003 , beta = 1.105606 , sigma_e =  0.009984055 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11014 , beta = 1.106212 , sigma_e =  0.009984055 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11014 , beta = 1.105001 , sigma_e =  0.009984055 , lik =  79619.94 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11014 , beta = 1.105606 , sigma_e =  0.009994044 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11014 , beta = 1.105606 , sigma_e =  0.009974076 , lik =  79619.91 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11352 , beta = 1.105353 , sigma_e =  0.009985671 , lik =  79619.94 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.12364 , beta = 1.105353 , sigma_e =  0.009985671 , lik =  79619.95 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.10341 , beta = 1.105353 , sigma_e =  0.009985671 , lik =  79619.93 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11352 , beta = 1.105959 , sigma_e =  0.009985671 , lik =  79619.93 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11352 , beta = 1.104748 , sigma_e =  0.009985671 , lik =  79619.95 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11352 , beta = 1.105353 , sigma_e =  0.009995661 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.11352 , beta = 1.105353 , sigma_e =  0.00997569 , lik =  79619.93 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.13932 , beta = 1.103662 , sigma_e =  0.009989344 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.14946 , beta = 1.103662 , sigma_e =  0.009989344 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.12918 , beta = 1.103662 , sigma_e =  0.009989344 , lik =  79619.96 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.13932 , beta = 1.104266 , sigma_e =  0.009989344 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.13932 , beta = 1.103059 , sigma_e =  0.009989344 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.13932 , beta = 1.103662 , sigma_e =  0.009999338 , lik =  79619.94 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.13932 , beta = 1.103662 , sigma_e =  0.00997936 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15817 , beta = 1.102485 , sigma_e =  0.009986621 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16833 , beta = 1.102485 , sigma_e =  0.009986621 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.14801 , beta = 1.102485 , sigma_e =  0.009986621 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15817 , beta = 1.101883 , sigma_e =  0.009986621 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15817 , beta = 1.102485 , sigma_e =  0.009996612 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15817 , beta = 1.102485 , sigma_e =  0.009976639 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.23391 , beta = 1.0978 , sigma_e =  0.009975736 , lik =  79619.88 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.24415 , beta = 1.0978 , sigma_e =  0.009975736 , lik =  79619.86 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.22368 , beta = 1.0978 , sigma_e =  0.009975736 , lik =  79619.87 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.23391 , beta = 1.098398 , sigma_e =  0.009975736 , lik =  79619.89 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.23391 , beta = 1.097203 , sigma_e =  0.009975736 , lik =  79619.84 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.23391 , beta = 1.0978 , sigma_e =  0.009985716 , lik =  79619.89 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.23391 , beta = 1.0978 , sigma_e =  0.009965765 , lik =  79619.82 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.18335 , beta = 1.10092 , sigma_e =  0.009982991 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.19354 , beta = 1.10092 , sigma_e =  0.009982991 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17317 , beta = 1.10092 , sigma_e =  0.009982991 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.18335 , beta = 1.101521 , sigma_e =  0.009982991 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.18335 , beta = 1.100319 , sigma_e =  0.009982991 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.18335 , beta = 1.10092 , sigma_e =  0.009992979 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.18335 , beta = 1.10092 , sigma_e =  0.009973013 , lik =  79619.96 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17112 , beta = 1.102287 , sigma_e =  0.009986153 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.1813 , beta = 1.102287 , sigma_e =  0.009986153 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16095 , beta = 1.102287 , sigma_e =  0.009986153 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17112 , beta = 1.102889 , sigma_e =  0.009986153 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17112 , beta = 1.101685 , sigma_e =  0.009986153 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17112 , beta = 1.102287 , sigma_e =  0.009996144 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17112 , beta = 1.102287 , sigma_e =  0.009976171 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16149 , beta = 1.102428 , sigma_e =  0.009984134 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17166 , beta = 1.102428 , sigma_e =  0.009984134 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15133 , beta = 1.102428 , sigma_e =  0.009984134 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16149 , beta = 1.103031 , sigma_e =  0.009984134 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16149 , beta = 1.101826 , sigma_e =  0.009984134 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16149 , beta = 1.102428 , sigma_e =  0.009994123 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16149 , beta = 1.102428 , sigma_e =  0.009974155 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16499 , beta = 1.102299 , sigma_e =  0.009984282 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17516 , beta = 1.102299 , sigma_e =  0.009984282 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15483 , beta = 1.102299 , sigma_e =  0.009984282 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16499 , beta = 1.102901 , sigma_e =  0.009984282 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16499 , beta = 1.101697 , sigma_e =  0.009984282 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16499 , beta = 1.102299 , sigma_e =  0.009994272 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16499 , beta = 1.102299 , sigma_e =  0.009974303 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102905 , sigma_e =  0.009984241 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.1017 , sigma_e =  0.009984241 , lik =  79619.99 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009994231 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009974262 , lik =  79619.98 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.18532 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102905 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.1017 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102302 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102302 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.14466 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102905 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.1017 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102302 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102302 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102905 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102905 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.103508 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102905 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102905 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.1017 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.1017 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.101099 , sigma_e =  0.009984241 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.1017 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.1017 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102302 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102302 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102905 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.1017 , sigma_e =  0.009994231 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.01000423 , lik =  79619.92 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.17514 , beta = 1.102302 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.15481 , beta = 1.102302 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102905 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.1017 , sigma_e =  0.009974262 , lik =  79619.97 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009984241 , lik =  79620 nz =  3 , nz.p =  2 
#> alpha =  0 , tau =  10.16497 , beta = 1.102302 , sigma_e =  0.009964293 , lik =  79619.92 nz =  3 , nz.p =  2

rbind(c(fit$coeff$random_effects[c("beta", "tau")], fit$coeff$measurement_error), 
      c(beta, tau, sigma.e))
#>          beta      tau    std. dev
#> [1,] 1.102302 10.16497 0.009984241
#> [2,] 1.100000 10.00000 0.010000000
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
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01612686 , lik =  -22448.37 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01612686 , lik =  -22448.37 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  10.57653 , beta = 1.05 , sigma_e =  0.01612686 , lik =  -91247.25 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01612686 , lik =  -118916.9 nz =  8 , nz.p =  7 
#> alpha =  1.964212 , tau =  7 , beta = 1.05 , sigma_e =  0.01612686 , lik =  -678604.3 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.586479 , sigma_e =  0.01612686 , lik =  9579.203 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.0243666 , lik =  1131.112 nz =  8 , nz.p =  7 
#> alpha =  0.8603959 , tau =  8.256501 , beta = 1.238475 , sigma_e =  0.01902164 , lik =  51637.92 nz =  8 , nz.p =  7 
#> alpha =  0.569447 , tau =  8.966956 , beta = 1.345043 , sigma_e =  0.02065841 , lik =  66091.95 nz =  8 , nz.p =  7 
#> alpha =  0.9344312 , tau =  9.116221 , beta = 1.367433 , sigma_e =  0.02100229 , lik =  56545.01 nz =  8 , nz.p =  7 
#> alpha =  1.014837 , tau =  8.533666 , beta = 1.28005 , sigma_e =  0.01966018 , lik =  46590.45 nz =  8 , nz.p =  7 
#> alpha =  0.818821 , tau =  5.685368 , beta = 1.519822 , sigma_e =  0.02334283 , lik =  64838.8 nz =  8 , nz.p =  7 
#> alpha =  0.9191313 , tau =  6.639795 , beta = 1.385613 , sigma_e =  0.02128152 , lik =  58478.03 nz =  8 , nz.p =  7 
#> alpha =  0.6805925 , tau =  7.904379 , beta = 1.762115 , sigma_e =  0.02706418 , lik =  64348.54 nz =  8 , nz.p =  7 
#> alpha =  0.8001128 , tau =  7.667881 , beta = 1.548186 , sigma_e =  0.02377846 , lik =  61912.53 nz =  8 , nz.p =  7 
#> alpha =  0.5253693 , tau =  8.298043 , beta = 2.167565 , sigma_e =  0.01868056 , lik =  65310.5 nz =  8 , nz.p =  7 
#> alpha =  0.3656565 , tau =  8.882326 , beta = 1.759677 , sigma_e =  0.02993442 , lik =  66148.82 nz =  8 , nz.p =  7 
#> alpha =  0.193927 , tau =  10.00555 , beta = 1.951189 , sigma_e =  0.0407832 , lik =  60283.37 nz =  8 , nz.p =  7 
#> alpha =  0.3495261 , tau =  6.749858 , beta = 2.171096 , sigma_e =  0.02648596 , lik =  66852.7 nz =  8 , nz.p =  7 
#> alpha =  0.2137695 , tau =  5.808106 , beta = 2.742492 , sigma_e =  0.02974336 , lik =  62512.48 nz =  8 , nz.p =  7 
#> alpha =  0.3676011 , tau =  7.308235 , beta = 1.788692 , sigma_e =  0.02037679 , lik =  71263.1 nz =  8 , nz.p =  7 
#> alpha =  0.2701603 , tau =  7.027241 , beta = 1.83563 , sigma_e =  0.01768098 , lik =  72795.83 nz =  8 , nz.p =  7 
#> alpha =  0.1960867 , tau =  11.05942 , beta = 2.266989 , sigma_e =  0.02114065 , lik =  68048.22 nz =  8 , nz.p =  7 
#> alpha =  0.2803066 , tau =  9.364591 , beta = 2.061106 , sigma_e =  0.02167091 , lik =  68566.01 nz =  8 , nz.p =  7 
#> alpha =  0.2377058 , tau =  7.952805 , beta = 1.551818 , sigma_e =  0.02803347 , lik =  67662.72 nz =  8 , nz.p =  7 
#> alpha =  0.2898318 , tau =  8.037744 , beta = 1.678469 , sigma_e =  0.02532825 , lik =  69100.82 nz =  8 , nz.p =  7 
#> alpha =  0.1673895 , tau =  7.044442 , beta = 2.48488 , sigma_e =  0.02750255 , lik =  66011.55 nz =  8 , nz.p =  7 
#> alpha =  0.4192973 , tau =  8.442002 , beta = 1.593485 , sigma_e =  0.02219044 , lik =  69506.58 nz =  8 , nz.p =  7 
#> alpha =  0.2753517 , tau =  6.96839 , beta = 1.955358 , sigma_e =  0.01683826 , lik =  71463.34 nz =  8 , nz.p =  7 
#> alpha =  0.2955862 , tau =  7.404249 , beta = 1.908155 , sigma_e =  0.01944311 , lik =  71424.1 nz =  8 , nz.p =  7 
#> alpha =  0.2618452 , tau =  9.287517 , beta = 1.546665 , sigma_e =  0.01587688 , lik =  72793.37 nz =  8 , nz.p =  7 
#> alpha =  0.2814514 , tau =  8.575281 , beta = 1.677209 , sigma_e =  0.01804379 , lik =  72397.48 nz =  8 , nz.p =  7 
#> alpha =  0.317728 , tau =  6.672212 , beta = 1.434725 , sigma_e =  0.017134 , lik =  73122.44 nz =  8 , nz.p =  7 
#> alpha =  0.3382724 , tau =  5.63197 , beta = 1.194648 , sigma_e =  0.01523525 , lik =  73734.19 nz =  8 , nz.p =  7 
#> alpha =  0.326882 , tau =  6.742716 , beta = 1.537696 , sigma_e =  0.01196268 , lik =  71126.57 nz =  8 , nz.p =  7 
#> alpha =  0.3171975 , tau =  7.045466 , beta = 1.572603 , sigma_e =  0.01443021 , lik =  72839.53 nz =  8 , nz.p =  7 
#> alpha =  0.2020869 , tau =  5.971856 , beta = 1.5782 , sigma_e =  0.01149504 , lik =  71860.97 nz =  8 , nz.p =  7 
#> alpha =  0.2425404 , tau =  6.511685 , beta = 1.592285 , sigma_e =  0.01354953 , lik =  73334.88 nz =  8 , nz.p =  7 
#> alpha =  0.2925034 , tau =  7.0412 , beta = 1.216145 , sigma_e =  0.01388769 , lik =  73577.53 nz =  8 , nz.p =  7 
#> alpha =  0.2881178 , tau =  7.022927 , beta = 1.366294 , sigma_e =  0.01457292 , lik =  73534.05 nz =  8 , nz.p =  7 
#> alpha =  0.3215204 , tau =  4.729098 , beta = 1.380666 , sigma_e =  0.01396073 , lik =  73654.24 nz =  8 , nz.p =  7 
#> alpha =  0.3054341 , tau =  5.598332 , beta = 1.422495 , sigma_e =  0.01441691 , lik =  73669.18 nz =  8 , nz.p =  7 
#> alpha =  0.327303 , tau =  5.706932 , beta = 1.052639 , sigma_e =  0.01155349 , lik =  72018.68 nz =  8 , nz.p =  7 
#> alpha =  0.2834351 , tau =  6.670975 , beta = 1.598158 , sigma_e =  0.01589674 , lik =  73471.51 nz =  8 , nz.p =  7 
#> alpha =  0.2664964 , tau =  5.569159 , beta = 1.247791 , sigma_e =  0.01471502 , lik =  73724.19 nz =  8 , nz.p =  7 
#> alpha =  0.2783562 , tau =  5.906356 , beta = 1.319029 , sigma_e =  0.0146433 , lik =  73778.97 nz =  8 , nz.p =  7 
#> alpha =  0.3682679 , tau =  5.795041 , beta = 1.108211 , sigma_e =  0.01616572 , lik =  73651.35 nz =  8 , nz.p =  7 
#> alpha =  0.3317561 , tau =  5.966446 , beta = 1.224097 , sigma_e =  0.01546776 , lik =  73747.2 nz =  8 , nz.p =  7 
#> alpha =  0.3356142 , tau =  5.409674 , beta = 1.011607 , sigma_e =  0.01362885 , lik =  73494.27 nz =  8 , nz.p =  7 
#> alpha =  0.3217315 , tau =  5.700667 , beta = 1.135667 , sigma_e =  0.01416353 , lik =  73729.08 nz =  8 , nz.p =  7 
#> alpha =  0.3378377 , tau =  4.710048 , beta = 1.298685 , sigma_e =  0.01572367 , lik =  73753.01 nz =  8 , nz.p =  7 
#> alpha =  0.3258846 , tau =  5.208116 , beta = 1.277225 , sigma_e =  0.01524309 , lik =  73788.6 nz =  8 , nz.p =  7 
#> alpha =  0.332032 , tau =  5.755223 , beta = 1.06235 , sigma_e =  0.01548787 , lik =  73734.58 nz =  8 , nz.p =  7 
#> alpha =  0.3192115 , tau =  5.673483 , beta = 1.297075 , sigma_e =  0.01633882 , lik =  73623.48 nz =  8 , nz.p =  7 
#> alpha =  0.3210997 , tau =  5.693859 , beta = 1.173764 , sigma_e =  0.01467857 , lik =  73792.62 nz =  8 , nz.p =  7 
#> alpha =  0.297341 , tau =  5.767845 , beta = 1.221908 , sigma_e =  0.01496499 , lik =  73801.19 nz =  8 , nz.p =  7 
#> alpha =  0.2898487 , tau =  5.649271 , beta = 1.450452 , sigma_e =  0.01452004 , lik =  73715.54 nz =  8 , nz.p =  7 
#> alpha =  0.3209429 , tau =  5.72855 , beta = 1.149621 , sigma_e =  0.01524003 , lik =  73797.33 nz =  8 , nz.p =  7 
#> alpha =  0.3197482 , tau =  5.80904 , beta = 1.22595 , sigma_e =  0.01520756 , lik =  73791.97 nz =  8 , nz.p =  7 
#> alpha =  0.3606424 , tau =  5.379976 , beta = 1.098047 , sigma_e =  0.01549931 , lik =  73755.42 nz =  8 , nz.p =  7 
#> alpha =  0.296975 , tau =  5.77012 , beta = 1.26396 , sigma_e =  0.01485276 , lik =  73801.25 nz =  8 , nz.p =  7 
#> alpha =  0.2968085 , tau =  6.356544 , beta = 1.14215 , sigma_e =  0.01473575 , lik =  73795.29 nz =  8 , nz.p =  7 
#> alpha =  0.3038248 , tau =  6.04764 , beta = 1.173726 , sigma_e =  0.01486098 , lik =  73806.26 nz =  8 , nz.p =  7 
#> alpha =  0.2963863 , tau =  5.791485 , beta = 1.16778 , sigma_e =  0.0146346 , lik =  73783.53 nz =  8 , nz.p =  7 
#> alpha =  0.3137405 , tau =  5.804646 , beta = 1.21095 , sigma_e =  0.01506225 , lik =  73805.65 nz =  8 , nz.p =  7 
#> alpha =  0.2924144 , tau =  5.954362 , beta = 1.233064 , sigma_e =  0.0153193 , lik =  73796.39 nz =  8 , nz.p =  7 
#> alpha =  0.299336 , tau =  5.88814 , beta = 1.218402 , sigma_e =  0.01515654 , lik =  73806.56 nz =  8 , nz.p =  7 
#> alpha =  0.2845137 , tau =  5.983711 , beta = 1.286785 , sigma_e =  0.01472253 , lik =  73792.11 nz =  8 , nz.p =  7 
#> alpha =  0.3114201 , tau =  5.791302 , beta = 1.183402 , sigma_e =  0.01510897 , lik =  73807.53 nz =  8 , nz.p =  7 
#> alpha =  0.3128333 , tau =  5.9526 , beta = 1.197444 , sigma_e =  0.01505066 , lik =  73794.09 nz =  8 , nz.p =  7 
#> alpha =  0.3011406 , tau =  5.813489 , beta = 1.215858 , sigma_e =  0.01498636 , lik =  73808.11 nz =  8 , nz.p =  7 
#> alpha =  0.3149689 , tau =  5.968103 , beta = 1.139331 , sigma_e =  0.01521879 , lik =  73794.48 nz =  8 , nz.p =  7 
#> alpha =  0.3013747 , tau =  5.818991 , beta = 1.231833 , sigma_e =  0.01494344 , lik =  73807.96 nz =  8 , nz.p =  7 
#> alpha =  0.2933806 , tau =  5.938464 , beta = 1.197951 , sigma_e =  0.01495966 , lik =  73803.74 nz =  8 , nz.p =  7 
#> alpha =  0.3085218 , tau =  5.837816 , beta = 1.207725 , sigma_e =  0.01503654 , lik =  73808.94 nz =  8 , nz.p =  7 
#> alpha =  0.3048203 , tau =  5.619915 , beta = 1.250577 , sigma_e =  0.01523367 , lik =  73802.28 nz =  8 , nz.p =  7 
#> alpha =  0.3040734 , tau =  5.93775 , beta = 1.192384 , sigma_e =  0.01495328 , lik =  73809.21 nz =  8 , nz.p =  7 
#> alpha =  0.3113407 , tau =  5.791554 , beta = 1.193851 , sigma_e =  0.01485615 , lik =  73806.39 nz =  8 , nz.p =  7 
#> alpha =  0.3022931 , tau =  5.863843 , beta = 1.212293 , sigma_e =  0.01508088 , lik =  73808.91 nz =  8 , nz.p =  7 
#> alpha =  0.2957201 , tau =  5.917789 , beta = 1.240749 , sigma_e =  0.01489183 , lik =  73806.58 nz =  8 , nz.p =  7 
#> alpha =  0.3074186 , tau =  5.822668 , beta = 1.197656 , sigma_e =  0.01505439 , lik =  73809.28 nz =  8 , nz.p =  7 
#> alpha =  0.3080136 , tau =  5.891118 , beta = 1.178989 , sigma_e =  0.01510142 , lik =  73806.51 nz =  8 , nz.p =  7 
#> alpha =  0.3030209 , tau =  5.83694 , beta = 1.218437 , sigma_e =  0.01498278 , lik =  73809.37 nz =  8 , nz.p =  7 
#> alpha =  0.3090215 , tau =  5.906196 , beta = 1.195469 , sigma_e =  0.01505672 , lik =  73807.09 nz =  8 , nz.p =  7 
#> alpha =  0.3030918 , tau =  5.836529 , beta = 1.210765 , sigma_e =  0.01500392 , lik =  73809.39 nz =  8 , nz.p =  7 
#> alpha =  0.3081685 , tau =  5.844553 , beta = 1.198421 , sigma_e =  0.01493177 , lik =  73808.79 nz =  8 , nz.p =  7 
#> alpha =  0.3037514 , tau =  5.859015 , beta = 1.20883 , sigma_e =  0.01504346 , lik =  73809.51 nz =  8 , nz.p =  7 
#> alpha =  0.3000707 , tau =  5.879129 , beta = 1.203422 , sigma_e =  0.01497856 , lik =  73808.62 nz =  8 , nz.p =  7 
#> alpha =  0.306387 , tau =  5.848117 , beta = 1.206656 , sigma_e =  0.01502202 , lik =  73809.55 nz =  8 , nz.p =  7 
#> alpha =  0.3053851 , tau =  5.74512 , beta = 1.224818 , sigma_e =  0.01508961 , lik =  73808.51 nz =  8 , nz.p =  7 
#> alpha =  0.3044008 , tau =  5.888995 , beta = 1.200388 , sigma_e =  0.01498725 , lik =  73809.64 nz =  8 , nz.p =  7 
#> alpha =  0.3008723 , tau =  5.885273 , beta = 1.220381 , sigma_e =  0.01496149 , lik =  73809.22 nz =  8 , nz.p =  7 
#> alpha =  0.3057688 , tau =  5.838256 , beta = 1.203326 , sigma_e =  0.01503111 , lik =  73809.6 nz =  8 , nz.p =  7 
#> alpha =  0.3063431 , tau =  5.871413 , beta = 1.193648 , sigma_e =  0.01505238 , lik =  73809.05 nz =  8 , nz.p =  7 
#> alpha =  0.3038481 , tau =  5.845539 , beta = 1.2122 , sigma_e =  0.01500015 , lik =  73809.64 nz =  8 , nz.p =  7 
#> alpha =  0.3065769 , tau =  5.875451 , beta = 1.201781 , sigma_e =  0.01502966 , lik =  73809.19 nz =  8 , nz.p =  7 
#> alpha =  0.3039593 , tau =  5.846235 , beta = 1.20852 , sigma_e =  0.01501035 , lik =  73809.63 nz =  8 , nz.p =  7 
#> alpha =  0.305995 , tau =  5.847792 , beta = 1.203593 , sigma_e =  0.01497695 , lik =  73809.59 nz =  8 , nz.p =  7 
#> alpha =  0.3054325 , tau =  5.850595 , beta = 1.204903 , sigma_e =  0.01499355 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3029844 , tau =  5.859683 , beta = 1.205064 , sigma_e =  0.01498695 , lik =  73809.52 nz =  8 , nz.p =  7 
#> alpha =  0.3055327 , tau =  5.851006 , beta = 1.206259 , sigma_e =  0.01501324 , lik =  73809.67 nz =  8 , nz.p =  7 
#> alpha =  0.3035031 , tau =  5.874703 , beta = 1.20957 , sigma_e =  0.01497076 , lik =  73809.6 nz =  8 , nz.p =  7 
#> alpha =  0.3040679 , tau =  5.86557 , beta = 1.208009 , sigma_e =  0.01498582 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3053535 , tau =  5.874438 , beta = 1.204173 , sigma_e =  0.01498166 , lik =  73809.52 nz =  8 , nz.p =  7 
#> alpha =  0.3043073 , tau =  5.853273 , beta = 1.207433 , sigma_e =  0.01500317 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048732 , tau =  5.817608 , beta = 1.21519 , sigma_e =  0.01501113 , lik =  73809.55 nz =  8 , nz.p =  7 
#> alpha =  0.3045188 , tau =  5.871067 , beta = 1.204067 , sigma_e =  0.01499322 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3043093 , tau =  5.851914 , beta = 1.209164 , sigma_e =  0.01499897 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3035242 , tau =  5.86596 , beta = 1.207166 , sigma_e =  0.01497667 , lik =  73809.63 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.860598 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.848889 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3053345 , tau =  5.854741 , beta = 1.206182 , sigma_e =  0.01500409 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3047245 , tau =  5.854741 , beta = 1.206792 , sigma_e =  0.01500409 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.207499 , sigma_e =  0.01500409 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.205476 , sigma_e =  0.01500409 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.0150191 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.0149891 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9405383 , sigma_e =  0.01346469 , lik =  66906.95 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.050654 , beta = 0.9405383 , sigma_e =  0.01346469 , lik =  66915.45 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.044559 , beta = 0.9405383 , sigma_e =  0.01346469 , lik =  66898.46 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9405383 , sigma_e =  0.01346469 , lik =  66908.37 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9405383 , sigma_e =  0.01346469 , lik =  66905.54 nz =  0 , nz.p =  0 
#> alpha =  0.185152 , tau =  3.047605 , beta = 0.9403532 , sigma_e =  0.01346469 , lik =  66911.16 nz =  0 , nz.p =  0 
#> alpha =  0.1847821 , tau =  3.047605 , beta = 0.9407231 , sigma_e =  0.01346469 , lik =  66902.76 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9411641 , sigma_e =  0.01346469 , lik =  66908.46 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9399131 , sigma_e =  0.01346469 , lik =  66905.47 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9405383 , sigma_e =  0.01347816 , lik =  66911.45 nz =  0 , nz.p =  0 
#> alpha =  0.184967 , tau =  3.047605 , beta = 0.9405383 , sigma_e =  0.01345123 , lik =  66902.4 nz =  0 , nz.p =  0 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.206212 , sigma_e =  0.01500275 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.857444 , beta = 1.206212 , sigma_e =  0.01500275 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.845741 , beta = 1.206212 , sigma_e =  0.01500275 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.206212 , sigma_e =  0.01500275 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.206212 , sigma_e =  0.01500275 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3052086 , tau =  5.85159 , beta = 1.205907 , sigma_e =  0.01500275 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3045988 , tau =  5.85159 , beta = 1.206517 , sigma_e =  0.01500275 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.207224 , sigma_e =  0.01500275 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.205202 , sigma_e =  0.01500275 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.206212 , sigma_e =  0.01501776 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3049036 , tau =  5.85159 , beta = 1.206212 , sigma_e =  0.01498776 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.851346 , beta = 1.206217 , sigma_e =  0.01500127 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.8572 , beta = 1.206217 , sigma_e =  0.01500127 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.845498 , beta = 1.206217 , sigma_e =  0.01500127 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.851346 , beta = 1.206217 , sigma_e =  0.01500127 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051967 , tau =  5.851346 , beta = 1.205912 , sigma_e =  0.01500127 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045869 , tau =  5.851346 , beta = 1.206522 , sigma_e =  0.01500127 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.851346 , beta = 1.207228 , sigma_e =  0.01500127 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.851346 , beta = 1.205206 , sigma_e =  0.01500127 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.851346 , beta = 1.206217 , sigma_e =  0.01501628 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048917 , tau =  5.851346 , beta = 1.206217 , sigma_e =  0.01498628 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.206235 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.856225 , beta = 1.206235 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.844524 , beta = 1.206235 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.206235 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.206235 , sigma_e =  0.01499535 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051491 , tau =  5.850372 , beta = 1.20593 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3045394 , tau =  5.850372 , beta = 1.20654 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.207247 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.205225 , sigma_e =  0.01499535 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.206235 , sigma_e =  0.01501035 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048441 , tau =  5.850372 , beta = 1.206235 , sigma_e =  0.01498036 , lik =  73809.66 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.206223 , sigma_e =  0.0149993 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.856875 , beta = 1.206223 , sigma_e =  0.0149993 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.206223 , sigma_e =  0.0149993 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.206223 , sigma_e =  0.0149993 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051809 , tau =  5.851021 , beta = 1.205918 , sigma_e =  0.0149993 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045711 , tau =  5.851021 , beta = 1.206528 , sigma_e =  0.0149993 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.207234 , sigma_e =  0.0149993 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.205212 , sigma_e =  0.0149993 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.206223 , sigma_e =  0.01501431 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048758 , tau =  5.851021 , beta = 1.206223 , sigma_e =  0.01498431 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.862928 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.857068 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.857068 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.67 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.839523 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.845366 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.845366 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.857068 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.845366 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3054956 , tau =  5.851214 , beta = 1.205609 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.206926 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.204904 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.857068 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.845366 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3042761 , tau =  5.851214 , beta = 1.206828 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.207536 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.205513 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.206926 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.207536 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.208243 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.204904 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.205513 , sigma_e =  0.01500047 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.204199 , sigma_e =  0.01500047 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01501548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01501548 , lik =  73809.69 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01501548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.0150305 , lik =  73809.63 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.857068 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.67 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.845366 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3051903 , tau =  5.851214 , beta = 1.205914 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3045805 , tau =  5.851214 , beta = 1.206524 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.207231 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.205209 , sigma_e =  0.01498548 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.01500047 , lik =  73809.71 nz =  8 , nz.p =  7 
#> alpha =  0.3048852 , tau =  5.851214 , beta = 1.206219 , sigma_e =  0.0149705 , lik =  73809.61 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73809.7 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  8.910184 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  72845.05 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73727.77 nz =  8 , nz.p =  7 
#> alpha =  0.4642166 , tau =  5.854741 , beta = 1.0473 , sigma_e =  0.01500409 , lik =  72811.42 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.734373 , sigma_e =  0.01500409 , lik =  73387.12 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.02283436 , lik =  70714.11 nz =  8 , nz.p =  7 
#> alpha =  0.3608222 , tau =  6.92563 , beta = 1.335711 , sigma_e =  0.009858951 , lik =  66007.86 nz =  8 , nz.p =  7 
#> alpha =  0.3181116 , tau =  6.105842 , beta = 1.236787 , sigma_e =  0.01850969 , lik =  72951 nz =  8 , nz.p =  7 
#> alpha =  0.2038252 , tau =  7.042947 , beta = 1.512976 , sigma_e =  0.0163187 , lik =  73597.11 nz =  8 , nz.p =  7 
#> alpha =  0.2503938 , tau =  6.725007 , beta = 1.411478 , sigma_e =  0.01597963 , lik =  73657.08 nz =  8 , nz.p =  7 
#> alpha =  0.2866484 , tau =  4.135213 , beta = 1.499508 , sigma_e =  0.0167351 , lik =  73275.33 nz =  8 , nz.p =  7 
#> alpha =  0.2911371 , tau =  5.010084 , beta = 1.420057 , sigma_e =  0.01628447 , lik =  73619.06 nz =  8 , nz.p =  7 
#> alpha =  0.2652888 , tau =  5.575435 , beta = 1.538049 , sigma_e =  0.01288809 , lik =  73104.87 nz =  8 , nz.p =  7 
#> alpha =  0.2776093 , tau =  5.703553 , beta = 1.45861 , sigma_e =  0.01410884 , lik =  73618.89 nz =  8 , nz.p =  7 
#> alpha =  0.2664365 , tau =  5.754025 , beta = 1.051676 , sigma_e =  0.01551282 , lik =  73576.05 nz =  8 , nz.p =  7 
#> alpha =  0.2756009 , tau =  5.779041 , beta = 1.182581 , sigma_e =  0.01538404 , lik =  73778.4 nz =  8 , nz.p =  7 
#> alpha =  0.291904 , tau =  5.937698 , beta = 1.128145 , sigma_e =  0.01707815 , lik =  73433.44 nz =  8 , nz.p =  7 
#> alpha =  0.2980023 , tau =  5.415971 , beta = 1.308859 , sigma_e =  0.01563118 , lik =  73761.33 nz =  8 , nz.p =  7 
#> alpha =  0.3050294 , tau =  5.854741 , beta = 1.206487 , sigma_e =  0.01500409 , lik =  73789.16 nz =  8 , nz.p =  7 
#> alpha =  0.2763647 , tau =  6.274805 , beta = 1.307726 , sigma_e =  0.01548418 , lik =  73766.1 nz =  8 , nz.p =  7 
#> alpha =  0.289942 , tau =  5.816767 , beta = 1.194546 , sigma_e =  0.01519288 , lik =  73801.98 nz =  8 , nz.p =  7 
#> alpha =  0.2909965 , tau =  5.778652 , beta = 1.327241 , sigma_e =  0.01454958 , lik =  73754.78 nz =  8 , nz.p =  7 
#> alpha =  0.2983945 , tau =  5.896073 , beta = 1.166305 , sigma_e =  0.01600757 , lik =  73712.3 nz =  8 , nz.p =  7 
#> alpha =  0.2928287 , tau =  5.807786 , beta = 1.284873 , sigma_e =  0.01490113 , lik =  73794.32 nz =  8 , nz.p =  7 
#> alpha =  0.2893454 , tau =  6.469144 , beta = 1.174916 , sigma_e =  0.01461757 , lik =  73769.76 nz =  8 , nz.p =  7 
#> alpha =  0.2914858 , tau =  6.188055 , beta = 1.206604 , sigma_e =  0.01486463 , lik =  73796.68 nz =  8 , nz.p =  7 
#> alpha =  0.3187187 , tau =  5.552704 , beta = 1.13378 , sigma_e =  0.01451727 , lik =  73777.64 nz =  8 , nz.p =  7 
#> alpha =  0.3075575 , tau =  5.725039 , beta = 1.17626 , sigma_e =  0.01475319 , lik =  73801.7 nz =  8 , nz.p =  7 
#> alpha =  0.2897226 , tau =  5.898009 , beta = 1.219841 , sigma_e =  0.01488106 , lik =  73802.34 nz =  8 , nz.p =  7 
#> alpha =  0.2934757 , tau =  5.887162 , beta = 1.216575 , sigma_e =  0.01491172 , lik =  73808.11 nz =  8 , nz.p =  7 
#> alpha =  0.3020605 , tau =  5.978064 , beta = 1.121246 , sigma_e =  0.01498814 , lik =  73790.5 nz =  8 , nz.p =  7 
#> alpha =  0.2951098 , tau =  5.849896 , beta = 1.241727 , sigma_e =  0.01492284 , lik =  73805.05 nz =  8 , nz.p =  7 
#> alpha =  0.3049556 , tau =  5.485985 , beta = 1.207302 , sigma_e =  0.01504845 , lik =  73797.71 nz =  8 , nz.p =  7 
#> alpha =  0.3015309 , tau =  5.653657 , beta = 1.207166 , sigma_e =  0.01500228 , lik =  73806.51 nz =  8 , nz.p =  7 
#> alpha =  0.286741 , tau =  5.899982 , beta = 1.250512 , sigma_e =  0.01526401 , lik =  73797.84 nz =  8 , nz.p =  7 
#> alpha =  0.3022158 , tau =  5.768283 , beta = 1.194697 , sigma_e =  0.01487927 , lik =  73807.92 nz =  8 , nz.p =  7 
#> alpha =  0.3092484 , tau =  5.787532 , beta = 1.232511 , sigma_e =  0.01469911 , lik =  73796.89 nz =  8 , nz.p =  7 
#> alpha =  0.2946526 , tau =  5.809445 , beta = 1.203851 , sigma_e =  0.0150679 , lik =  73807.99 nz =  8 , nz.p =  7 
#> alpha =  0.3036446 , tau =  5.738807 , beta = 1.170724 , sigma_e =  0.01502313 , lik =  73799.27 nz =  8 , nz.p =  7 
#> alpha =  0.2972208 , tau =  5.821923 , beta = 1.223633 , sigma_e =  0.01494785 , lik =  73808.67 nz =  8 , nz.p =  7 
#> alpha =  0.295472 , tau =  6.008071 , beta = 1.210889 , sigma_e =  0.01492186 , lik =  73807.29 nz =  8 , nz.p =  7 
#> alpha =  0.2969752 , tau =  5.917437 , beta = 1.20997 , sigma_e =  0.01494192 , lik =  73808.95 nz =  8 , nz.p =  7 
#> alpha =  0.2927468 , tau =  5.949121 , beta = 1.229597 , sigma_e =  0.01507053 , lik =  73808.36 nz =  8 , nz.p =  7 
#> alpha =  0.2950859 , tau =  5.903387 , beta = 1.22084 , sigma_e =  0.01502249 , lik =  73809.22 nz =  8 , nz.p =  7 
#> alpha =  0.3004377 , tau =  5.944993 , beta = 1.227364 , sigma_e =  0.01486391 , lik =  73805.99 nz =  8 , nz.p =  7 
#> alpha =  0.2960883 , tau =  5.843039 , beta = 1.209661 , sigma_e =  0.01501664 , lik =  73809.26 nz =  8 , nz.p =  7 
#> alpha =  0.3027137 , tau =  5.848886 , beta = 1.211604 , sigma_e =  0.01506177 , lik =  73808.54 nz =  8 , nz.p =  7 
#> alpha =  0.3003773 , tau =  5.858432 , beta = 1.212872 , sigma_e =  0.01502412 , lik =  73809.31 nz =  8 , nz.p =  7 
#> alpha =  0.3001652 , tau =  5.929234 , beta = 1.200422 , sigma_e =  0.01505599 , lik =  73809.24 nz =  8 , nz.p =  7 
#> alpha =  0.2994264 , tau =  5.902223 , beta = 1.20619 , sigma_e =  0.01502888 , lik =  73809.56 nz =  8 , nz.p =  7 
#> alpha =  0.3014027 , tau =  5.827526 , beta = 1.212471 , sigma_e =  0.01509696 , lik =  73808.35 nz =  8 , nz.p =  7 
#> alpha =  0.298076 , tau =  5.89483 , beta = 1.210596 , sigma_e =  0.01498053 , lik =  73809.48 nz =  8 , nz.p =  7 
#> alpha =  0.3045583 , tau =  5.838007 , beta = 1.197476 , sigma_e =  0.01499921 , lik =  73809.39 nz =  8 , nz.p =  7 
#> alpha =  0.302162 , tau =  5.854284 , beta = 1.203327 , sigma_e =  0.01500502 , lik =  73809.58 nz =  8 , nz.p =  7 
#> alpha =  0.3060025 , tau =  5.902841 , beta = 1.206066 , sigma_e =  0.0150004 , lik =  73809.02 nz =  8 , nz.p =  7 
#> alpha =  0.2985363 , tau =  5.857933 , beta = 1.208789 , sigma_e =  0.01501258 , lik =  73809.55 nz =  8 , nz.p =  7 
#> alpha =  0.3008924 , tau =  5.887131 , beta = 1.201331 , sigma_e =  0.01498833 , lik =  73809.48 nz =  8 , nz.p =  7 
#> alpha =  0.3007635 , tau =  5.879943 , beta = 1.204205 , sigma_e =  0.01499727 , lik =  73809.61 nz =  8 , nz.p =  7 
#> alpha =  0.3043061 , tau =  5.844866 , beta = 1.200985 , sigma_e =  0.01503866 , lik =  73809.41 nz =  8 , nz.p =  7 
#> alpha =  0.2996214 , tau =  5.882299 , beta = 1.208205 , sigma_e =  0.01499504 , lik =  73809.61 nz =  8 , nz.p =  7 
#> alpha =  0.3042781 , tau =  5.891455 , beta = 1.202558 , sigma_e =  0.01499954 , lik =  73809.64 nz =  8 , nz.p =  7 
#> alpha =  0.3028324 , tau =  5.883056 , beta = 1.204126 , sigma_e =  0.0150028 , lik =  73809.66 nz =  8 , nz.p =  7 
#> alpha =  0.3047493 , tau =  5.839643 , beta = 1.204336 , sigma_e =  0.01497286 , lik =  73809.53 nz =  8 , nz.p =  7 
#> alpha =  0.3007483 , tau =  5.886515 , beta = 1.205734 , sigma_e =  0.01501486 , lik =  73809.63 nz =  8 , nz.p =  7 
#> alpha =  0.3014242 , tau =  5.900406 , beta = 1.208187 , sigma_e =  0.0150006 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3016085 , tau =  5.888841 , beta = 1.20697 , sigma_e =  0.0150017 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3047849 , tau =  5.874916 , beta = 1.202786 , sigma_e =  0.01501325 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3034858 , tau =  5.876761 , beta = 1.204149 , sigma_e =  0.0150087 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3047239 , tau =  5.875997 , beta = 1.206783 , sigma_e =  0.01501559 , lik =  73809.66 nz =  8 , nz.p =  7 
#> alpha =  0.3037289 , tau =  5.876984 , beta = 1.206139 , sigma_e =  0.01501101 , lik =  73809.68 nz =  8 , nz.p =  7 
#> alpha =  0.3059438 , tau =  5.865634 , beta = 1.2054 , sigma_e =  0.01499647 , lik =  73809.72 nz =  8 , nz.p =  7 
#> alpha =  0.308575 , tau =  5.855221 , beta = 1.205208 , sigma_e =  0.01498728 , lik =  73809.74 nz =  8 , nz.p =  7 
#> alpha =  0.30613 , tau =  5.857959 , beta = 1.207471 , sigma_e =  0.01500231 , lik =  73809.75 nz =  8 , nz.p =  7 
#> alpha =  0.3077923 , tau =  5.84545 , beta = 1.209146 , sigma_e =  0.01500207 , lik =  73809.77 nz =  8 , nz.p =  7 
#> alpha =  0.3098776 , tau =  5.834918 , beta = 1.205441 , sigma_e =  0.01500355 , lik =  73809.77 nz =  8 , nz.p =  7 
#> alpha =  0.3140968 , tau =  5.808141 , beta = 1.204609 , sigma_e =  0.01500447 , lik =  73809.8 nz =  8 , nz.p =  7 
#> alpha =  0.3122235 , tau =  5.819507 , beta = 1.208498 , sigma_e =  0.01499487 , lik =  73809.75 nz =  8 , nz.p =  7 
#> alpha =  0.3100158 , tau =  5.833768 , beta = 1.207418 , sigma_e =  0.01499833 , lik =  73809.76 nz =  8 , nz.p =  7 
#> alpha =  0.3145407 , tau =  5.802132 , beta = 1.206968 , sigma_e =  0.01498749 , lik =  73809.84 nz =  8 , nz.p =  7 
#> alpha =  0.3200901 , tau =  5.765065 , beta = 1.20729 , sigma_e =  0.01497575 , lik =  73809.87 nz =  8 , nz.p =  7 
#> alpha =  0.3192961 , tau =  5.788327 , beta = 1.206916 , sigma_e =  0.01498307 , lik =  73809.91 nz =  8 , nz.p =  7 
#> alpha =  0.3266778 , tau =  5.755403 , beta = 1.206962 , sigma_e =  0.01497256 , lik =  73809.99 nz =  8 , nz.p =  7 
#> alpha =  0.3229074 , tau =  5.748181 , beta = 1.208986 , sigma_e =  0.01499398 , lik =  73809.99 nz =  8 , nz.p =  7 
#> alpha =  0.3192631 , tau =  5.774756 , beta = 1.208073 , sigma_e =  0.01499231 , lik =  73809.93 nz =  8 , nz.p =  7 
#> alpha =  0.3266876 , tau =  5.735312 , beta = 1.207323 , sigma_e =  0.0149812 , lik =  73810.04 nz =  8 , nz.p =  7 
#> alpha =  0.3353568 , tau =  5.686709 , beta = 1.207044 , sigma_e =  0.01497264 , lik =  73810.13 nz =  8 , nz.p =  7 
#> alpha =  0.3405324 , tau =  5.66116 , beta = 1.204267 , sigma_e =  0.01496571 , lik =  73810.13 nz =  8 , nz.p =  7 
#> alpha =  0.3581862 , tau =  5.571205 , beta = 1.200828 , sigma_e =  0.01494756 , lik =  73810.21 nz =  8 , nz.p =  7 
#> alpha =  0.3516987 , tau =  5.603392 , beta = 1.207604 , sigma_e =  0.01494058 , lik =  73810.05 nz =  8 , nz.p =  7 
#> alpha =  0.3418959 , tau =  5.653893 , beta = 1.207108 , sigma_e =  0.01495653 , lik =  73810.11 nz =  8 , nz.p =  7 
#> alpha =  0.3543323 , tau =  5.601461 , beta = 1.204815 , sigma_e =  0.01496155 , lik =  73810.39 nz =  8 , nz.p =  7 
#> alpha =  0.3728035 , tau =  5.521409 , beta = 1.202593 , sigma_e =  0.01495445 , lik =  73810.53 nz =  8 , nz.p =  7 
#> alpha =  0.3720196 , tau =  5.528189 , beta = 1.200057 , sigma_e =  0.01492759 , lik =  73810.2 nz =  8 , nz.p =  7 
#> alpha =  0.3590822 , tau =  5.582385 , beta = 1.202804 , sigma_e =  0.01494416 , lik =  73810.28 nz =  8 , nz.p =  7 
#> alpha =  0.3819037 , tau =  5.454252 , beta = 1.199966 , sigma_e =  0.01493759 , lik =  73810.37 nz =  8 , nz.p =  7 
#> alpha =  0.3672783 , tau =  5.528029 , beta = 1.202327 , sigma_e =  0.01494633 , lik =  73810.39 nz =  8 , nz.p =  7 
#> alpha =  0.350006 , tau =  5.615634 , beta = 1.205278 , sigma_e =  0.01495478 , lik =  73810.24 nz =  8 , nz.p =  7 
#> alpha =  0.3894344 , tau =  5.443199 , beta = 1.197046 , sigma_e =  0.0149263 , lik =  73810.34 nz =  8 , nz.p =  7 
#> alpha =  0.3751479 , tau =  5.503081 , beta = 1.200141 , sigma_e =  0.01493787 , lik =  73810.4 nz =  8 , nz.p =  7 
#> alpha =  0.371424 , tau =  5.528773 , beta = 1.204556 , sigma_e =  0.01494747 , lik =  73810.27 nz =  8 , nz.p =  7 
#> alpha =  0.3680694 , tau =  5.539351 , beta = 1.203644 , sigma_e =  0.01494749 , lik =  73810.37 nz =  8 , nz.p =  7 
#> alpha =  0.3878333 , tau =  5.455105 , beta = 1.198643 , sigma_e =  0.01493735 , lik =  73810.4 nz =  8 , nz.p =  7 
#> alpha =  0.3780094 , tau =  5.494802 , beta = 1.200583 , sigma_e =  0.0149417 , lik =  73810.44 nz =  8 , nz.p =  7 
#> alpha =  0.3858782 , tau =  5.452994 , beta = 1.200606 , sigma_e =  0.01494698 , lik =  73810.5 nz =  8 , nz.p =  7 
#> alpha =  0.3789973 , tau =  5.485058 , beta = 1.201284 , sigma_e =  0.01494627 , lik =  73810.49 nz =  8 , nz.p =  7 
#> alpha =  0.3836388 , tau =  5.460928 , beta = 1.198794 , sigma_e =  0.01494344 , lik =  73810.55 nz =  8 , nz.p =  7 
#> alpha =  0.3916688 , tau =  5.422134 , beta = 1.196163 , sigma_e =  0.01494141 , lik =  73810.6 nz =  8 , nz.p =  7 
#> alpha =  0.3944813 , tau =  5.429938 , beta = 1.19744 , sigma_e =  0.01494264 , lik =  73810.55 nz =  8 , nz.p =  7 
#> alpha =  0.3874973 , tau =  5.454297 , beta = 1.198801 , sigma_e =  0.01494356 , lik =  73810.55 nz =  8 , nz.p =  7 
#> alpha =  0.3940477 , tau =  5.42544 , beta = 1.198791 , sigma_e =  0.014953 , lik =  73810.55 nz =  8 , nz.p =  7 
#> alpha =  0.3892353 , tau =  5.444747 , beta = 1.199189 , sigma_e =  0.01494922 , lik =  73810.57 nz =  8 , nz.p =  7 
#> alpha =  0.3956694 , tau =  5.413761 , beta = 1.197783 , sigma_e =  0.01495218 , lik =  73810.68 nz =  8 , nz.p =  7 
#> alpha =  0.4048064 , tau =  5.373689 , beta = 1.196152 , sigma_e =  0.01495742 , lik =  73810.74 nz =  8 , nz.p =  7 
#> alpha =  0.395097 , tau =  5.423393 , beta = 1.196175 , sigma_e =  0.01495108 , lik =  73810.74 nz =  8 , nz.p =  7 
#> alpha =  0.3997887 , tau =  5.408653 , beta = 1.193886 , sigma_e =  0.01495313 , lik =  73810.81 nz =  8 , nz.p =  7 
#> alpha =  0.4205467 , tau =  5.312169 , beta = 1.189434 , sigma_e =  0.01494308 , lik =  73810.52 nz =  8 , nz.p =  7 
#> alpha =  0.3842055 , tau =  5.468339 , beta = 1.199735 , sigma_e =  0.01495161 , lik =  73810.64 nz =  8 , nz.p =  7 
#> alpha =  0.3932626 , tau =  5.416903 , beta = 1.196718 , sigma_e =  0.01495848 , lik =  73810.83 nz =  8 , nz.p =  7 
#> alpha =  0.3926547 , tau =  5.410398 , beta = 1.196357 , sigma_e =  0.0149664 , lik =  73810.92 nz =  8 , nz.p =  7 
#> alpha =  0.3999597 , tau =  5.388513 , beta = 1.193758 , sigma_e =  0.01495877 , lik =  73810.87 nz =  8 , nz.p =  7 
#> alpha =  0.3972512 , tau =  5.402516 , beta = 1.195141 , sigma_e =  0.01495638 , lik =  73810.84 nz =  8 , nz.p =  7 
#> alpha =  0.4008191 , tau =  5.397539 , beta = 1.195858 , sigma_e =  0.01497353 , lik =  73811.01 nz =  8 , nz.p =  7 
#> alpha =  0.4054741 , tau =  5.385284 , beta = 1.195653 , sigma_e =  0.01498962 , lik =  73811.13 nz =  8 , nz.p =  7 
#> alpha =  0.4175069 , tau =  5.319269 , beta = 1.190055 , sigma_e =  0.01497853 , lik =  73811.08 nz =  8 , nz.p =  7 
#> alpha =  0.4089203 , tau =  5.356151 , beta = 1.192689 , sigma_e =  0.01497179 , lik =  73811.03 nz =  8 , nz.p =  7 
#> alpha =  0.4011866 , tau =  5.390965 , beta = 1.191867 , sigma_e =  0.01498115 , lik =  73811.04 nz =  8 , nz.p =  7 
#> alpha =  0.4020885 , tau =  5.386641 , beta = 1.192936 , sigma_e =  0.01497522 , lik =  73811.05 nz =  8 , nz.p =  7 
#> alpha =  0.4071549 , tau =  5.347387 , beta = 1.193722 , sigma_e =  0.0149943 , lik =  73811.21 nz =  8 , nz.p =  7 
#> alpha =  0.4108887 , tau =  5.317015 , beta = 1.193607 , sigma_e =  0.01501494 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4113932 , tau =  5.338772 , beta = 1.193767 , sigma_e =  0.01501114 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4085045 , tau =  5.351164 , beta = 1.193784 , sigma_e =  0.01499803 , lik =  73811.23 nz =  8 , nz.p =  7 
#> alpha =  0.4269344 , tau =  5.288907 , beta = 1.189556 , sigma_e =  0.01502141 , lik =  73811.24 nz =  8 , nz.p =  7 
#> alpha =  0.4180936 , tau =  5.319022 , beta = 1.191463 , sigma_e =  0.01500764 , lik =  73811.26 nz =  8 , nz.p =  7 
#> alpha =  0.4234782 , tau =  5.285457 , beta = 1.192752 , sigma_e =  0.01502556 , lik =  73811.19 nz =  8 , nz.p =  7 
#> alpha =  0.4180264 , tau =  5.310574 , beta = 1.192864 , sigma_e =  0.01501296 , lik =  73811.25 nz =  8 , nz.p =  7 
#> alpha =  0.4080421 , tau =  5.3489 , beta = 1.196878 , sigma_e =  0.01503604 , lik =  73811.24 nz =  8 , nz.p =  7 
#> alpha =  0.410388 , tau =  5.341477 , beta = 1.195192 , sigma_e =  0.01502164 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4221805 , tau =  5.266097 , beta = 1.190993 , sigma_e =  0.01503774 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.43079 , tau =  5.207497 , beta = 1.188456 , sigma_e =  0.01506186 , lik =  73811.15 nz =  8 , nz.p =  7 
#> alpha =  0.4111264 , tau =  5.322247 , beta = 1.193169 , sigma_e =  0.01502428 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4077192 , tau =  5.328093 , beta = 1.193294 , sigma_e =  0.01502994 , lik =  73811.25 nz =  8 , nz.p =  7 
#> alpha =  0.4083065 , tau =  5.315083 , beta = 1.195217 , sigma_e =  0.01503627 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4107316 , tau =  5.316067 , beta = 1.194297 , sigma_e =  0.0150291 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4146895 , tau =  5.2864 , beta = 1.193171 , sigma_e =  0.01503995 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4122148 , tau =  5.325631 , beta = 1.19362 , sigma_e =  0.01501834 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4164445 , tau =  5.277448 , beta = 1.191098 , sigma_e =  0.01502811 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.414922 , tau =  5.293383 , beta = 1.192129 , sigma_e =  0.0150265 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4175654 , tau =  5.29229 , beta = 1.192089 , sigma_e =  0.01503945 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4158861 , tau =  5.29846 , beta = 1.192476 , sigma_e =  0.01503332 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4039625 , tau =  5.356572 , beta = 1.195136 , sigma_e =  0.01501488 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4175504 , tau =  5.288572 , beta = 1.192088 , sigma_e =  0.01503202 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4179664 , tau =  5.295226 , beta = 1.191078 , sigma_e =  0.01502468 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4161459 , tau =  5.300429 , beta = 1.191893 , sigma_e =  0.01502578 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4180466 , tau =  5.275698 , beta = 1.191073 , sigma_e =  0.01503843 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4141757 , tau =  5.313015 , beta = 1.192761 , sigma_e =  0.01502206 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4155335 , tau =  5.296905 , beta = 1.192012 , sigma_e =  0.01502614 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4136285 , tau =  5.311327 , beta = 1.192537 , sigma_e =  0.01502503 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4168475 , tau =  5.294497 , beta = 1.191991 , sigma_e =  0.0150289 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.416016 , tau =  5.299444 , beta = 1.192185 , sigma_e =  0.01502955 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4178678 , tau =  5.290402 , beta = 1.191793 , sigma_e =  0.01502794 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4200037 , tau =  5.27997 , beta = 1.191408 , sigma_e =  0.0150294 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4168852 , tau =  5.3022 , beta = 1.19224 , sigma_e =  0.01502756 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4165469 , tau =  5.300876 , beta = 1.192183 , sigma_e =  0.0150272 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4192081 , tau =  5.281289 , beta = 1.191245 , sigma_e =  0.0150337 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4163463 , tau =  5.300652 , beta = 1.192038 , sigma_e =  0.01502649 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4166972 , tau =  5.297685 , beta = 1.192087 , sigma_e =  0.01502805 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4162813 , tau =  5.30016 , beta = 1.192184 , sigma_e =  0.01502838 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4178753 , tau =  5.291073 , beta = 1.191716 , sigma_e =  0.01503045 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.417245 , tau =  5.295036 , beta = 1.191906 , sigma_e =  0.01502899 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.306238 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.300935 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.300935 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.285056 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.290343 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.290343 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.32 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.300935 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.290343 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.418042 , tau =  5.295636 , beta = 1.191153 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.19268 , sigma_e =  0.01502757 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.190462 , sigma_e =  0.01502757 , lik =  73811.32 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01501255 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.300935 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.290343 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4163732 , tau =  5.295636 , beta = 1.192822 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.193515 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.191297 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01501255 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.19268 , sigma_e =  0.01502757 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.193515 , sigma_e =  0.01502757 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.194209 , sigma_e =  0.01502757 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01504261 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01501255 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.32 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.190462 , sigma_e =  0.01502757 , lik =  73811.32 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.191297 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.189772 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01501255 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4176242 , tau =  5.295636 , beta = 1.191571 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01504261 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01504261 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01505766 , lik =  73811.23 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.300935 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.290343 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.28 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01501255 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4167898 , tau =  5.295636 , beta = 1.192405 , sigma_e =  0.01501255 , lik =  73811.29 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.193098 , sigma_e =  0.01501255 , lik =  73811.27 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.19088 , sigma_e =  0.01501255 , lik =  73811.3 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01502757 , lik =  73811.31 nz =  8 , nz.p =  7 
#> alpha =  0.4172068 , tau =  5.295636 , beta = 1.191988 , sigma_e =  0.01499755 , lik =  73811.22 nz =  8 , nz.p =  7
#> Warning in rspde_lme(y ~ -1, loc = "loc", repl = "rep", data = data, model =
#> op, : All optimization methods failed to provide a numerically
#> positive-definite Hessian. The optimization method with largest likelihood was
#> chosen. You can try to obtain a positive-definite Hessian by setting
#> 'improve_hessian' to TRUE.

# Compare estimated and true parameter values
rbind(c(fit$coeff$random_effects[c("alpha", "beta", "tau", "kappa")], fit$coeff$measurement_error), 
      c(alpha, beta, tau, kappa, sigma.e))
#>          alpha     beta      tau    kappa   std. dev
#> [1,] 0.4172068 1.191988 5.295636 13.21793 0.01502757
#> [2,] 0.3000000 1.200000 7.000000 15.00000 0.01500000
```

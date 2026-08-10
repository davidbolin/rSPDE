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
#> tau 0.124947 0.0266582  0.0797501 0.122644   0.184025 0.118174
#> nu  0.969556 0.0751728  0.8257020 0.968103   1.120290 0.963343
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
#> 1       tau  0.2 0.12494708 0.11817367
#> 2        nu  0.8 0.96955644 0.96334344
#> 3   sigma.e  0.1 0.09779271 0.09812017
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
#>         mean        sd 0.025quant 0.5quant 0.975quant     mode
#> tau 0.177359 0.0084898   0.161014 0.177273   0.194339 0.177229
#> nu  0.924774 0.0114848   0.902582 0.924623   0.947682 0.924125
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
#>   parameter true       mean      mode
#> 1       tau  0.2 0.17735871 0.1772287
#> 2        nu  0.9 0.92477411 0.9241250
#> 3   sigma.e  0.1 0.09994077 0.1000328
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
#>     Pre = 0.147, Running = 22.5, Post = 0.0573, Total = 22.7 
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
#> Precision for the Gaussian observations     109.79 100.31
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
#>             mean          sd 0.025quant    0.5quant  0.975quant        mode
#> tau    0.0025326 0.000121496 0.00230355  0.00252837  0.00278021  0.00252039
#> kappa 10.4939000 0.906296000 8.80894000 10.46260000 12.36570000 10.40910000
tau <- op$tau
result_df <- data.frame(
  parameter = c("tau", "kappa"),
  true = c(tau, kappa), mean = c(result_fit$summary.tau$mean,
                                     result_fit$summary.kappa$mean),
  mode = c(result_fit$summary.tau$mode, result_fit$summary.kappa$mode)
)
print(result_df)
#>   parameter    true         mean         mode
#> 1       tau  0.0025  0.002532603  0.002520395
#> 2     kappa 10.0000 10.493875165 10.409057938
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
#>          beta      tau    std. dev
#> [1,] 1.089129 10.02531 0.009988161
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
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01600249 , lik =  -19908.86 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01600249 , lik =  -19908.86 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  10.58472 , beta = 1.05 , sigma_e =  0.01600249 , lik =  -89098.26 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.05 , sigma_e =  0.01600249 , lik =  -117276.2 nz =  8 , nz.p =  7 
#> alpha =  1.965733 , tau =  7 , beta = 1.05 , sigma_e =  0.01600249 , lik =  -692264.5 nz =  8 , nz.p =  7 
#> alpha =  1.3 , tau =  7 , beta = 1.587708 , sigma_e =  0.01600249 , lik =  13011.82 nz =  8 , nz.p =  7 
#> alpha =  1.533825 , tau =  8.259058 , beta = 1.238859 , sigma_e =  0.01058294 , lik =  -316292.2 nz =  8 , nz.p =  7 
#> alpha =  1.471695 , tau =  7.92451 , beta = 1.188676 , sigma_e =  0.01301359 , lik =  -173869.7 nz =  8 , nz.p =  7 
#> alpha =  0.9034661 , tau =  8.679212 , beta = 1.301882 , sigma_e =  0.01473233 , lik =  38462.9 nz =  8 , nz.p =  7 
#> alpha =  0.6124992 , tau =  9.664322 , beta = 1.449648 , sigma_e =  0.01413557 , lik =  59173.68 nz =  8 , nz.p =  7 
#> alpha =  0.8498335 , tau =  8.300142 , beta = 1.245021 , sigma_e =  0.01872529 , lik =  53235.62 nz =  8 , nz.p =  7 
#> alpha =  0.9748882 , tau =  8.204597 , beta = 1.230689 , sigma_e =  0.01709703 , lik =  38582.59 nz =  8 , nz.p =  7 
#> alpha =  0.8116418 , tau =  10.059 , beta = 1.50885 , sigma_e =  0.01621565 , lik =  57110.82 nz =  8 , nz.p =  7 
#> alpha =  0.9130807 , tau =  9.187358 , beta = 1.378104 , sigma_e =  0.01616209 , lik =  48328 nz =  8 , nz.p =  7 
#> alpha =  0.6722531 , tau =  6.518143 , beta = 1.744331 , sigma_e =  0.0163017 , lik =  66780.89 nz =  8 , nz.p =  7 
#> alpha =  0.4834233 , tau =  5.115006 , beta = 2.264848 , sigma_e =  0.0164534 , lik =  67175.21 nz =  8 , nz.p =  7 
#> alpha =  0.4525708 , tau =  8.693694 , beta = 2.412764 , sigma_e =  0.0164839 , lik =  63010.91 nz =  8 , nz.p =  7 
#> alpha =  0.5891837 , tau =  8.235269 , beta = 1.930709 , sigma_e =  0.01636221 , lik =  64931.31 nz =  8 , nz.p =  7 
#> alpha =  0.3297687 , tau =  9.277613 , beta = 1.861801 , sigma_e =  0.01663113 , lik =  71292.54 nz =  8 , nz.p =  7 
#> alpha =  0.1660896 , tau =  10.68084 , beta = 2.079932 , sigma_e =  0.01695465 , lik =  69409.08 nz =  8 , nz.p =  7 
#> alpha =  0.3454373 , tau =  8.180674 , beta = 2.593097 , sigma_e =  0.01355433 , lik =  61379.93 nz =  8 , nz.p =  7 
#> alpha =  0.4326238 , tau =  8.210379 , beta = 2.128674 , sigma_e =  0.01469487 , lik =  63961.15 nz =  8 , nz.p =  7 
#> alpha =  0.2812141 , tau =  6.22362 , beta = 2.498379 , sigma_e =  0.01504753 , lik =  66502.04 nz =  8 , nz.p =  7 
#> alpha =  0.3665379 , tau =  7.017315 , beta = 2.189611 , sigma_e =  0.01533142 , lik =  66294.69 nz =  8 , nz.p =  7 
#> alpha =  0.2729409 , tau =  5.432193 , beta = 3.070726 , sigma_e =  0.01769892 , lik =  52265.48 nz =  8 , nz.p =  7 
#> alpha =  0.5004334 , tau =  8.368013 , beta = 1.715603 , sigma_e =  0.01495277 , lik =  67909.75 nz =  8 , nz.p =  7 
#> alpha =  0.4095982 , tau =  6.442984 , beta = 1.948681 , sigma_e =  0.01714418 , lik =  70331.18 nz =  8 , nz.p =  7 
#> alpha =  0.4152371 , tau =  6.845514 , beta = 1.991841 , sigma_e =  0.01649601 , lik =  68331.06 nz =  8 , nz.p =  7 
#> alpha =  0.2600523 , tau =  5.822858 , beta = 2.193804 , sigma_e =  0.01568738 , lik =  69414 nz =  8 , nz.p =  7 
#> alpha =  0.3190496 , tau =  6.349968 , beta = 2.128994 , sigma_e =  0.01585343 , lik =  68661.09 nz =  8 , nz.p =  7 
#> alpha =  0.5280606 , tau =  7.502537 , beta = 1.555844 , sigma_e =  0.01734448 , lik =  70207.97 nz =  8 , nz.p =  7 
#> alpha =  0.4510992 , tau =  7.160063 , beta = 1.760645 , sigma_e =  0.0167393 , lik =  71298.76 nz =  8 , nz.p =  7 
#> alpha =  0.2987934 , tau =  10.44319 , beta = 1.615733 , sigma_e =  0.01597215 , lik =  71420.43 nz =  8 , nz.p =  7 
#> alpha =  0.2349053 , tau =  14.92199 , beta = 1.387093 , sigma_e =  0.01573684 , lik =  70266.08 nz =  8 , nz.p =  7 
#> alpha =  0.2348347 , tau =  6.974957 , beta = 1.981178 , sigma_e =  0.01804488 , lik =  71527.16 nz =  8 , nz.p =  7 
#> alpha =  0.1608682 , tau =  6.367971 , beta = 2.05535 , sigma_e =  0.01982304 , lik =  70953.32 nz =  8 , nz.p =  7 
#> alpha =  0.4338389 , tau =  10.77432 , beta = 1.496687 , sigma_e =  0.0181906 , lik =  68912.26 nz =  8 , nz.p =  7 
#> alpha =  0.2955477 , tau =  6.791248 , beta = 2.011807 , sigma_e =  0.01627888 , lik =  69689.35 nz =  8 , nz.p =  7 
#> alpha =  0.3817346 , tau =  9.237967 , beta = 1.664748 , sigma_e =  0.01752963 , lik =  70955.37 nz =  8 , nz.p =  7 
#> alpha =  0.3580784 , tau =  8.554009 , beta = 1.749861 , sigma_e =  0.01720821 , lik =  71614.18 nz =  8 , nz.p =  7 
#> alpha =  0.2610232 , tau =  10.90616 , beta = 1.658495 , sigma_e =  0.01666972 , lik =  72241.52 nz =  8 , nz.p =  7 
#> alpha =  0.2083721 , tau =  14.18941 , beta = 1.532298 , sigma_e =  0.01643744 , lik =  71953.86 nz =  8 , nz.p =  7 
#> alpha =  0.2952674 , tau =  8.080107 , beta = 1.658526 , sigma_e =  0.01719989 , lik =  72720.49 nz =  8 , nz.p =  7 
#> alpha =  0.1821924 , tau =  10.98845 , beta = 1.663603 , sigma_e =  0.01727531 , lik =  72814.57 nz =  8 , nz.p =  7 
#> alpha =  0.1157869 , tau =  13.61277 , beta = 1.577512 , sigma_e =  0.01754971 , lik =  72800.53 nz =  8 , nz.p =  7 
#> alpha =  0.225681 , tau =  7.688851 , beta = 1.870765 , sigma_e =  0.01868204 , lik =  72546.6 nz =  8 , nz.p =  7 
#> alpha =  0.2420829 , tau =  8.300503 , beta = 1.806799 , sigma_e =  0.01796425 , lik =  72734.01 nz =  8 , nz.p =  7 
#> alpha =  0.2906595 , tau =  12.3418 , beta = 1.472458 , sigma_e =  0.01650653 , lik =  70984.88 nz =  8 , nz.p =  7 
#> alpha =  0.2476953 , tau =  8.044531 , beta = 1.841771 , sigma_e =  0.01764735 , lik =  72700.09 nz =  8 , nz.p =  7 
#> alpha =  0.164492 , tau =  9.819348 , beta = 1.677307 , sigma_e =  0.01748444 , lik =  73189.98 nz =  8 , nz.p =  7 
#> alpha =  0.111488 , tau =  10.52057 , beta = 1.614246 , sigma_e =  0.01762421 , lik =  73303.86 nz =  8 , nz.p =  7 
#> alpha =  0.1605696 , tau =  7.591493 , beta = 1.772667 , sigma_e =  0.01845575 , lik =  72923.52 nz =  8 , nz.p =  7 
#> alpha =  0.1813081 , tau =  8.311193 , beta = 1.748487 , sigma_e =  0.01799206 , lik =  73065.87 nz =  8 , nz.p =  7 
#> alpha =  0.1493565 , tau =  10.4271 , beta = 1.577565 , sigma_e =  0.01756875 , lik =  73202.83 nz =  8 , nz.p =  7 
#> alpha =  0.1694914 , tau =  9.772323 , beta = 1.639466 , sigma_e =  0.01758837 , lik =  73156.38 nz =  8 , nz.p =  7 
#> alpha =  0.09539595 , tau =  11.49597 , beta = 1.658008 , sigma_e =  0.01817945 , lik =  73051.77 nz =  8 , nz.p =  7 
#> alpha =  0.1265323 , tau =  10.526 , beta = 1.674218 , sigma_e =  0.01792945 , lik =  73112.3 nz =  8 , nz.p =  7 
#> alpha =  0.08975866 , tau =  12.30614 , beta = 1.507339 , sigma_e =  0.01739244 , lik =  73333.43 nz =  8 , nz.p =  7 
#> alpha =  0.05465533 , tau =  14.98411 , beta = 1.368679 , sigma_e =  0.01711339 , lik =  73388.33 nz =  8 , nz.p =  7 
#> alpha =  0.07368023 , tau =  10.52368 , beta = 1.512337 , sigma_e =  0.01801806 , lik =  73214.22 nz =  8 , nz.p =  7 
#> alpha =  0.09239443 , tau =  10.638 , beta = 1.553441 , sigma_e =  0.01782943 , lik =  73281.45 nz =  8 , nz.p =  7 
#> alpha =  0.05653997 , tau =  15.35746 , beta = 1.379041 , sigma_e =  0.01723756 , lik =  73307.67 nz =  8 , nz.p =  7 
#> alpha =  0.07566081 , tau =  13.17211 , beta = 1.46457 , sigma_e =  0.01742316 , lik =  73336.07 nz =  8 , nz.p =  7 
#> alpha =  0.06595371 , tau =  13.26404 , beta = 1.374496 , sigma_e =  0.01710068 , lik =  73515.24 nz =  8 , nz.p =  7 
#> alpha =  0.04761658 , tau =  14.88957 , beta = 1.252045 , sigma_e =  0.01670078 , lik =  73645.53 nz =  8 , nz.p =  7 
#> alpha =  0.03536959 , tau =  15.44057 , beta = 1.304765 , sigma_e =  0.01710164 , lik =  73517.98 nz =  8 , nz.p =  7 
#> alpha =  0.0507024 , tau =  13.99711 , beta = 1.372859 , sigma_e =  0.01721724 , lik =  73485.29 nz =  8 , nz.p =  7 
#> alpha =  0.07444488 , tau =  12.05916 , beta = 1.471178 , sigma_e =  0.01750667 , lik =  73399.17 nz =  8 , nz.p =  7 
#> alpha =  0.0274576 , tau =  18.75946 , beta = 1.169572 , sigma_e =  0.01672118 , lik =  73455.62 nz =  8 , nz.p =  7 
#> alpha =  0.0389766 , tau =  16.23398 , beta = 1.263694 , sigma_e =  0.0169425 , lik =  73499.03 nz =  8 , nz.p =  7 
#> alpha =  0.03103079 , tau =  16.28966 , beta = 1.209394 , sigma_e =  0.01672591 , lik =  73634.15 nz =  8 , nz.p =  7 
#> alpha =  0.03877594 , tau =  15.44713 , beta = 1.267332 , sigma_e =  0.01689756 , lik =  73579.5 nz =  8 , nz.p =  7 
#> alpha =  0.03425275 , tau =  14.80638 , beta = 1.229408 , sigma_e =  0.01687338 , lik =  73657.51 nz =  8 , nz.p =  7 
#> alpha =  0.02711608 , tau =  14.71831 , beta = 1.167383 , sigma_e =  0.01675463 , lik =  73689.48 nz =  8 , nz.p =  7 
#> alpha =  0.01678997 , tau =  19.92412 , beta = 1.055746 , sigma_e =  0.01620717 , lik =  73619.21 nz =  8 , nz.p =  7 
#> alpha =  0.02436386 , tau =  17.57372 , beta = 1.14121 , sigma_e =  0.01652271 , lik =  73664.5 nz =  8 , nz.p =  7 
#> alpha =  0.02657404 , tau =  15.2773 , beta = 1.165423 , sigma_e =  0.01657962 , lik =  73737.48 nz =  8 , nz.p =  7 
#> alpha =  0.02194242 , tau =  14.82031 , beta = 1.120578 , sigma_e =  0.01640111 , lik =  73743.98 nz =  8 , nz.p =  7 
#> alpha =  0.0241932 , tau =  15.80105 , beta = 1.0697 , sigma_e =  0.01615283 , lik =  73784.57 nz =  8 , nz.p =  7 
#> alpha =  0.02000899 , tau =  15.98443 , beta = 0.9793212 , sigma_e =  0.01569834 , lik =  73668.45 nz =  0 , nz.p =  0 
#> alpha =  0.02496265 , tau =  14.79608 , beta = 1.092011 , sigma_e =  0.01628696 , lik =  73746.3 nz =  8 , nz.p =  7 
#> alpha =  0.02635823 , tau =  15.15612 , beta = 1.1194 , sigma_e =  0.01639561 , lik =  73762.5 nz =  8 , nz.p =  7 
#> alpha =  0.0128407 , tau =  16.30268 , beta = 1.011613 , sigma_e =  0.01619159 , lik =  73678.84 nz =  8 , nz.p =  7 
#> alpha =  0.01781891 , tau =  15.9373 , beta = 1.064964 , sigma_e =  0.01631741 , lik =  73735.59 nz =  8 , nz.p =  7 
#> alpha =  0.02214393 , tau =  13.28306 , beta = 1.075361 , sigma_e =  0.01628444 , lik =  73557.28 nz =  8 , nz.p =  7 
#> alpha =  0.02378884 , tau =  16.38597 , beta = 1.124073 , sigma_e =  0.01646282 , lik =  73743.44 nz =  8 , nz.p =  7 
#> alpha =  0.01888133 , tau =  16.55577 , beta = 1.038036 , sigma_e =  0.01594656 , lik =  73795.03 nz =  8 , nz.p =  7 
#> alpha =  0.01575561 , tau =  17.55881 , beta = 0.9829564 , sigma_e =  0.01555726 , lik =  73719.57 nz =  0 , nz.p =  0 
#> alpha =  0.02019624 , tau =  15.83298 , beta = 1.079108 , sigma_e =  0.016294 , lik =  73758.44 nz =  8 , nz.p =  7 
#> alpha =  0.0206303 , tau =  14.8931 , beta = 1.047645 , sigma_e =  0.0160145 , lik =  73729.52 nz =  8 , nz.p =  7 
#> alpha =  0.02295653 , tau =  15.99928 , beta = 1.10402 , sigma_e =  0.01634958 , lik =  73768.79 nz =  8 , nz.p =  7 
#> alpha =  0.0227742 , tau =  16.97842 , beta = 1.044639 , sigma_e =  0.01605453 , lik =  73780.97 nz =  8 , nz.p =  7 
#> alpha =  0.02256334 , tau =  16.41108 , beta = 1.062761 , sigma_e =  0.01614048 , lik =  73784.25 nz =  8 , nz.p =  7 
#> alpha =  0.02586487 , tau =  16.1221 , beta = 1.07675 , sigma_e =  0.01609899 , lik =  73803.65 nz =  8 , nz.p =  7 
#> alpha =  0.0292705 , tau =  16.26864 , beta = 1.075006 , sigma_e =  0.01600237 , lik =  73815.83 nz =  8 , nz.p =  7 
#> alpha =  0.02066665 , tau =  17.32613 , beta = 1.02394 , sigma_e =  0.01584461 , lik =  73795.14 nz =  8 , nz.p =  7 
#> alpha =  0.02196249 , tau =  16.75611 , beta = 1.04634 , sigma_e =  0.0159806 , lik =  73799.11 nz =  8 , nz.p =  7 
#> alpha =  0.02331257 , tau =  16.71931 , beta = 1.01599 , sigma_e =  0.0157448 , lik =  73804.5 nz =  8 , nz.p =  7 
#> alpha =  0.02322305 , tau =  16.53632 , beta = 1.036776 , sigma_e =  0.01589387 , lik =  73805.08 nz =  8 , nz.p =  7 
#> alpha =  0.02399381 , tau =  16.34941 , beta = 1.043556 , sigma_e =  0.01585085 , lik =  73819.42 nz =  8 , nz.p =  7 
#> alpha =  0.0247427 , tau =  16.31866 , beta = 1.034123 , sigma_e =  0.01570798 , lik =  73831.42 nz =  8 , nz.p =  7 
#> alpha =  0.02257932 , tau =  17.20097 , beta = 1.023345 , sigma_e =  0.01566279 , lik =  73821.41 nz =  8 , nz.p =  7 
#> alpha =  0.02297241 , tau =  16.83977 , beta = 1.034568 , sigma_e =  0.01578389 , lik =  73819.58 nz =  8 , nz.p =  7 
#> alpha =  0.03107729 , tau =  16.66982 , beta = 1.046469 , sigma_e =  0.01575185 , lik =  73825.13 nz =  8 , nz.p =  7 
#> alpha =  0.02743724 , tau =  16.64123 , beta = 1.044881 , sigma_e =  0.0158003 , lik =  73824.37 nz =  8 , nz.p =  7 
#> alpha =  0.03069719 , tau =  16.43653 , beta = 1.038942 , sigma_e =  0.01562791 , lik =  73843.32 nz =  8 , nz.p =  7 
#> alpha =  0.03629168 , tau =  16.27904 , beta = 1.034018 , sigma_e =  0.0154545 , lik =  73853.37 nz =  8 , nz.p =  7 
#> alpha =  0.03470693 , tau =  16.55086 , beta = 1.047613 , sigma_e =  0.01553798 , lik =  73847.31 nz =  8 , nz.p =  7 
#> alpha =  0.03139004 , tau =  16.54722 , beta = 1.045268 , sigma_e =  0.0156262 , lik =  73844.24 nz =  8 , nz.p =  7 
#> alpha =  0.02947811 , tau =  16.93932 , beta = 1.002243 , sigma_e =  0.0152519 , lik =  73835.19 nz =  8 , nz.p =  7 
#> alpha =  0.02942607 , tau =  16.7691 , beta = 1.019573 , sigma_e =  0.01543615 , lik =  73842.17 nz =  8 , nz.p =  7 
#> alpha =  0.0424842 , tau =  15.85903 , beta = 1.047395 , sigma_e =  0.01549198 , lik =  73857.32 nz =  8 , nz.p =  7 
#> alpha =  0.05827546 , tau =  15.22785 , beta = 1.054891 , sigma_e =  0.01540728 , lik =  73847.16 nz =  8 , nz.p =  7 
#> alpha =  0.03497912 , tau =  16.04121 , beta = 1.027335 , sigma_e =  0.01530223 , lik =  73867.64 nz =  8 , nz.p =  7 
#> alpha =  0.03711007 , tau =  15.73585 , beta = 1.017739 , sigma_e =  0.01508226 , lik =  73872.21 nz =  8 , nz.p =  7 
#> alpha =  0.05166813 , tau =  16.14975 , beta = 1.027767 , sigma_e =  0.01509747 , lik =  73824.52 nz =  8 , nz.p =  7 
#> alpha =  0.02974343 , tau =  16.27627 , beta = 1.034195 , sigma_e =  0.01555308 , lik =  73851.29 nz =  8 , nz.p =  7 
#> alpha =  0.0436316 , tau =  15.52954 , beta = 1.052587 , sigma_e =  0.01540978 , lik =  73877.12 nz =  8 , nz.p =  7 
#> alpha =  0.05312949 , tau =  14.94456 , beta = 1.0682 , sigma_e =  0.01539662 , lik =  73885.7 nz =  8 , nz.p =  7 
#> alpha =  0.04386166 , tau =  15.10473 , beta = 1.03301 , sigma_e =  0.01525293 , lik =  73880.14 nz =  8 , nz.p =  7 
#> alpha =  0.04136831 , tau =  15.45396 , beta = 1.036861 , sigma_e =  0.0153237 , lik =  73877.09 nz =  8 , nz.p =  7 
#> alpha =  0.05978438 , tau =  14.9074 , beta = 1.041317 , sigma_e =  0.01511981 , lik =  73873.75 nz =  8 , nz.p =  7 
#> alpha =  0.05020981 , tau =  15.23843 , beta = 1.041378 , sigma_e =  0.01522698 , lik =  73879.43 nz =  8 , nz.p =  7 
#> alpha =  0.05579102 , tau =  14.51619 , beta = 1.04737 , sigma_e =  0.01512626 , lik =  73897.92 nz =  8 , nz.p =  7 
#> alpha =  0.06917401 , tau =  13.7077 , beta = 1.051116 , sigma_e =  0.01496476 , lik =  73904.86 nz =  8 , nz.p =  7 
#> alpha =  0.05792022 , tau =  14.05658 , beta = 1.037075 , sigma_e =  0.01488211 , lik =  73898.61 nz =  8 , nz.p =  7 
#> alpha =  0.05360181 , tau =  14.48702 , beta = 1.04011 , sigma_e =  0.01503229 , lik =  73898.21 nz =  8 , nz.p =  7 
#> alpha =  0.0792321 , tau =  13.54135 , beta = 1.071251 , sigma_e =  0.01520493 , lik =  73893.97 nz =  8 , nz.p =  7 
#> alpha =  0.06554639 , tau =  14.05948 , beta = 1.059584 , sigma_e =  0.01517417 , lik =  73900.9 nz =  8 , nz.p =  7 
#> alpha =  0.06514911 , tau =  13.54016 , beta = 1.058547 , sigma_e =  0.01503949 , lik =  73915.2 nz =  8 , nz.p =  7 
#> alpha =  0.07421103 , tau =  12.76338 , beta = 1.066187 , sigma_e =  0.01494662 , lik =  73923.02 nz =  8 , nz.p =  7 
#> alpha =  0.09203763 , tau =  12.77003 , beta = 1.074862 , sigma_e =  0.01489257 , lik =  73904.01 nz =  8 , nz.p =  7 
#> alpha =  0.07647074 , tau =  13.31748 , beta = 1.066684 , sigma_e =  0.01498185 , lik =  73913.93 nz =  8 , nz.p =  7 
#> alpha =  0.08789793 , tau =  12.32529 , beta = 1.039894 , sigma_e =  0.01459331 , lik =  73893.87 nz =  8 , nz.p =  7 
#> alpha =  0.07750289 , tau =  12.93358 , beta = 1.048667 , sigma_e =  0.01479012 , lik =  73908.46 nz =  8 , nz.p =  7 
#> alpha =  0.09058978 , tau =  12.67466 , beta = 1.078518 , sigma_e =  0.01506043 , lik =  73910.86 nz =  8 , nz.p =  7 
#> alpha =  0.08100595 , tau =  13.00685 , beta = 1.06875 , sigma_e =  0.01501565 , lik =  73915.72 nz =  8 , nz.p =  7 
#> alpha =  0.07037989 , tau =  13.5928 , beta = 1.060092 , sigma_e =  0.01505642 , lik =  73910.69 nz =  8 , nz.p =  7 
#> alpha =  0.08312964 , tau =  12.55655 , beta = 1.072947 , sigma_e =  0.01495095 , lik =  73927.37 nz =  8 , nz.p =  7 
#> alpha =  0.09113021 , tau =  12.01775 , beta = 1.083608 , sigma_e =  0.01494405 , lik =  73935.2 nz =  8 , nz.p =  7 
#> alpha =  0.07916237 , tau =  12.92284 , beta = 1.090467 , sigma_e =  0.01519027 , lik =  73925.52 nz =  8 , nz.p =  7 
#> alpha =  0.07874419 , tau =  12.92553 , beta = 1.079745 , sigma_e =  0.01508923 , lik =  73925.46 nz =  8 , nz.p =  7 
#> alpha =  0.09137264 , tau =  12.05002 , beta = 1.089912 , sigma_e =  0.01497451 , lik =  73937.39 nz =  8 , nz.p =  7 
#> alpha =  0.1041118 , tau =  11.34559 , beta = 1.104094 , sigma_e =  0.01493372 , lik =  73945.51 nz =  8 , nz.p =  7 
#> alpha =  0.09514179 , tau =  11.5355 , beta = 1.098949 , sigma_e =  0.01502969 , lik =  73941.96 nz =  8 , nz.p =  7 
#> alpha =  0.09008496 , tau =  11.95729 , beta = 1.090906 , sigma_e =  0.01501772 , lik =  73938.31 nz =  8 , nz.p =  7 
#> alpha =  0.09578366 , tau =  11.25724 , beta = 1.109636 , sigma_e =  0.01500147 , lik =  73944.11 nz =  8 , nz.p =  7 
#> alpha =  0.09185394 , tau =  11.67124 , beta = 1.099218 , sigma_e =  0.01500501 , lik =  73942.4 nz =  8 , nz.p =  7 
#> alpha =  0.1157895 , tau =  10.91001 , beta = 1.128167 , sigma_e =  0.01509285 , lik =  73949.7 nz =  8 , nz.p =  7 
#> alpha =  0.1446337 , tau =  10.08683 , beta = 1.157222 , sigma_e =  0.0151665 , lik =  73938.46 nz =  8 , nz.p =  7 
#> alpha =  0.1263802 , tau =  10.06976 , beta = 1.115695 , sigma_e =  0.01481259 , lik =  73949.58 nz =  8 , nz.p =  7 
#> alpha =  0.1124317 , tau =  10.71776 , beta = 1.110829 , sigma_e =  0.01490612 , lik =  73949.2 nz =  8 , nz.p =  7 
#> alpha =  0.1251244 , tau =  10.08864 , beta = 1.139895 , sigma_e =  0.01500353 , lik =  73958.56 nz =  8 , nz.p =  7 
#> alpha =  0.1466161 , tau =  9.243515 , beta = 1.167979 , sigma_e =  0.01503336 , lik =  73963.93 nz =  8 , nz.p =  7 
#> alpha =  0.1424811 , tau =  9.619221 , beta = 1.150388 , sigma_e =  0.01491949 , lik =  73961.65 nz =  8 , nz.p =  7 
#> alpha =  0.1287985 , tau =  10.06616 , beta = 1.138131 , sigma_e =  0.01494696 , lik =  73959.09 nz =  8 , nz.p =  7 
#> alpha =  0.1658912 , tau =  9.256089 , beta = 1.151578 , sigma_e =  0.01491484 , lik =  73942.95 nz =  8 , nz.p =  7 
#> alpha =  0.1098814 , tau =  10.71965 , beta = 1.122022 , sigma_e =  0.01497976 , lik =  73952.38 nz =  8 , nz.p =  7 
#> alpha =  0.1559497 , tau =  8.977882 , beta = 1.169095 , sigma_e =  0.01500095 , lik =  73964.54 nz =  8 , nz.p =  7 
#> alpha =  0.1908652 , tau =  7.986331 , beta = 1.199639 , sigma_e =  0.01503468 , lik =  73966.46 nz =  8 , nz.p =  7 
#> alpha =  0.1515042 , tau =  9.219615 , beta = 1.195458 , sigma_e =  0.01521392 , lik =  73958.53 nz =  8 , nz.p =  7 
#> alpha =  0.1447899 , tau =  9.425175 , beta = 1.174637 , sigma_e =  0.01511258 , lik =  73962.3 nz =  8 , nz.p =  7 
#> alpha =  0.1808265 , tau =  8.025882 , beta = 1.197285 , sigma_e =  0.01493922 , lik =  73966.04 nz =  8 , nz.p =  7 
#> alpha =  0.1617573 , tau =  8.666143 , beta = 1.180703 , sigma_e =  0.01497748 , lik =  73965.48 nz =  8 , nz.p =  7 
#> alpha =  0.2325637 , tau =  7.275483 , beta = 1.22746 , sigma_e =  0.01503569 , lik =  73948.11 nz =  8 , nz.p =  7 
#> alpha =  0.1325344 , tau =  9.729736 , beta = 1.150734 , sigma_e =  0.01499372 , lik =  73960.6 nz =  8 , nz.p =  7 
#> alpha =  0.1928136 , tau =  8.015702 , beta = 1.204253 , sigma_e =  0.01502169 , lik =  73964.42 nz =  8 , nz.p =  7 
#> alpha =  0.175564 , tau =  8.413592 , beta = 1.19159 , sigma_e =  0.01501469 , lik =  73966.22 nz =  8 , nz.p =  7 
#> alpha =  0.1949757 , tau =  7.684919 , beta = 1.223551 , sigma_e =  0.01513489 , lik =  73962.81 nz =  8 , nz.p =  7 
#> alpha =  0.1802705 , tau =  8.128576 , beta = 1.205088 , sigma_e =  0.01508075 , lik =  73965.14 nz =  8 , nz.p =  7 
#> alpha =  0.2094469 , tau =  7.392305 , beta = 1.207253 , sigma_e =  0.01492892 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2519079 , tau =  6.54674 , beta = 1.217677 , sigma_e =  0.01483793 , lik =  73964.95 nz =  8 , nz.p =  7 
#> alpha =  0.2385544 , tau =  6.893014 , beta = 1.228186 , sigma_e =  0.0149658 , lik =  73965.84 nz =  8 , nz.p =  7 
#> alpha =  0.2112205 , tau =  7.417647 , beta = 1.215007 , sigma_e =  0.01498266 , lik =  73967.13 nz =  8 , nz.p =  7 
#> alpha =  0.2067136 , tau =  7.556706 , beta = 1.19896 , sigma_e =  0.01487988 , lik =  73967.31 nz =  8 , nz.p =  7 
#> alpha =  0.1997597 , tau =  7.695787 , beta = 1.200792 , sigma_e =  0.01492984 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2144519 , tau =  7.525898 , beta = 1.207959 , sigma_e =  0.01501708 , lik =  73964.63 nz =  8 , nz.p =  7 
#> alpha =  0.1900575 , tau =  7.859101 , beta = 1.199203 , sigma_e =  0.01493453 , lik =  73967.14 nz =  8 , nz.p =  7 
#> alpha =  0.1952618 , tau =  7.839713 , beta = 1.200252 , sigma_e =  0.01498217 , lik =  73967.28 nz =  8 , nz.p =  7 
#> alpha =  0.1872715 , tau =  8.046689 , beta = 1.196424 , sigma_e =  0.01497221 , lik =  73967.23 nz =  8 , nz.p =  7 
#> alpha =  0.204546 , tau =  7.542519 , beta = 1.204044 , sigma_e =  0.01492938 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2054102 , tau =  7.555437 , beta = 1.207889 , sigma_e =  0.01495623 , lik =  73967.49 nz =  8 , nz.p =  7 
#> alpha =  0.2069755 , tau =  7.610369 , beta = 1.20444 , sigma_e =  0.0149734 , lik =  73966.87 nz =  8 , nz.p =  7 
#> alpha =  0.1941527 , tau =  7.796166 , beta = 1.200596 , sigma_e =  0.01494424 , lik =  73967.44 nz =  8 , nz.p =  7 
#> alpha =  0.2131086 , tau =  7.339507 , beta = 1.208619 , sigma_e =  0.01492455 , lik =  73967.5 nz =  8 , nz.p =  7 
#> alpha =  0.206333 , tau =  7.510252 , beta = 1.205736 , sigma_e =  0.01493645 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2089469 , tau =  7.405011 , beta = 1.207305 , sigma_e =  0.0148964 , lik =  73967.43 nz =  8 , nz.p =  7 
#> alpha =  0.2054382 , tau =  7.511373 , beta = 1.205585 , sigma_e =  0.0149178 , lik =  73967.55 nz =  8 , nz.p =  7 
#> alpha =  0.2149438 , tau =  7.33635 , beta = 1.208728 , sigma_e =  0.01492364 , lik =  73967.5 nz =  8 , nz.p =  7 
#> alpha =  0.2095461 , tau =  7.448697 , beta = 1.206808 , sigma_e =  0.01492879 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2047903 , tau =  7.527137 , beta = 1.201329 , sigma_e =  0.01490072 , lik =  73967.48 nz =  8 , nz.p =  7 
#> alpha =  0.205255 , tau =  7.548352 , beta = 1.206244 , sigma_e =  0.01494233 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2046893 , tau =  7.586182 , beta = 1.203884 , sigma_e =  0.01494893 , lik =  73967.55 nz =  8 , nz.p =  7 
#> alpha =  0.2048763 , tau =  7.56741 , beta = 1.204309 , sigma_e =  0.01494114 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2032249 , tau =  7.61036 , beta = 1.203158 , sigma_e =  0.01493214 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2055515 , tau =  7.535154 , beta = 1.205093 , sigma_e =  0.01493537 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.203989 , tau =  7.575812 , beta = 1.200185 , sigma_e =  0.01491407 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.574218 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.559085 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046147 , tau =  7.566648 , beta = 1.201995 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042059 , tau =  7.566648 , beta = 1.202403 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.203106 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.201293 , sigma_e =  0.01492348 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01493841 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01490857 , lik =  73967.55 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.438898 , sigma_e =  0.0384372 , lik =  61702.74 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.989593 , beta = 1.438898 , sigma_e =  0.0384372 , lik =  61702.24 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.971631 , beta = 1.438898 , sigma_e =  0.0384372 , lik =  61703.23 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.438898 , sigma_e =  0.0384372 , lik =  61702.62 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.438898 , sigma_e =  0.0384372 , lik =  61702.85 nz =  8 , nz.p =  7 
#> alpha =  0.2281544 , tau =  8.980607 , beta = 1.43867 , sigma_e =  0.0384372 , lik =  61702.34 nz =  8 , nz.p =  7 
#> alpha =  0.2276985 , tau =  8.980607 , beta = 1.439126 , sigma_e =  0.0384372 , lik =  61703.14 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.440066 , sigma_e =  0.0384372 , lik =  61702.67 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.437732 , sigma_e =  0.0384372 , lik =  61702.81 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.438898 , sigma_e =  0.03847566 , lik =  61681.79 nz =  8 , nz.p =  7 
#> alpha =  0.2279263 , tau =  8.980607 , beta = 1.438898 , sigma_e =  0.03839878 , lik =  61723.67 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.574767 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.559632 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046241 , tau =  7.567196 , beta = 1.202082 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042152 , tau =  7.567196 , beta = 1.202491 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.203193 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.20138 , sigma_e =  0.01492945 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01494439 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01491453 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.27426 , sigma_e =  0.02046192 , lik =  72167.24 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.019726 , beta = 1.27426 , sigma_e =  0.02046192 , lik =  72166.94 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.003702 , beta = 1.27426 , sigma_e =  0.02046192 , lik =  72167.53 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.27426 , sigma_e =  0.02046192 , lik =  72167.17 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.27426 , sigma_e =  0.02046192 , lik =  72167.3 nz =  8 , nz.p =  7 
#> alpha =  0.2121847 , tau =  8.01171 , beta = 1.274048 , sigma_e =  0.02046192 , lik =  72167.04 nz =  8 , nz.p =  7 
#> alpha =  0.2117608 , tau =  8.01171 , beta = 1.274472 , sigma_e =  0.02046192 , lik =  72167.43 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.275247 , sigma_e =  0.02046192 , lik =  72167.1 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.273275 , sigma_e =  0.02046192 , lik =  72167.38 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.27426 , sigma_e =  0.0204824 , lik =  72156.81 nz =  8 , nz.p =  7 
#> alpha =  0.2119726 , tau =  8.01171 , beta = 1.27426 , sigma_e =  0.02044147 , lik =  72177.64 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.225569 , sigma_e =  0.01658357 , lik =  73744.19 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.720273 , beta = 1.225569 , sigma_e =  0.01658357 , lik =  73744.08 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.704848 , beta = 1.225569 , sigma_e =  0.01658357 , lik =  73744.3 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.225569 , sigma_e =  0.01658357 , lik =  73744.17 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.225569 , sigma_e =  0.01658357 , lik =  73744.22 nz =  8 , nz.p =  7 
#> alpha =  0.2071139 , tau =  7.712557 , beta = 1.225362 , sigma_e =  0.01658357 , lik =  73744.12 nz =  8 , nz.p =  7 
#> alpha =  0.2067001 , tau =  7.712557 , beta = 1.225775 , sigma_e =  0.01658357 , lik =  73744.26 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.226502 , sigma_e =  0.01658357 , lik =  73744.13 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.224637 , sigma_e =  0.01658357 , lik =  73744.25 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.225569 , sigma_e =  0.01660016 , lik =  73740.09 nz =  8 , nz.p =  7 
#> alpha =  0.2069069 , tau =  7.712557 , beta = 1.225569 , sigma_e =  0.01656699 , lik =  73748.26 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.20997 , sigma_e =  0.01546163 , lik =  73941.73 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.622962 , beta = 1.20997 , sigma_e =  0.01546163 , lik =  73941.69 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.607731 , beta = 1.20997 , sigma_e =  0.01546163 , lik =  73941.76 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.20997 , sigma_e =  0.01546163 , lik =  73941.72 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.20997 , sigma_e =  0.01546163 , lik =  73941.74 nz =  8 , nz.p =  7 
#> alpha =  0.2054507 , tau =  7.615342 , beta = 1.209765 , sigma_e =  0.01546163 , lik =  73941.7 nz =  8 , nz.p =  7 
#> alpha =  0.2050402 , tau =  7.615342 , beta = 1.210175 , sigma_e =  0.01546163 , lik =  73941.75 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.210886 , sigma_e =  0.01546163 , lik =  73941.71 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.209056 , sigma_e =  0.01546163 , lik =  73941.75 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.20997 , sigma_e =  0.0154771 , lik =  73940.26 nz =  8 , nz.p =  7 
#> alpha =  0.2052453 , tau =  7.615342 , beta = 1.20997 , sigma_e =  0.01544617 , lik =  73943.15 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.583211 , beta = 1.204839 , sigma_e =  0.01510478 , lik =  73964.66 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.590798 , beta = 1.204839 , sigma_e =  0.01510478 , lik =  73964.65 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.575631 , beta = 1.204839 , sigma_e =  0.01510478 , lik =  73964.67 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.583211 , beta = 1.204839 , sigma_e =  0.01510478 , lik =  73964.66 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.583211 , beta = 1.204839 , sigma_e =  0.01510478 , lik =  73964.67 nz =  8 , nz.p =  7 
#> alpha =  0.2048992 , tau =  7.583211 , beta = 1.204634 , sigma_e =  0.01510478 , lik =  73964.65 nz =  8 , nz.p =  7 
#> alpha =  0.2044899 , tau =  7.583211 , beta = 1.205044 , sigma_e =  0.01510478 , lik =  73964.67 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.583211 , beta = 1.20393 , sigma_e =  0.01510478 , lik =  73964.66 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.583211 , beta = 1.204839 , sigma_e =  0.01511989 , lik =  73964.15 nz =  8 , nz.p =  7 
#> alpha =  0.2046945 , tau =  7.583211 , beta = 1.204839 , sigma_e =  0.01508968 , lik =  73965.14 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.576502 , beta = 1.203769 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.584082 , beta = 1.203769 , sigma_e =  0.01503112 , lik =  73966.58 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.568929 , beta = 1.203769 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.576502 , beta = 1.203769 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.204784 , tau =  7.576502 , beta = 1.203565 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.2043748 , tau =  7.576502 , beta = 1.203974 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.576502 , beta = 1.204678 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.576502 , beta = 1.202861 , sigma_e =  0.01503112 , lik =  73966.59 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.576502 , beta = 1.203769 , sigma_e =  0.01504616 , lik =  73966.28 nz =  8 , nz.p =  7 
#> alpha =  0.2045793 , tau =  7.576502 , beta = 1.203769 , sigma_e =  0.0150161 , lik =  73966.86 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.202599 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.576734 , beta = 1.202599 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.561596 , beta = 1.202599 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.202599 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.202599 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2046579 , tau =  7.569161 , beta = 1.202395 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.204249 , tau =  7.569161 , beta = 1.202804 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.203507 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.201693 , sigma_e =  0.01495088 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.202599 , sigma_e =  0.01496584 , lik =  73967.46 nz =  8 , nz.p =  7 
#> alpha =  0.2044533 , tau =  7.569161 , beta = 1.202599 , sigma_e =  0.01493594 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.575422 , beta = 1.202391 , sigma_e =  0.01493659 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.560287 , beta = 1.202391 , sigma_e =  0.01493659 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.567851 , beta = 1.202391 , sigma_e =  0.01493659 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.567851 , beta = 1.202391 , sigma_e =  0.01493659 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046353 , tau =  7.567851 , beta = 1.202186 , sigma_e =  0.01493659 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042265 , tau =  7.567851 , beta = 1.202595 , sigma_e =  0.01493659 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.567851 , beta = 1.203298 , sigma_e =  0.01493659 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.567851 , beta = 1.201484 , sigma_e =  0.01493659 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.567851 , beta = 1.202391 , sigma_e =  0.01495153 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2044308 , tau =  7.567851 , beta = 1.202391 , sigma_e =  0.01492166 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.574767 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.559632 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046241 , tau =  7.567196 , beta = 1.202082 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042153 , tau =  7.567196 , beta = 1.202491 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.203194 , sigma_e =  0.01492945 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.20138 , sigma_e =  0.01492945 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01494439 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044196 , tau =  7.567196 , beta = 1.202286 , sigma_e =  0.01491453 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.202348 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.575109 , beta = 1.202348 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.559974 , beta = 1.202348 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.202348 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.202348 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.20463 , tau =  7.567538 , beta = 1.202144 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042212 , tau =  7.567538 , beta = 1.202553 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.203256 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.201442 , sigma_e =  0.01492867 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.202348 , sigma_e =  0.01494361 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044255 , tau =  7.567538 , beta = 1.202348 , sigma_e =  0.01491375 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.202466 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.575474 , beta = 1.202466 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.560338 , beta = 1.202466 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.202466 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.202466 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046374 , tau =  7.567903 , beta = 1.202261 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042286 , tau =  7.567903 , beta = 1.20267 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.203373 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.20156 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.202466 , sigma_e =  0.01494337 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044329 , tau =  7.567903 , beta = 1.202466 , sigma_e =  0.01491351 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.202554 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.575287 , beta = 1.202554 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.560152 , beta = 1.202554 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.202554 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.202554 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042281 , tau =  7.567716 , beta = 1.202758 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.203461 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.201647 , sigma_e =  0.01492845 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.202554 , sigma_e =  0.01494339 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044324 , tau =  7.567716 , beta = 1.202554 , sigma_e =  0.01491353 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.202904 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.574538 , beta = 1.202904 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.559404 , beta = 1.202904 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.202904 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.202904 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046351 , tau =  7.566968 , beta = 1.202699 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042263 , tau =  7.566968 , beta = 1.203108 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.203811 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.201997 , sigma_e =  0.01492852 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.202904 , sigma_e =  0.01494345 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044306 , tau =  7.566968 , beta = 1.202904 , sigma_e =  0.0149136 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.20292 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.57369 , beta = 1.20292 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.558558 , beta = 1.20292 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.20292 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.20292 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046303 , tau =  7.56612 , beta = 1.202716 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042214 , tau =  7.56612 , beta = 1.203125 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.203828 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.202013 , sigma_e =  0.01492844 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.20292 , sigma_e =  0.01494337 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044258 , tau =  7.56612 , beta = 1.20292 , sigma_e =  0.01491352 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.580217 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.57264 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.57264 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.549956 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.55751 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.55751 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.57264 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.55751 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2048359 , tau =  7.565071 , beta = 1.202593 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.203705 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.20189 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.57264 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.55751 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.204114 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.202299 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.203705 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.204114 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.204819 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01491372 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.20189 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.202299 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.201189 , sigma_e =  0.01492865 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01494358 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01494358 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01494358 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01495853 , lik =  73967.5 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.57264 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.55751 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2046312 , tau =  7.565071 , beta = 1.202797 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2042223 , tau =  7.565071 , beta = 1.203206 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.20391 , sigma_e =  0.01491372 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.202095 , sigma_e =  0.01491372 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01492865 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044266 , tau =  7.565071 , beta = 1.203002 , sigma_e =  0.01489882 , lik =  73967.5 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  11.5217 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73013.93 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73931.3 nz =  8 , nz.p =  7 
#> alpha =  0.3112544 , tau =  7.566648 , beta = 1.095355 , sigma_e =  0.01492348 , lik =  73410.7 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.676079 , sigma_e =  0.01492348 , lik =  73621.19 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.02272392 , lik =  70873.82 nz =  8 , nz.p =  7 
#> alpha =  0.2418509 , tau =  8.95259 , beta = 1.330817 , sigma_e =  0.009800703 , lik =  66713.18 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  9.337057 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73739.98 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01492348 , lik =  73958.65 nz =  8 , nz.p =  7 
#> alpha =  0.2522371 , tau =  7.566648 , beta = 1.154372 , sigma_e =  0.01492348 , lik =  73862.55 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.414323 , sigma_e =  0.01492348 , lik =  73856.91 nz =  8 , nz.p =  7 
#> alpha =  0.2044102 , tau =  7.566648 , beta = 1.202199 , sigma_e =  0.01841521 , lik =  73107.54 nz =  8 , nz.p =  7 
#> alpha =  0.2223438 , tau =  8.230498 , beta = 1.263806 , sigma_e =  0.01209383 , lik =  72424 nz =  8 , nz.p =  7 
#> alpha =  0.2087532 , tau =  7.727414 , beta = 1.217118 , sigma_e =  0.01657767 , lik =  73744.38 nz =  8 , nz.p =  7 
#> alpha =  0.2242215 , tau =  6.183713 , beta = 1.270256 , sigma_e =  0.01556437 , lik =  73898.59 nz =  8 , nz.p =  7 
#> alpha =  0.2190956 , tau =  6.854716 , beta = 1.252647 , sigma_e =  0.01540161 , lik =  73944.41 nz =  8 , nz.p =  7 
#> alpha =  0.2238447 , tau =  7.122087 , beta = 1.268962 , sigma_e =  0.0136049 , lik =  73704.8 nz =  8 , nz.p =  7 
#> alpha =  0.2124279 , tau =  7.571422 , beta = 1.229742 , sigma_e =  0.01577853 , lik =  73896.39 nz =  8 , nz.p =  7 
#> alpha =  0.232146 , tau =  7.275242 , beta = 1.03495 , sigma_e =  0.01545351 , lik =  73855.12 nz =  8 , nz.p =  7 
#> alpha =  0.2110169 , tau =  7.49272 , beta = 1.307006 , sigma_e =  0.01505426 , lik =  73918.57 nz =  8 , nz.p =  7 
#> alpha =  0.1751693 , tau =  7.246726 , beta = 1.316327 , sigma_e =  0.01550753 , lik =  73903.02 nz =  8 , nz.p =  7 
#> alpha =  0.1918871 , tau =  7.325416 , beta = 1.27767 , sigma_e =  0.01535941 , lik =  73942.76 nz =  8 , nz.p =  7 
#> alpha =  0.199705 , tau =  7.1471 , beta = 1.265367 , sigma_e =  0.01451008 , lik =  73943.33 nz =  8 , nz.p =  7 
#> alpha =  0.2028125 , tau =  7.250898 , beta = 1.256482 , sigma_e =  0.0148173 , lik =  73961.44 nz =  8 , nz.p =  7 
#> alpha =  0.1978766 , tau =  7.128017 , beta = 1.17431 , sigma_e =  0.01511195 , lik =  73934.68 nz =  8 , nz.p =  7 
#> alpha =  0.2010829 , tau =  7.217494 , beta = 1.205477 , sigma_e =  0.01509751 , lik =  73958.91 nz =  8 , nz.p =  7 
#> alpha =  0.2217163 , tau =  7.247657 , beta = 1.169786 , sigma_e =  0.01471017 , lik =  73958.66 nz =  8 , nz.p =  7 
#> alpha =  0.1951081 , tau =  7.919972 , beta = 1.164248 , sigma_e =  0.0144028 , lik =  73935.06 nz =  8 , nz.p =  7 
#> alpha =  0.2128355 , tau =  7.106782 , beta = 1.229503 , sigma_e =  0.01514559 , lik =  73962.35 nz =  8 , nz.p =  7 
#> alpha =  0.2125314 , tau =  6.997089 , beta = 1.222915 , sigma_e =  0.01495234 , lik =  73957.89 nz =  8 , nz.p =  7 
#> alpha =  0.2064109 , tau =  7.420053 , beta = 1.207323 , sigma_e =  0.01493069 , lik =  73965.32 nz =  8 , nz.p =  7 
#> alpha =  0.1904153 , tau =  7.374067 , beta = 1.270302 , sigma_e =  0.01525972 , lik =  73949.67 nz =  8 , nz.p =  7 
#> alpha =  0.213439 , tau =  7.279055 , beta = 1.194885 , sigma_e =  0.01484567 , lik =  73966.22 nz =  8 , nz.p =  7 
#> alpha =  0.2150225 , tau =  7.430084 , beta = 1.230473 , sigma_e =  0.01476851 , lik =  73949.91 nz =  8 , nz.p =  7 
#> alpha =  0.2044807 , tau =  7.270064 , beta = 1.21166 , sigma_e =  0.01501458 , lik =  73965.92 nz =  8 , nz.p =  7 
#> alpha =  0.2138885 , tau =  7.403668 , beta = 1.16333 , sigma_e =  0.01512763 , lik =  73960.47 nz =  8 , nz.p =  7 
#> alpha =  0.2055265 , tau =  7.288792 , beta = 1.232556 , sigma_e =  0.01489428 , lik =  73965.81 nz =  8 , nz.p =  7 
#> alpha =  0.2009868 , tau =  7.6306 , beta = 1.190401 , sigma_e =  0.014701 , lik =  73961.84 nz =  8 , nz.p =  7 
#> alpha =  0.2098094 , tau =  7.234265 , beta = 1.219525 , sigma_e =  0.0150332 , lik =  73966.14 nz =  8 , nz.p =  7 
#> alpha =  0.2086008 , tau =  7.234686 , beta = 1.216938 , sigma_e =  0.01495346 , lik =  73965.65 nz =  8 , nz.p =  7 
#> alpha =  0.2080511 , tau =  7.280589 , beta = 1.214522 , sigma_e =  0.01494777 , lik =  73967.02 nz =  8 , nz.p =  7 
#> alpha =  0.210524 , tau =  7.361656 , beta = 1.185 , sigma_e =  0.01501153 , lik =  73965.67 nz =  8 , nz.p =  7 
#> alpha =  0.2067646 , tau =  7.30694 , beta = 1.220492 , sigma_e =  0.0149235 , lik =  73967.05 nz =  8 , nz.p =  7 
#> alpha =  0.2125429 , tau =  7.395572 , beta = 1.208864 , sigma_e =  0.01485506 , lik =  73965.69 nz =  8 , nz.p =  7 
#> alpha =  0.2064671 , tau =  7.30124 , beta = 1.210987 , sigma_e =  0.01497454 , lik =  73967.12 nz =  8 , nz.p =  7 
#> alpha =  0.2058185 , tau =  7.459615 , beta = 1.197845 , sigma_e =  0.01481347 , lik =  73966.23 nz =  8 , nz.p =  7 
#> alpha =  0.2068091 , tau =  7.402628 , beta = 1.203205 , sigma_e =  0.0148681 , lik =  73967.21 nz =  8 , nz.p =  7 
#> alpha =  0.1997809 , tau =  7.463805 , beta = 1.225472 , sigma_e =  0.01500965 , lik =  73965.15 nz =  8 , nz.p =  7 
#> alpha =  0.2099393 , tau =  7.324809 , beta = 1.202588 , sigma_e =  0.0148865 , lik =  73967.32 nz =  8 , nz.p =  7 
#> alpha =  0.2056966 , tau =  7.480328 , beta = 1.201285 , sigma_e =  0.01488267 , lik =  73967.05 nz =  8 , nz.p =  7 
#> alpha =  0.2062827 , tau =  7.429885 , beta = 1.204572 , sigma_e =  0.01489891 , lik =  73967.41 nz =  8 , nz.p =  7 
#> alpha =  0.2067833 , tau =  7.503264 , beta = 1.189198 , sigma_e =  0.01489703 , lik =  73967.14 nz =  8 , nz.p =  7 
#> alpha =  0.2072081 , tau =  7.591553 , beta = 1.189808 , sigma_e =  0.01481548 , lik =  73965.91 nz =  8 , nz.p =  7 
#> alpha =  0.2066521 , tau =  7.37276 , beta = 1.205649 , sigma_e =  0.01493461 , lik =  73967.52 nz =  8 , nz.p =  7 
#> alpha =  0.2068387 , tau =  7.335487 , beta = 1.218332 , sigma_e =  0.01490758 , lik =  73967.06 nz =  8 , nz.p =  7 
#> alpha =  0.2067972 , tau =  7.460964 , beta = 1.196394 , sigma_e =  0.01489967 , lik =  73967.43 nz =  8 , nz.p =  7 
#> alpha =  0.2068083 , tau =  7.458592 , beta = 1.201358 , sigma_e =  0.01494926 , lik =  73967.48 nz =  8 , nz.p =  7 
#> alpha =  0.2068085 , tau =  7.444562 , beta = 1.201819 , sigma_e =  0.01492893 , lik =  73967.56 nz =  8 , nz.p =  7 
#> alpha =  0.2025039 , tau =  7.586887 , beta = 1.201611 , sigma_e =  0.01494779 , lik =  73967.5 nz =  8 , nz.p =  7 
#> alpha =  0.2043377 , tau =  7.520501 , beta = 1.201873 , sigma_e =  0.01493244 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2053141 , tau =  7.515947 , beta = 1.19861 , sigma_e =  0.01494877 , lik =  73967.38 nz =  8 , nz.p =  7 
#> alpha =  0.2060401 , tau =  7.451308 , beta = 1.203077 , sigma_e =  0.01491136 , lik =  73967.53 nz =  8 , nz.p =  7 
#> alpha =  0.204503 , tau =  7.480763 , beta = 1.209479 , sigma_e =  0.01495271 , lik =  73967.46 nz =  8 , nz.p =  7 
#> alpha =  0.2050742 , tau =  7.475808 , beta = 1.206198 , sigma_e =  0.01493943 , lik =  73967.55 nz =  8 , nz.p =  7 
#> alpha =  0.2040201 , tau =  7.612406 , beta = 1.200426 , sigma_e =  0.01491964 , lik =  73967.51 nz =  8 , nz.p =  7 
#> alpha =  0.205991 , tau =  7.431955 , beta = 1.20434 , sigma_e =  0.01493087 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2046066 , tau =  7.524329 , beta = 1.203493 , sigma_e =  0.01495072 , lik =  73967.53 nz =  8 , nz.p =  7 
#> alpha =  0.2056808 , tau =  7.469496 , beta = 1.203182 , sigma_e =  0.01492119 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2058133 , tau =  7.497136 , beta = 1.199179 , sigma_e =  0.01491535 , lik =  73967.54 nz =  8 , nz.p =  7 
#> alpha =  0.2052587 , tau =  7.481135 , beta = 1.20444 , sigma_e =  0.01493341 , lik =  73967.57 nz =  8 , nz.p =  7 
#> alpha =  0.2034743 , tau =  7.543376 , beta = 1.20458 , sigma_e =  0.01492763 , lik =  73967.55 nz =  8 , nz.p =  7 
#> alpha =  0.2059698 , tau =  7.469143 , beta = 1.202514 , sigma_e =  0.0149286 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2042734 , tau =  7.571273 , beta = 1.201345 , sigma_e =  0.01492478 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2047015 , tau =  7.5362 , beta = 1.202093 , sigma_e =  0.0149263 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2060727 , tau =  7.488372 , beta = 1.203899 , sigma_e =  0.01492075 , lik =  73967.55 nz =  8 , nz.p =  7 
#> alpha =  0.2047701 , tau =  7.512456 , beta = 1.202379 , sigma_e =  0.01492952 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.204364 , tau =  7.556821 , beta = 1.202268 , sigma_e =  0.01493534 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2046924 , tau =  7.534895 , beta = 1.202496 , sigma_e =  0.0149318 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.204558 , tau =  7.566706 , beta = 1.20024 , sigma_e =  0.01492248 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.204733 , tau =  7.545222 , beta = 1.201288 , sigma_e =  0.01492521 , lik =  73967.58 nz =  8 , nz.p =  7 
#> alpha =  0.2033612 , tau =  7.609639 , beta = 1.201663 , sigma_e =  0.01492592 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2020693 , tau =  7.680875 , beta = 1.201229 , sigma_e =  0.01492458 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2039886 , tau =  7.604764 , beta = 1.201517 , sigma_e =  0.01492357 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.203599 , tau =  7.651342 , beta = 1.201086 , sigma_e =  0.01492059 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2036167 , tau =  7.626923 , beta = 1.2014 , sigma_e =  0.0149245 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2030764 , tau =  7.672693 , beta = 1.201052 , sigma_e =  0.01492359 , lik =  73967.59 nz =  8 , nz.p =  7 
#> alpha =  0.2031406 , tau =  7.650691 , beta = 1.202247 , sigma_e =  0.01492531 , lik =  73967.6 nz =  8 , nz.p =  7 
#> alpha =  0.202349 , tau =  7.703977 , beta = 1.202722 , sigma_e =  0.01492535 , lik =  73967.61 nz =  8 , nz.p =  7 
#> alpha =  0.2022472 , tau =  7.729487 , beta = 1.20113 , sigma_e =  0.01491614 , lik =  73967.6 nz =  8 , nz.p =  7 
#> alpha =  0.2028558 , tau =  7.680372 , beta = 1.201473 , sigma_e =  0.01492006 , lik =  73967.6 nz =  8 , nz.p =  7 
#> alpha =  0.201909 , tau =  7.743117 , beta = 1.201135 , sigma_e =  0.01492308 , lik =  73967.61 nz =  8 , nz.p =  7 
#> alpha =  0.2006699 , tau =  7.832889 , beta = 1.200597 , sigma_e =  0.01492288 , lik =  73967.62 nz =  8 , nz.p =  7 
#> alpha =  0.2018718 , tau =  7.788941 , beta = 1.20125 , sigma_e =  0.01491943 , lik =  73967.62 nz =  8 , nz.p =  7 
#> alpha =  0.2022432 , tau =  7.743723 , beta = 1.201354 , sigma_e =  0.01492105 , lik =  73967.61 nz =  8 , nz.p =  7 
#> alpha =  0.2009502 , tau =  7.801914 , beta = 1.201886 , sigma_e =  0.0149243 , lik =  73967.62 nz =  8 , nz.p =  7 
#> alpha =  0.2016091 , tau =  7.763996 , beta = 1.201689 , sigma_e =  0.01492337 , lik =  73967.62 nz =  8 , nz.p =  7 
#> alpha =  0.199876 , tau =  7.898237 , beta = 1.201759 , sigma_e =  0.01492031 , lik =  73967.64 nz =  8 , nz.p =  7 
#> alpha =  0.1980316 , tau =  8.037493 , beta = 1.201917 , sigma_e =  0.01491822 , lik =  73967.66 nz =  8 , nz.p =  7 
#> alpha =  0.1987035 , tau =  7.987163 , beta = 1.201866 , sigma_e =  0.01492402 , lik =  73967.66 nz =  8 , nz.p =  7 
#> alpha =  0.1966593 , tau =  8.145124 , beta = 1.202037 , sigma_e =  0.014926 , lik =  73967.69 nz =  8 , nz.p =  7 
#> alpha =  0.1969414 , tau =  8.142034 , beta = 1.200349 , sigma_e =  0.01491898 , lik =  73967.67 nz =  8 , nz.p =  7 
#> alpha =  0.1982796 , tau =  8.030239 , beta = 1.20095 , sigma_e =  0.01492057 , lik =  73967.66 nz =  8 , nz.p =  7 
#> alpha =  0.1967203 , tau =  8.178233 , beta = 1.200578 , sigma_e =  0.01491791 , lik =  73967.66 nz =  8 , nz.p =  7 
#> alpha =  0.1977694 , tau =  8.082485 , beta = 1.20091 , sigma_e =  0.01491951 , lik =  73967.66 nz =  8 , nz.p =  7 
#> alpha =  0.1942205 , tau =  8.313977 , beta = 1.201025 , sigma_e =  0.0149228 , lik =  73967.7 nz =  8 , nz.p =  7 
#> alpha =  0.1905043 , tau =  8.589621 , beta = 1.200829 , sigma_e =  0.01492449 , lik =  73967.69 nz =  8 , nz.p =  7 
#> alpha =  0.1928474 , tau =  8.466814 , beta = 1.201839 , sigma_e =  0.01491932 , lik =  73967.69 nz =  8 , nz.p =  7 
#> alpha =  0.1947739 , tau =  8.303678 , beta = 1.201553 , sigma_e =  0.01492021 , lik =  73967.69 nz =  8 , nz.p =  7 
#> alpha =  0.1941244 , tau =  8.35952 , beta = 1.200429 , sigma_e =  0.01492478 , lik =  73967.7 nz =  8 , nz.p =  7 
#> alpha =  0.1921999 , tau =  8.52534 , beta = 1.199668 , sigma_e =  0.01492806 , lik =  73967.7 nz =  8 , nz.p =  7 
#> alpha =  0.1929408 , tau =  8.425809 , beta = 1.201228 , sigma_e =  0.01492561 , lik =  73967.72 nz =  8 , nz.p =  7 
#> alpha =  0.1905709 , tau =  8.602902 , beta = 1.201352 , sigma_e =  0.01492866 , lik =  73967.73 nz =  8 , nz.p =  7 
#> alpha =  0.1912203 , tau =  8.550439 , beta = 1.202186 , sigma_e =  0.01493 , lik =  73967.71 nz =  8 , nz.p =  7 
#> alpha =  0.1926348 , tau =  8.446456 , beta = 1.201741 , sigma_e =  0.01492725 , lik =  73967.72 nz =  8 , nz.p =  7 
#> alpha =  0.1899173 , tau =  8.672301 , beta = 1.200373 , sigma_e =  0.01492348 , lik =  73967.67 nz =  8 , nz.p =  7 
#> alpha =  0.1949517 , tau =  8.273835 , beta = 1.201636 , sigma_e =  0.01492537 , lik =  73967.71 nz =  8 , nz.p =  7 
#> alpha =  0.1918256 , tau =  8.494467 , beta = 1.200922 , sigma_e =  0.01493134 , lik =  73967.73 nz =  8 , nz.p =  7 
#> alpha =  0.1925585 , tau =  8.446362 , beta = 1.201083 , sigma_e =  0.01492855 , lik =  73967.72 nz =  8 , nz.p =  7 
#> alpha =  0.1914199 , tau =  8.557142 , beta = 1.201408 , sigma_e =  0.01493215 , lik =  73967.73 nz =  8 , nz.p =  7 
#> alpha =  0.1900347 , tau =  8.681379 , beta = 1.201587 , sigma_e =  0.01493683 , lik =  73967.74 nz =  8 , nz.p =  7 
#> alpha =  0.1930572 , tau =  8.428801 , beta = 1.200944 , sigma_e =  0.01492733 , lik =  73967.72 nz =  8 , nz.p =  7 
#> alpha =  0.1883474 , tau =  8.794622 , beta = 1.200946 , sigma_e =  0.01493519 , lik =  73967.7 nz =  8 , nz.p =  7 
#> alpha =  0.1932792 , tau =  8.401067 , beta = 1.201479 , sigma_e =  0.01492783 , lik =  73967.72 nz =  8 , nz.p =  7 
#> alpha =  0.1908676 , tau =  8.596338 , beta = 1.200778 , sigma_e =  0.01493355 , lik =  73967.74 nz =  8 , nz.p =  7 
#> alpha =  0.1899901 , tau =  8.672273 , beta = 1.200293 , sigma_e =  0.0149367 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1892339 , tau =  8.713027 , beta = 1.201299 , sigma_e =  0.0149372 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1901826 , tau =  8.641085 , beta = 1.201216 , sigma_e =  0.01493474 , lik =  73967.74 nz =  8 , nz.p =  7 
#> alpha =  0.187424 , tau =  8.870218 , beta = 1.200672 , sigma_e =  0.01494047 , lik =  73967.71 nz =  8 , nz.p =  7 
#> alpha =  0.1917985 , tau =  8.515976 , beta = 1.20129 , sigma_e =  0.01493098 , lik =  73967.74 nz =  8 , nz.p =  7 
#> alpha =  0.1888336 , tau =  8.781565 , beta = 1.201398 , sigma_e =  0.01493681 , lik =  73967.74 nz =  8 , nz.p =  7 
#> alpha =  0.1895772 , tau =  8.708893 , beta = 1.201283 , sigma_e =  0.01493544 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1896798 , tau =  8.713452 , beta = 1.200951 , sigma_e =  0.01494221 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1899022 , tau =  8.685682 , beta = 1.201052 , sigma_e =  0.01493882 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1900728 , tau =  8.647432 , beta = 1.200463 , sigma_e =  0.01493618 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1900919 , tau =  8.630509 , beta = 1.199902 , sigma_e =  0.01493586 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1876528 , tau =  8.862623 , beta = 1.200187 , sigma_e =  0.01494398 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1886807 , tau =  8.77466 , beta = 1.200468 , sigma_e =  0.01494073 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1894919 , tau =  8.692418 , beta = 1.199884 , sigma_e =  0.01494164 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1895133 , tau =  8.696534 , beta = 1.200234 , sigma_e =  0.01494009 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1899477 , tau =  8.681712 , beta = 1.199441 , sigma_e =  0.01494103 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1903056 , tau =  8.666096 , beta = 1.198511 , sigma_e =  0.01494295 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1896083 , tau =  8.668625 , beta = 1.199186 , sigma_e =  0.01493556 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1896262 , tau =  8.67981 , beta = 1.199627 , sigma_e =  0.01493722 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1891535 , tau =  8.712814 , beta = 1.199577 , sigma_e =  0.01494128 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1887365 , tau =  8.733156 , beta = 1.199219 , sigma_e =  0.01494357 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1904886 , tau =  8.594834 , beta = 1.198897 , sigma_e =  0.01493838 , lik =  73967.75 nz =  8 , nz.p =  7 
#> alpha =  0.1891311 , tau =  8.729354 , beta = 1.200077 , sigma_e =  0.01494014 , lik =  73967.76 nz =  8 , nz.p =  7 
#> alpha =  0.1894987 , tau =  8.685122 , beta = 1.199074 , sigma_e =  0.01493904 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1895024 , tau =  8.687973 , beta = 1.199364 , sigma_e =  0.0149393 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1886873 , tau =  8.774827 , beta = 1.199188 , sigma_e =  0.01494465 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1890375 , tau =  8.738523 , beta = 1.199367 , sigma_e =  0.01494245 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1896083 , tau =  8.679111 , beta = 1.198731 , sigma_e =  0.01494128 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1894889 , tau =  8.691644 , beta = 1.199067 , sigma_e =  0.014941 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1886582 , tau =  8.725693 , beta = 1.199081 , sigma_e =  0.0149405 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1880167 , tau =  8.747766 , beta = 1.198899 , sigma_e =  0.01494023 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1883349 , tau =  8.754872 , beta = 1.198605 , sigma_e =  0.01494552 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1886569 , tau =  8.736046 , beta = 1.198861 , sigma_e =  0.01494344 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1891984 , tau =  8.703854 , beta = 1.198937 , sigma_e =  0.01494154 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1886057 , tau =  8.704971 , beta = 1.198746 , sigma_e =  0.01494078 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1887136 , tau =  8.713347 , beta = 1.198901 , sigma_e =  0.0149412 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1878294 , tau =  8.765839 , beta = 1.198561 , sigma_e =  0.01494469 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1882463 , tau =  8.746307 , beta = 1.198763 , sigma_e =  0.01494334 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1875855 , tau =  8.774717 , beta = 1.198837 , sigma_e =  0.01494371 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1867842 , tau =  8.810365 , beta = 1.198784 , sigma_e =  0.0149448 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.187264 , tau =  8.776121 , beta = 1.198384 , sigma_e =  0.01494218 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1876311 , tau =  8.76536 , beta = 1.198593 , sigma_e =  0.01494252 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.186935 , tau =  8.78498 , beta = 1.198633 , sigma_e =  0.01494194 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.187364 , tau =  8.772721 , beta = 1.198691 , sigma_e =  0.01494231 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.188075 , tau =  8.744038 , beta = 1.198799 , sigma_e =  0.01494202 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1871465 , tau =  8.775098 , beta = 1.198922 , sigma_e =  0.01493991 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.187317 , tau =  8.772782 , beta = 1.198832 , sigma_e =  0.01494111 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1865324 , tau =  8.809463 , beta = 1.198732 , sigma_e =  0.01494174 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1857659 , tau =  8.842358 , beta = 1.198696 , sigma_e =  0.0149416 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1856935 , tau =  8.843672 , beta = 1.198549 , sigma_e =  0.01494408 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1862716 , tau =  8.819597 , beta = 1.198638 , sigma_e =  0.01494312 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1853042 , tau =  8.857379 , beta = 1.198835 , sigma_e =  0.01494241 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1841516 , tau =  8.90375 , beta = 1.198948 , sigma_e =  0.01494235 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1848815 , tau =  8.885196 , beta = 1.198926 , sigma_e =  0.01494316 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1853927 , tau =  8.860035 , beta = 1.198855 , sigma_e =  0.01494285 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1837753 , tau =  8.939678 , beta = 1.19863 , sigma_e =  0.01494648 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1846123 , tau =  8.898246 , beta = 1.198708 , sigma_e =  0.01494484 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1829394 , tau =  8.955939 , beta = 1.198703 , sigma_e =  0.01494227 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1810468 , tau =  9.029626 , beta = 1.198639 , sigma_e =  0.01494101 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.182158 , tau =  8.996783 , beta = 1.198983 , sigma_e =  0.01494176 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1830355 , tau =  8.958259 , beta = 1.19888 , sigma_e =  0.01494234 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1806642 , tau =  9.060652 , beta = 1.198934 , sigma_e =  0.0149443 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1781662 , tau =  9.171811 , beta = 1.19901 , sigma_e =  0.01494565 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1803708 , tau =  9.054385 , beta = 1.199183 , sigma_e =  0.01493909 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.181216 , tau =  9.025571 , beta = 1.19905 , sigma_e =  0.01494093 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1775289 , tau =  9.178899 , beta = 1.198937 , sigma_e =  0.01494078 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1739631 , tau =  9.329372 , beta = 1.198856 , sigma_e =  0.01493959 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1742207 , tau =  9.332467 , beta = 1.19885 , sigma_e =  0.01494049 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.176652 , tau =  9.22339 , beta = 1.198915 , sigma_e =  0.01494095 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1739783 , tau =  9.328321 , beta = 1.198814 , sigma_e =  0.01494075 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1700272 , tau =  9.498643 , beta = 1.19862 , sigma_e =  0.01494025 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1722796 , tau =  9.416173 , beta = 1.199209 , sigma_e =  0.01494141 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1744308 , tau =  9.318012 , beta = 1.1991 , sigma_e =  0.01494131 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1706245 , tau =  9.498848 , beta = 1.198602 , sigma_e =  0.01494421 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1778833 , tau =  9.163512 , beta = 1.199078 , sigma_e =  0.01494037 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.172626 , tau =  9.373841 , beta = 1.198871 , sigma_e =  0.01493554 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1739947 , tau =  9.32292 , beta = 1.198919 , sigma_e =  0.01493807 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1713787 , tau =  9.447278 , beta = 1.198721 , sigma_e =  0.0149399 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1762343 , tau =  9.233645 , beta = 1.199006 , sigma_e =  0.01494025 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1729822 , tau =  9.375523 , beta = 1.198828 , sigma_e =  0.01494002 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1754156 , tau =  9.268912 , beta = 1.198966 , sigma_e =  0.01494019 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1751648 , tau =  9.271021 , beta = 1.198692 , sigma_e =  0.01493852 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.174981 , tau =  9.282747 , beta = 1.198794 , sigma_e =  0.01493922 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.172306 , tau =  9.390196 , beta = 1.198805 , sigma_e =  0.01493818 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1755553 , tau =  9.264812 , beta = 1.198895 , sigma_e =  0.01494026 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1735757 , tau =  9.342425 , beta = 1.198743 , sigma_e =  0.01493896 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1740339 , tau =  9.323992 , beta = 1.1988 , sigma_e =  0.01493927 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.174869 , tau =  9.287073 , beta = 1.198812 , sigma_e =  0.01493931 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1753238 , tau =  9.265995 , beta = 1.198789 , sigma_e =  0.01493917 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1738072 , tau =  9.335459 , beta = 1.19888 , sigma_e =  0.01493973 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1732233 , tau =  9.361928 , beta = 1.19892 , sigma_e =  0.01493998 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1747171 , tau =  9.300231 , beta = 1.198741 , sigma_e =  0.01494154 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1727612 , tau =  9.37514 , beta = 1.198711 , sigma_e =  0.01493951 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1734555 , tau =  9.347436 , beta = 1.19876 , sigma_e =  0.0149397 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.324529 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.315209 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.315209 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.287305 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.296597 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.296597 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.85 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.315209 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.296597 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1748857 , tau =  9.305898 , beta = 1.198436 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.199484 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.197738 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01492574 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.315209 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.296597 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.19896 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1741875 , tau =  9.305898 , beta = 1.199134 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.199834 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.198087 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.19896 , sigma_e =  0.01492574 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.8 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.199484 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.199834 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.200534 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01492574 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.197738 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.198087 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.19704 , sigma_e =  0.01494067 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01492574 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.77 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01495562 , lik =  73967.79 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01495562 , lik =  73967.78 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01497058 , lik =  73967.7 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.315209 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.81 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.296597 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.85 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01492574 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1747109 , tau =  9.305898 , beta = 1.198611 , sigma_e =  0.01492574 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1743618 , tau =  9.305898 , beta = 1.19896 , sigma_e =  0.01492574 , lik =  73967.84 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.199659 , sigma_e =  0.01492574 , lik =  73967.82 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.197912 , sigma_e =  0.01492574 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01494067 , lik =  73967.83 nz =  8 , nz.p =  7 
#> alpha =  0.1745362 , tau =  9.305898 , beta = 1.198785 , sigma_e =  0.01491082 , lik =  73967.79 nz =  8 , nz.p =  7
#> Warning in rspde_lme(y ~ -1, loc = "loc", repl = "rep", data = data, model =
#> op, : All optimization methods failed to provide a numerically
#> positive-definite Hessian. The optimization method with largest likelihood was
#> chosen. You can try to obtain a positive-definite Hessian by setting
#> 'improve_hessian' to TRUE.

# Compare estimated and true parameter values
rbind(c(fit$coeff$random_effects[c("alpha", "beta", "tau", "kappa")], fit$coeff$measurement_error), 
      c(alpha, beta, tau, kappa, sigma.e))
#>          alpha     beta      tau    kappa   std. dev
#> [1,] 0.1745362 1.198785 9.305898 22.09978 0.01494067
#> [2,] 0.3000000 1.200000 7.000000 15.00000 0.01500000
```

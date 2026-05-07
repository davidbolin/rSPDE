# inlabru implementation of the rational SPDE approach

## Introduction

In this vignette we will present the [`inlabru`](http://inlabru.org/)
implementation of the covariance-based rational SPDE approach. For
further technical details on the covariance-based approach, see the
[Rational approximation with the `rSPDE`
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette and [Bolin et al.
(2023)](https://doi.org/10.1080/10618600.2023.2231051).

We begin by providing a step-by-step illustration on how to use our
implementation. To this end we will consider a real world data set that
consists of precipitation measurements from the Paraná region in Brazil.

After the initial model fitting, we will show how to change some
parameters of the model. In the end, we will also provide an example in
which we have replicates.

The examples in this vignette are the same as those in the [R-INLA
implementation of the rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md)
vignette. As in that case, it is important to mention that one can
improve the performance by using the PARDISO solver. Please, go to
<https://www.pardiso-project.org/r-inla/#license> to apply for a
license. Also, use
[`inla.pardiso()`](https://rdrr.io/pkg/INLA/man/pardiso.html) for
instructions on how to enable the PARDISO sparse library.

## An example with real data

To illustrate our implementation of `rSPDE` in
[`inlabru`](http://inlabru.org/) we will consider a dataset available in
[`R-INLA`](https://www.r-inla.org). This data has also been used to
illustrate the SPDE approach, see for instance the book [Advanced
Spatial Modeling with Stochastic Partial Differential Equations Using R
and
INLA](https://www.routledge.com/Advanced-Spatial-Modeling-with-Stochastic-Partial-Differential-Equations/Krainski-Gomez-Rubio-Bakka-Lenzi-Castro-Camilo-Simpson-Lindgren-Rue/p/book/9780367570644)
and also the vignette [Spatial Statistics using R-INLA and Gaussian
Markov random
fields](https://sites.stat.washington.edu/peter/591/INLA.html). See also
[Lindgren et al.
(2011)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2011.00777.x)
for theoretical details on the standard SPDE approach.

The data consist of precipitation measurements from the Paraná region in
Brazil and were provided by the Brazilian National Water Agency. The
data were collected at 616 gauge stations in Paraná state, south of
Brazil, for each day in 2011.

### An rSPDE model for precipitation

We will follow the vignette [Spatial Statistics using R-INLA and
Gaussian Markov random
fields](https://sites.stat.washington.edu/peter/591/INLA.html). As
precipitation data are always positive, we will assume it is Gamma
distributed. [`R-INLA`](https://www.r-inla.org) uses the following
parameterization of the Gamma distribution,
``` math
\Gamma(\mu, \phi): \pi (y) = \frac{1}{\Gamma(\phi)} \left(\frac{\phi}{\mu}\right)^{\phi} y^{\phi - 1} \exp\left(-\frac{\phi y}{\mu}\right) .
```
In this parameterization, the distribution has expected value
$`E(x) = \mu`$ and variance $`V(x) = \mu^2/\phi`$, where $`1/\phi`$ is a
dispersion parameter.

In this example $`\mu`$ will be modelled using a stochastic model that
includes both covariates and spatial structure, resulting in the latent
Gaussian model for the precipitation measurements
``` math
\begin{align} y_i\mid \mu(s_i), \theta &\sim \Gamma(\mu(s_i),c\phi)\\ \log (\mu(s)) &= \eta(s) = \sum_k f_k(c_k(s))+u(s)\\ \theta &\sim \pi(\theta) \end{align},
```

where $`y_i`$ denotes the measurement taken at location $`s_i`$,
$`c_k(s)`$ are covariates, $`u(s)`$ is a mean-zero Gaussian Matérn
field, and $`\theta`$ is a vector containing all parameters of the
model, including smoothness of the field. That is, by using the `rSPDE`
model we will also be able to estimate the smoothness of the latent
field.

### Examining the data

We will be using [`inlabru`](http://inlabru.org/). The `inlabru` package
is available on CRAN and also on
[GitHub](https://github.com/inlabru-org/inlabru).

We begin by loading some libraries we need to get the data and build the
plots.

``` r

library(ggplot2)
library(INLA)
library(inlabru)
library(splancs)
library(viridis)
```

Let us load the data and the border of the region

``` r

data(PRprec)
data(PRborder)
```

The data frame contains daily measurements at 616 stations for the year
2011, as well as coordinates and altitude information for the
measurement stations. We will not analyze the full spatio-temporal data
set, but instead look at the total precipitation in January, which we
calculate as

``` r

Y <- rowMeans(PRprec[, 3 + 1:31])
```

In the next snippet of code, we extract the coordinates and altitudes
and remove the locations with missing values.

``` r

ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])
alt <- PRprec$Altitude[ind]
```

Let us build a plot for the precipitations:

``` r

ggplot() +
  geom_point(aes(
    x = coords[, 1], y = coords[, 2],
    colour = Y
  ), size = 2, alpha = 1) +
  geom_path(aes(x = PRborder[, 1], y = PRborder[, 2])) +
  geom_path(aes(x = PRborder[1034:1078, 1], y = PRborder[
    1034:1078,
    2
  ]), colour = "red") + 
  scale_color_viridis()
```

![](rspde_inlabru_files/figure-html/plot_precipitations-1.png)

The red line in the figure shows the coast line, and we expect the
distance to the coast to be a good covariate for precipitation.

This covariate is not available, so let us calculate it for each
observation location:

``` r

seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)
```

Now, let us plot the precipitation as a function of the possible
covariates:

``` r

par(mfrow = c(2, 2))
plot(coords[, 1], Y, cex = 0.5, xlab = "Longitude")
plot(coords[, 2], Y, cex = 0.5, xlab = "Latitude")
plot(seaDist, Y, cex = 0.5, xlab = "Distance to sea")
plot(alt, Y, cex = 0.5, xlab = "Altitude")
```

![](rspde_inlabru_files/figure-html/plot_prec_as_func-1.png)

``` r

par(mfrow = c(1, 1))
```

### Creating the rSPDE model

To use the [`inlabru`](http://inlabru.org/) implementation of the
`rSPDE` model we need to load the functions:

``` r

library(rSPDE)
```

To create a `rSPDE` model, one would the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function in a similar fashion as one would use the
[`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html)
function.

#### Mesh

We can use
[`fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
function from the `fmesher` package for creating the mesh. Let us create
a mesh which is based on a non-convex hull to avoid adding many small
triangles outside the domain of interest:

``` r

library(fmesher)

prdomain <- fm_nonconvex_hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)
plot(prmesh, asp = 1, main = "")
lines(PRborder, col = 3)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")
```

![](rspde_inlabru_files/figure-html/mesh_creation-1.png)

#### Setting up the data frame

In place of a `inla.stack`, we can set up a
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) to use
[`inlabru`](http://inlabru.org/). We refer the reader to vignettes in
<https://inlabru-org.github.io/inlabru/index.html> for further details.

``` r

library(sf)
```

    ## Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

``` r

prdata <- data.frame(long = coords[,1], lat = coords[,2], 
                        seaDist = inla.group(seaDist), y = Y)
prdata <- st_as_sf(prdata, coords = c("long", "lat"), crs = 4326)
```

#### Setting up the rSPDE model

To set up an `rSPDE`model, all we need is the mesh. By default it will
assume that we want to estimate the smoothness parameter $`\nu`$ and to
do a covariance-based rational approximation of order 2.

Later in this vignette we will also see other options for setting up
`rSPDE` models such as keeping the smoothness parameter fixed and/or
increasing the order of the covariance-based rational approximation.

Therefore, to set up a model all we have to do is use the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function:

``` r

rspde_model <- rspde.matern(mesh = prmesh)
```

Notice that this function is very reminiscent of
[`R-INLA`](https://www.r-inla.org)’s
[`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html)
function.

We will assume the following linkage between model components and
observations
``` math
\eta(s) \sim A x(s) + A \text{ Intercept} + \text{seaDist}.
```
$`\eta(s)`$ will then be used in the observation-likelihood,
``` math
y_i\mid \eta(s_i),\theta \sim \Gamma(\exp(\eta (s_i)), c\phi).
```

### Model fitting

We will build a model using the distance to the sea $`x_i`$ as a
covariate through an improper CAR(1) model with
$`\beta_{ij}=1(i\sim j)`$, which [`R-INLA`](https://www.r-inla.org)
calls a random walk of order 1. We will fit it in `inlabru`’s style:

``` r

cmp <- y ~ Intercept(1) + distSea(seaDist, model="rw1") +
field(geometry, model = rspde_model)
```

To fit the model we simply use the
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)
function:

``` r

rspde_fit <- bru(cmp, data = prdata,
  family = "Gamma",
  options = list(
    control.inla = list(int.strategy = "eb"),
    verbose = FALSE,
    num.threads = "1:1")
)
```

### inlabru results

We can look at some summaries of the posterior distributions for the
parameters, for example the fixed effects (i.e. the intercept) and the
hyper-parameters (i.e. dispersion in the gamma likelihood, the precision
of the RW1, and the parameters of the spatial field):

``` r

summary(rspde_fit)
```

    ## inlabru version: 2.14.1 
    ## INLA version: 26.05.02 
    ## Latent components:
    ## Intercept: main = linear(1)
    ## distSea: main = rw1(seaDist)
    ## field: main = cgeneric(geometry)
    ## Observation models:
    ##   Model tag: <No tag>
    ##     Family: 'Gamma'
    ##     Data class: 'sf', 'data.frame'
    ##     Response class: 'numeric'
    ##     Predictor: y ~ Intercept + distSea + field
    ##     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
    ##     Used components: effect[Intercept, distSea, field], latent[] 
    ## Time used:
    ##     Pre = 0.279, Running = 6.37, Post = 0.0841, Total = 6.74 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.942 0.041      1.862    1.942      2.021 1.942   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     distSea RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.46    1.040      12.50
    ## Precision for distSea                          7137.59 3796.179    2433.36
    ## Theta1 for field                                 -4.62    3.918     -13.76
    ## Theta2 for field                                  2.15    0.851       0.86
    ## Theta3 for field                                  2.78    3.406      -2.15
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations    14.43      16.59   14.388
    ## Precision for distSea                           6283.39   16902.74 4899.309
    ## Theta1 for field                                  -4.07       1.04   -1.305
    ## Theta2 for field                                   2.04       4.11    1.496
    ## Theta3 for field                                   2.31      10.72   -0.093
    ## 
    ## Marginal log-Likelihood:  -1254.65 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

Let $`\theta_1 = \textrm{Theta1}`$, $`\theta_2=\textrm{Theta2}`$ and
$`\theta_3=\textrm{Theta3}`$. In terms of the SPDE
``` math
(\kappa^2 I - \Delta)^{\alpha/2}(\tau u) = \mathcal{W},
```
where $`\alpha = \nu + d/2`$, we have that
``` math
\tau = \exp(\theta_1),\quad \kappa = \exp(\theta_2), 
```
and by default
``` math
\nu = 4\Big(\frac{\exp(\theta_3)}{1+\exp(\theta_3)}\Big).
```
The number 4 comes from the upper bound for $`\nu`$, which is discussed
in [R-INLA implementation of the rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md)
vignette.

In general, we have
``` math
\nu = \nu_{UB}\Big(\frac{\exp(\theta_3)}{1+\exp(\theta_3)}\Big),
```
where $`\nu_{UB}`$ is the value of the upper bound for the smoothness
parameter $`\nu`$.

Another choice for prior for $`\nu`$ is a truncated lognormal
distribution and is also discussed in [R-INLA implementation of the
rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md)
vignette.

### inlabru results in the original scale

We can obtain outputs with respect to parameters in the original scale
by using the function
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.md):

``` r

result_fit <- rspde.result(rspde_fit, "field", 
                rspde_model)
```

    ## Warning in rspde.result(rspde_fit, "field", rspde_model): the mean or mode of
    ## nu is very close to nu.upper.bound, please consider increasing nu.upper.bound,
    ## and refitting the model.

``` r

summary(result_fit)
```

    ##            mean        sd  0.025quant 0.5quant 0.975quant        mode
    ## tau    0.344287  1.086370 1.02715e-06 0.019329    2.86476 6.52182e-10
    ## kappa 13.142900 17.887100 2.38494e+00 7.495160   59.43760 3.69218e+00
    ## nu     1.496280  0.588138 2.10226e-01 1.799590    1.99994 1.99999e+00

We can also plot the posterior densities. To this end we will use the
[`gg_df()`](https://davidbolin.github.io/rSPDE/reference/gg_df.md)
function, which creates `ggplot2` user-friendly data frames:

``` r

posterior_df_fit <- gg_df(result_fit)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](rspde_inlabru_files/figure-html/plot_post-1.png)

We can also obtain the summary on a different parameterization by
setting the `parameterization` argument on the
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.md)
function:

``` r

result_fit_matern <- rspde.result(rspde_fit, "field", 
                rspde_model, parameterization = "matern")
```

    ## Warning in rspde.result(rspde_fit, "field", rspde_model, parameterization =
    ## "matern"): the mean or mode of nu is very close to nu.upper.bound, please
    ## consider increasing nu.upper.bound, and refitting the model.

``` r

summary(result_fit_matern)
```

    ##              mean        sd 0.025quant 0.5quant 0.975quant      mode
    ## std.dev 11.877100 99.445100 -0.1738420 2.714810  68.252200 -0.297612
    ## range    0.415181  0.249353  0.0294138 0.401855   0.953867  0.420189
    ## nu       1.496280  0.588138  0.2102260 1.799590   1.999940  1.999990

In a similar manner, we can obtain posterior plots on the `matern`
parameterization:

``` r

posterior_df_fit_matern <- gg_df(result_fit_matern)

ggplot(posterior_df_fit_matern) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](rspde_inlabru_files/figure-html/plot_post_matern-1.png)

### Predictions

Let us now obtain predictions (i.e. do kriging) of the expected
precipitation on a dense grid in the region.

We begin by creating the grid in which we want to do the predictions. To
this end, we can use the
[`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.html)
function:

``` r

nxy <- c(150, 100)
projgrid <- fm_evaluator(prmesh,
  xlim = range(PRborder[, 1]),
  ylim = range(PRborder[, 2]), dims = nxy
)
```

This lattice contains 150 × 100 locations. One can easily change the
resolution of the kriging prediction by changing `nxy`. Let us find the
cells that are outside the region of interest so that we do not plot the
estimates there.

``` r

xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[, 1], PRborder[, 2]))
```

Let us plot the locations that we will do prediction:

``` r

coord.prd <- projgrid$lattice$loc[xy.in, ]
plot(coord.prd, type = "p", cex = 0.1)
lines(PRborder)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")
```

![](rspde_inlabru_files/figure-html/plot_prd-1.png)

Let us now create a
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) of the
coordinates:

``` r

coord.prd.df <- data.frame(x1 = coord.prd[,1],
                            x2 = coord.prd[,2])
coord.prd.df <- st_as_sf(coord.prd.df, coords = c("x1", "x2"), 
                  crs = 4326)
```

Since we are using distance to the sea as a covariate, we also have to
calculate this covariate for the prediction locations. Finally, we add
the prediction location to our prediction
[`data.frame()`](https://rdrr.io/r/base/data.frame.html), namely,
`coord.prd.df`:

``` r

seaDist.prd <- apply(spDists(coord.prd,
  PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)
coord.prd.df$seaDist <- seaDist.prd
```

``` r

pred_obs <- predict(rspde_fit, coord.prd.df, 
        ~exp(Intercept + field + distSea))
```

Finally, we plot the results. First the predicted mean:

``` r

ggplot() + gg(pred_obs, geom = "tile",
    aes(fill = mean)) +
  geom_raster() +
  scale_fill_viridis()
```

![](rspde_inlabru_files/figure-html/unnamed-chunk-4-1.png)

Then, the std. deviations:

``` r

ggplot() + gg(pred_obs, geom = "tile",
    aes(fill = sd)) +
  geom_raster() +
  scale_fill_viridis()
```

![](rspde_inlabru_files/figure-html/plot_pred_sd_bru-1.png)

## An example with replicates

For this example we will simulate a data with replicates. We will use
the same example considered in the [Rational approximation with the
`rSPDE`
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette (the only difference is the way the data is organized). We also
refer the reader to this vignette for a description of the function
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md),
along with its methods (for instance, the
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method).

### Simulating the data

Let us consider a simple Gaussian linear model with 30 independent
replicates of a latent spatial field $`x(\mathbf{s})`$, observed at the
same $`m`$ locations, $`\{\mathbf{s}_1 , \ldots , \mathbf{s}_m \}`$, for
each replicate. For each $`i = 1,\ldots,m,`$ we have

``` math
\begin{align} 
y_i &= x_1(\mathbf{s}_i)+\varepsilon_i,\\
\vdots &= \vdots\\

y_{i+29m} &= x_{30}(\mathbf{s}_i) + \varepsilon_{i+29m},
\end{align}
```

where $`\varepsilon_1,\ldots,\varepsilon_{30m}`$ are iid normally
distributed with mean 0 and standard deviation 0.1.

We use the basis function representation of $`x(\cdot)`$ to define the
$`A`$ matrix linking the point locations to the mesh. We also need to
account for the fact that we have 30 replicates at the same locations.
To this end, the $`A`$ matrix we need can be generated by
[`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
function. The reason being that we are sampling $`x(\cdot)`$ directly
and not the latent vector described in the introduction of the [Rational
approximation with the `rSPDE`
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette.

We begin by creating the mesh:

``` r

m <- 200
loc_2d_mesh <- matrix(runif(m * 2), m, 2)
mesh_2d <- fm_mesh_2d(
  loc = loc_2d_mesh,
  cutoff = 0.05,
  offset = c(0.1, 0.4),
  max.edge = c(0.05, 0.5)
)
plot(mesh_2d, main = "")
points(loc_2d_mesh[, 1], loc_2d_mesh[, 2])
```

![](rspde_inlabru_files/figure-html/unnamed-chunk-5-1.png)

We then compute the $`A`$ matrix, which is needed for simulation, and
connects the observation locations to the mesh. To this end we will use
the
[`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
helper function, which is a wrapper that uses the functions
[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html),
[`fm_block()`](https://inlabru-org.github.io/fmesher/reference/fm_block.html)
and
[`fm_row_kron()`](https://inlabru-org.github.io/fmesher/reference/fm_row_kron.html)
from the `fmesher` package.

``` r

n.rep <- 30
A <- spde.make.A(
  mesh = mesh_2d,
  loc = loc_2d_mesh,
  index = rep(1:m, times = n.rep),
  repl = rep(1:n.rep, each = m)
)
```

Notice that for the simulated data, we should use the $`A`$ matrix from
[`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
function instead of the
[`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md).

We will now simulate a latent process with standard deviation
$`\sigma=1`$ and range $`0.1`$. We will use $`\nu=0.5`$ so that the
model has an exponential covariance function. To this end we create a
model object with the
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
function:

``` r

nu <- 0.5
sigma <- 1
range <- 0.1
kappa <- sqrt(8 * nu) / range
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi) * gamma(nu + 1)))
d <- 2
operator_information <- matern.operators(
  mesh = mesh_2d,
  nu = nu,
  range = range,
  sigma = sigma,
  m = 2,
  parameterization = "matern"
)
```

More details on this function can be found at the [Rational
approximation with the rSPDE
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette.

To simulate the latent process all we need to do is to use the
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method on the
`operator_information` object. We then obtain the simulated data $`y`$
by connecting with the $`A`$ matrix and adding the gaussian noise.

``` r

set.seed(1)
u <- simulate(operator_information, nsim = n.rep)
y <- as.vector(A %*% as.vector(u)) +
  rnorm(m * n.rep) * 0.1
```

The first replicate of the simulated random field as well as the
observation locations are shown in the following figure.

``` r

proj <- fm_evaluator(mesh_2d, dims = c(100, 100))

df_field <- data.frame(x = proj$lattice$loc[,1],
                        y = proj$lattice$loc[,2],
                        field = as.vector(fm_evaluate(proj, 
                        field = as.vector(u[, 1]))),
                        type = "field")

df_loc <- data.frame(x = loc_2d_mesh[, 1],
                      y = loc_2d_mesh[, 2],
                      field = y[1:m],
                      type = "locations")
df_plot <- rbind(df_field, df_loc)

ggplot(df_plot) + aes(x = x, y = y, fill = field) +
        facet_wrap(~type) + xlim(0,1) + ylim(0,1) + 
        geom_raster(data = df_field) +
        geom_point(data = df_loc, aes(colour = field),
        show.legend = FALSE) + 
        scale_fill_viridis() + scale_colour_viridis()
```

![](rspde_inlabru_files/figure-html/unnamed-chunk-9-1.png)

### Fitting the inlabru rSPDE model

Let us then use the rational SPDE approach to fit the data.

We begin by creating the model object.

``` r

rspde_model.rep <- rspde.matern(mesh = mesh_2d,
          parameterization = "spde") 
```

Let us now create the
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) and the vector
with the replicates indexes:

``` r

rep.df <- data.frame(y = y, x1 = rep(loc_2d_mesh[,1], n.rep),
                      x2 = rep(loc_2d_mesh[,2], n.rep))
rep.df <- st_as_sf(rep.df, coords = c("x1", "x2"))
repl <- rep(1:n.rep, each=m)
```

Let us create the component and fit. It is extremely important not to
forget the `replicate` when fitting model with the
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)
function. It will not produce warning and might fit some meaningless
model.

``` r

cmp.rep <-
  y ~ -1 + field(geometry,
    model = rspde_model.rep,
    replicate = repl
  )


rspde_fit.rep <-
  bru(cmp.rep,
    data = rep.df,
    family = "gaussian",
    options = list(num.threads = "1:1")
  )
```

We can get the summary:

``` r

summary(rspde_fit.rep)
```

    ## inlabru version: 2.14.1 
    ## INLA version: 26.05.02 
    ## Latent components:
    ## field: main = cgeneric(geometry), replicate = iid(repl)
    ## Observation models:
    ##   Model tag: <No tag>
    ##     Family: 'gaussian'
    ##     Data class: 'sf', 'data.frame'
    ##     Response class: 'numeric'
    ##     Predictor: y ~ field
    ##     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
    ##     Used components: effect[field], latent[] 
    ## Time used:
    ##     Pre = 0.174, Running = 86.6, Post = 2.92, Total = 89.7 
    ## Random effects:
    ##   Name     Model
    ##     field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                           mean    sd 0.025quant 0.5quant
    ## Precision for the Gaussian observations 94.718 4.699     86.209   94.470
    ## Theta1 for field                        -3.235 0.064     -3.386   -3.226
    ## Theta2 for field                         3.164 0.031      3.104    3.164
    ## Theta3 for field                        -0.633 0.029     -0.678   -0.636
    ##                                         0.975quant   mode
    ## Precision for the Gaussian observations    104.690 93.665
    ## Theta1 for field                            -3.145 -3.182
    ## Theta2 for field                             3.227  3.162
    ## Theta3 for field                            -0.566 -0.655
    ## 
    ## Marginal log-Likelihood:  -4498.14 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

and the summary in the user’s scale:

``` r

result_fit_rep <- rspde.result(rspde_fit.rep, "field", rspde_model.rep)
summary(result_fit_rep)
```

    ##            mean         sd 0.025quant   0.5quant 0.975quant       mode
    ## tau    0.039456 0.00246436  0.0338999  0.0398204  0.0430781  0.0417688
    ## kappa 23.679600 0.73348000 22.2988000 23.6568000 25.1784000 23.5982000
    ## nu     0.693669 0.01324030  0.6736920  0.6916500  0.7240520  0.6839210

``` r

result_df <- data.frame(
  parameter = c("tau", "kappa", "nu"),
  true = c(tau, kappa, nu),
  mean = c(
    result_fit_rep$summary.tau$mean,
    result_fit_rep$summary.kappa$mean,
    result_fit_rep$summary.nu$mean
  ),
  mode = c(
    result_fit_rep$summary.tau$mode,
    result_fit_rep$summary.kappa$mode,
    result_fit_rep$summary.nu$mode
  )
)
print(result_df)
```

    ##   parameter        true        mean        mode
    ## 1       tau  0.08920621  0.03945597  0.04176884
    ## 2     kappa 20.00000000 23.67961478 23.59818584
    ## 3        nu  0.50000000  0.69366870  0.68392136

Let us also obtain the summary on the `matern` parameterization:

``` r

result_fit_rep_matern <- rspde.result(rspde_fit.rep, "field", rspde_model.rep, 
                          parameterization = "matern")
summary(result_fit_rep_matern)
```

    ##              mean         sd 0.025quant  0.5quant 0.975quant     mode
    ## std.dev 1.0870200 0.01270990  1.0626200 1.0867500   1.112820 1.086590
    ## range   0.0990926 0.00356549  0.0924013 0.0990079   0.106238 0.098348
    ## nu      0.6936690 0.01324030  0.6736920 0.6916500   0.724052 0.683921

``` r

result_df_matern <- data.frame(
  parameter = c("std_dev", "range", "nu"),
  true = c(sigma, range, nu),
  mean = c(
    result_fit_rep_matern$summary.std.dev$mean,
    result_fit_rep_matern$summary.range$mean,
    result_fit_rep_matern$summary.nu$mean
  ),
  mode = c(
    result_fit_rep$summary.std.dev$mode,
    result_fit_rep$summary.range$mode,
    result_fit_rep$summary.nu$mode
  )
)
print(result_df_matern)
```

    ##   parameter true       mean      mode
    ## 1   std_dev  1.0 1.08701761 0.6839214
    ## 2     range  0.1 0.09909259 0.6839214
    ## 3        nu  0.5 0.69366870 0.6839214

## An example with a non-stationary model

Our goal now is to show how one can fit model with non-stationary
$`\sigma`$ (std. deviation) and non-stationary $`\rho`$ (a range
parameter). One can also use the parameterization in terms of
non-stationary SPDE parameters $`\kappa`$ and $`\tau`$.

For this example we will consider simulated data.

### Simulating the data

Let us consider a simple Gaussian linear model with a latent spatial
field $`x(\mathbf{s})`$, defined on the rectangle
$`(0,10) \times (0,5)`$, where the std. deviation and range parameter
satisfy the following log-linear regressions:
``` math
\begin{align}
\log(\sigma(\mathbf{s})) &= \theta_1 + \theta_3 b(\mathbf{s}),\\
\log(\rho(\mathbf{s})) &= \theta_2 + \theta_3 b(\mathbf{s}),
\end{align}
```
where $`b(\mathbf{s}) = (s_1-5)/10`$. We assume the data is observed at
$`m`$ locations, $`\{\mathbf{s}_1 , \ldots , \mathbf{s}_m \}`$. For each
$`i = 1,\ldots,m,`$ we have

``` math
y_i = x_1(\mathbf{s}_i)+\varepsilon_i,
```

where $`\varepsilon_1,\ldots,\varepsilon_{m}`$ are iid normally
distributed with mean 0 and standard deviation 0.1.

We begin by defining the domain and creating the mesh:

``` r

rec_domain <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)

mesh <- fm_mesh_2d(loc.domain = rec_domain, cutoff = 0.1, 
  max.edge = c(0.5, 1.5), offset = c(0.5, 1.5))
```

We follow the same structure as `INLA`. However, `INLA` only allows one
to specify `B.tau` and `B.kappa` matrices, and, in `INLA`, if one wants
to parameterize in terms of range and standard deviation one needs to do
it manually. Here we provide the option to directly provide the matrices
`B.sigma` and `B.range`.

The usage of the matrices `B.tau` and `B.kappa` are identical to the
corresponding ones in
[`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html)
function. The matrices `B.sigma` and `B.range` work in the same way, but
they parameterize the stardard deviation and range, respectively.

The columns of the `B` matrices correspond to the same parameter. The
first column does not have any parameter to be estimated, it is a
constant column.

So, for instance, if one wants to share a parameter with both `sigma`
and `range` (or with both `tau` and `kappa`), one simply let the
corresponding column to be nonzero on both `B.sigma` and `B.range` (or
on `B.tau` and `B.kappa`).

We will assume $`\nu = 0.8`$, $`\theta_1 = 0, \theta_2 = 1`$ and
$`\theta_3=1`$. Let us now build the model to obtain the sample with the
[`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md)
function:

``` r

nu <- 0.8
true_theta <- c(0,1, 1)
B.sigma = cbind(0, 1, 0, (mesh$loc[,1] - 5) / 10)
B.range = cbind(0, 0, 1, (mesh$loc[,1] - 5) / 10)

# SPDE model
op_cov_ns <- spde.matern.operators(mesh = mesh, 
  theta = true_theta,
  nu = nu,
  B.sigma = B.sigma, 
  B.range = B.range, m = 2,
  parameterization = "matern")
```

Let us now sample the data with the
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method:

``` r

u <- as.vector(simulate(op_cov_ns, seed = 123))
```

Let us now obtain 600 random locations on the rectangle and compute the
$`A`$ matrix:

``` r

m <- 600
loc_mesh <- cbind(runif(m) * 10, runif(m) * 5)

A <- spde.make.A(
  mesh = mesh,
  loc = loc_mesh
)
```

We can now generate the response vector `y`:

``` r

y <- as.vector(A %*% as.vector(u)) + rnorm(m) * 0.1
```

### Fitting the inlabru rSPDE model

Let us then use the rational SPDE approach to fit the data.

We begin by creating the model object. We are creating a new one so that
we do not start the estimation at the true values.

``` r

rspde_model_nonstat <- rspde.matern(mesh = mesh,
  B.sigma = B.sigma,
  B.range = B.range,
  parameterization = "matern") 
```

Let us now create the
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) and the vector
with the replicates indexes:

``` r

nonstat_df <- data.frame(y = y, x1 = loc_mesh[,1],
                      x2 = loc_mesh[,2])
nonstat_df <- st_as_sf(nonstat_df, coords = c("x1", "x2"))
```

Let us create the component and fit. It is extremely important not to
forget the `replicate` when fitting model with the
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)
function. It will not produce warning and might fit some meaningless
model.

``` r

cmp_nonstat <-
  y ~ -1 + field(geometry,
    model = rspde_model_nonstat
  )


rspde_fit_nonstat <-
  bru(cmp_nonstat,
    data = nonstat_df,
    family = "gaussian",
    options = list(verbose = FALSE,
                   num.threads = "1:1")
  )
```

We can get the summary:

``` r

summary(rspde_fit_nonstat)
```

    ## inlabru version: 2.14.1 
    ## INLA version: 26.05.02 
    ## Latent components:
    ## field: main = cgeneric(geometry)
    ## Observation models:
    ##   Model tag: <No tag>
    ##     Family: 'gaussian'
    ##     Data class: 'sf', 'data.frame'
    ##     Response class: 'numeric'
    ##     Predictor: y ~ field
    ##     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
    ##     Used components: effect[field], latent[] 
    ## Time used:
    ##     Pre = 0.155, Running = 20.3, Post = 0.134, Total = 20.6 
    ## Random effects:
    ##   Name     Model
    ##     field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                            mean    sd 0.025quant 0.5quant
    ## Precision for the Gaussian observations 104.575 7.624     90.683  104.193
    ## Theta1 for field                         -0.072 0.079     -0.223   -0.074
    ## Theta2 for field                          0.808 0.098      0.616    0.808
    ## Theta3 for field                          1.200 0.158      0.911    1.193
    ## Theta4 for field                          0.019 0.059     -0.088    0.016
    ##                                         0.975quant    mode
    ## Precision for the Gaussian observations    120.690 103.204
    ## Theta1 for field                             0.089  -0.082
    ## Theta2 for field                             1.003   0.806
    ## Theta3 for field                             1.531   1.161
    ## Theta4 for field                             0.144   0.002
    ## 
    ## Marginal log-Likelihood:  2.06 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

We can obtain outputs with respect to parameters in the original scale
by using the function
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.md):

``` r

result_fit_nonstat <- rspde.result(rspde_fit_nonstat, "field", rspde_model_nonstat)
summary(result_fit_nonstat)
```

    ##                     mean        sd 0.025quant   0.5quant 0.975quant       mode
    ## Theta1.matern -0.0724459 0.0792918  -0.223314 -0.0741678  0.0887946 -0.0817976
    ## Theta2.matern  0.8082280 0.0981722   0.616486  0.8077060  1.0030300  0.8055140
    ## Theta3.matern  1.2000400 0.1578990   0.911145  1.1932600  1.5309400  1.1614500
    ## nu             1.0092600 0.0294487   0.956224  1.0073600  1.0712600  1.0009800

Let us compare the mean to the true values of the parameters:

``` r

summ_res_nonstat <- summary(result_fit_nonstat)
result_df <- data.frame(
  parameter = result_fit_nonstat$params,
  true = c(true_theta, nu),
  mean = summ_res_nonstat[,1],
  mode = summ_res_nonstat[,6]
)
print(result_df)
```

    ##       parameter true       mean       mode
    ## 1 Theta1.matern  0.0 -0.0724459 -0.0817976
    ## 2 Theta2.matern  1.0  0.8082280  0.8055140
    ## 3 Theta3.matern  1.0  1.2000400  1.1614500
    ## 4            nu  0.8  1.0092600  1.0009800

We can also plot the posterior densities. To this end we will use the
[`gg_df()`](https://davidbolin.github.io/rSPDE/reference/gg_df.md)
function, which creates `ggplot2` user-friendly data frames:

``` r

posterior_df_fit <- gg_df(result_fit_nonstat)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](rspde_inlabru_files/figure-html/plot_post_nonstat-1.png)

## Comparing the results by cross-validation

We can compare the models fitted by `inlabru` by using the function
[`cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.md).
To illustrate, we will consider the nonstationary model
`rspde_fit_nonstat` fitted in the previous example and a stationary fit
of the same dataset.

Let us, then, fit a stationary model with the previous dataset. We start
by defining the stationary model:

``` r

rspde_model_stat <- rspde.matern(mesh = mesh)
```

Then, `inlabru`’s component:

``` r

cmp_stat <-
  y ~ -1 + field(geometry,
    model = rspde_model_stat
  )
```

We can now fit the model:

``` r

rspde_fit_stat <-
  bru(cmp_stat,
    data = nonstat_df,
    family = "gaussian",
    options = list(verbose = FALSE,
                   num.threads = "1:1")
  )
```

To perform cross-validation, we create a list with the fitted models,
and we pass this list to the
[`cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.md)
function. It is also important to create a named list, so that the
output has meaningful names for the models. We will perform a
`leave percentage out` cross-validation, with the default that fits the
model on 20% of the data, to predict 80% of the data.

Let us create the models list:

``` r

models <- list(stationary = rspde_fit_stat, 
                nonstationary = rspde_fit_nonstat)
```

We will now run the cross-validation on the models above. We set the
`cv_type` to `lpo` to perform the leave percentage out cross-validation,
there are also the `k-fold` (default) and `loo` options to perform
k-fold and leave one out cross-validations, respectively. Observe that
by default we are performing a pseudo cross-validation, that is, we will
not refit the model for each fold, however only the training data will
be used to perform the prediction.

``` r

cv_result <- cross_validation(models, cv_type = "lpo", print = FALSE)
```

We can now look at the results by printing `cv_result`. Observe that the
best model with respect to each score is displayed in the last row.

``` r

cv_result
```

    ##           Model               mse               mae               dss
    ## 1    stationary 0.135844439112595 0.270711197652152 -1.20847740273866
    ## 2 nonstationary 0.132761955224834  0.26883965900393 -1.17628089962464
    ##            Best     nonstationary     nonstationary        stationary
    ##                crps             scrps
    ## 1 0.190113166035239 0.491644123284751
    ## 2 0.188069159815355 0.485678685130971
    ##       nonstationary     nonstationary

The
[`cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.md)
function also has the following useful options:

- `return_score_folds` option, so that the scores for each fold can be
  returned in order to create confidence regions for the scores.
- `return_train_test` To return the train and test indexes that were
  used to perform the cross-validation.
- `true_CV` To perform true cross-validation, that is, the data will be
  fit again for each fold, which is more costly.
- `train_test_indexes` In which the user can provide the indexes for the
  train and test sets.

More details can be found in the manual page of the
[`cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.md)
function.

## Further options of the `inlabru` implementation

There are several additional options that are available. For instance,
it is possible to change the order of the rational approximation, the
upper bound for the smoothness parameter (which may speed up the fit),
change the priors, change the type of the rational approximation, among
others. These options are described in the “Further options of the
`rSPDE`-`INLA` implementation” section of the [R-INLA implementation of
the rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md)
vignette. Observe that all these options are passed to the model through
the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function, and therefore the resulting model object can directly be used
in the
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.html)
function, in an identical manner to the examples above.

## References

Bolin, David, Alexandre B. Simas, and Zhen Xiong. 2023.
“Covariance-Based Rational Approximations of Fractional SPDEs for
Computationally Efficient Bayesian Inference.” *Journal of Computational
and Graphical Statistics*, ahead of print.
<https://doi.org/10.1080/10618600.2022.2139648>.

Lindgren, Finn, Håvard Rue, and Johan Lindström. 2011. “An Explicit Link
Between Gaussian Fields and Gaussian Markov Random Fields: The
Stochastic Partial Differential Equation Approach.” *Journal of the
Royal Statistical Society. Series B. Statistical Methodology* 73 (4):
423–98.

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

    ## Linking to GEOS 3.14.1, GDAL 3.12.2, PROJ 9.7.1; sf_use_s2() is TRUE

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

    ## inlabru version: 2.15.0 
    ## INLA version: 26.08.22 
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
    ##     Pre = 0.144, Running = 4.34, Post = 0.0951, Total = 4.58 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.041       1.86    1.941      2.022 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     distSea RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                    mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.442    1.039     12.496
    ## Precision for distSea                          7399.163 4016.397   2209.293
    ## Theta1 for field                                 -0.968    0.400     -1.800
    ## Theta2 for field                                  1.400    0.274      0.858
    ## Theta3 for field                                 -0.387    0.384     -1.108
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations   14.408     16.585   14.346
    ## Precision for distSea                          6546.062  17525.489 5021.593
    ## Theta1 for field                                 -0.953     -0.227   -0.884
    ## Theta2 for field                                  1.401      1.935    1.407
    ## Theta3 for field                                 -0.399      0.403   -0.450
    ## 
    ## Marginal log-Likelihood:  -1256.56 
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
\nu = 2\Big(\frac{\exp(\theta_3)}{1+\exp(\theta_3)}\Big).
```
The number 2 is the default upper bound for $`\nu`$, discussed in
[R-INLA implementation of the rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.md).

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
summary(result_fit)
```

    ##           mean       sd 0.025quant 0.5quant 0.975quant     mode
    ## tau   0.410534 0.161814   0.166759 0.387484   0.792517 0.337161
    ## kappa 4.208640 1.159930   2.369300 4.062260   6.896100 3.786620
    ## nu    0.814581 0.179577   0.498529 0.801650   1.194810 0.763370

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
summary(result_fit_matern)
```

    ##             mean        sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 0.315845 0.0388049   0.244795 0.313921   0.396473 0.309260
    ## range   0.640005 0.2050680   0.337204 0.606712   1.141910 0.542825
    ## nu      0.814581 0.1795770   0.498529 0.801650   1.194810 0.763370

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

    ## inlabru version: 2.15.0 
    ## INLA version: 26.08.22 
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
    ##     Pre = 0.163, Running = 42.4, Post = 2.95, Total = 45.5 
    ## Random effects:
    ##   Name     Model
    ##     field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                           mean    sd 0.025quant 0.5quant
    ## Precision for the Gaussian observations 91.187 4.584     82.468   91.080
    ## Theta1 for field                        -3.259 0.047     -3.362   -3.256
    ## Theta2 for field                         3.094 0.032      3.033    3.094
    ## Theta3 for field                        -0.578 0.027     -0.626   -0.579
    ##                                         0.975quant   mode
    ## Precision for the Gaussian observations     100.51 90.888
    ## Theta1 for field                             -3.18 -3.239
    ## Theta2 for field                              3.16  3.094
    ## Theta3 for field                             -0.52 -0.586
    ## 
    ## Marginal log-Likelihood:  -4533.21 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

and the summary in the user’s scale:

``` r

result_fit_rep <- rspde.result(rspde_fit.rep, "field", rspde_model.rep)
summary(result_fit_rep)
```

    ##             mean         sd 0.025quant   0.5quant 0.975quant       mode
    ## tau    0.0384701 0.00178786  0.0347072  0.0386004  0.0416442  0.0391484
    ## kappa 22.0865000 0.69626000 20.7622000 22.0706000 23.4963000 22.0328000
    ## nu     0.7190350 0.01240090  0.6970970  0.7180950  0.7454970  0.7149930

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
    ## 1       tau  0.08920621  0.03847007  0.03914841
    ## 2     kappa 20.00000000 22.08647014 22.03280766
    ## 3        nu  0.50000000  0.71903547  0.71499270

Let us also obtain the summary on the `matern` parameterization:

``` r

result_fit_rep_matern <- rspde.result(rspde_fit.rep, "field", rspde_model.rep, 
                          parameterization = "matern")
summary(result_fit_rep_matern)
```

    ##             mean         sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 1.099680 0.01195240   1.076660 1.099590   1.123630 1.099370
    ## range   0.108483 0.00381006   0.101274 0.108376   0.116113 0.107420
    ## nu      0.719035 0.01240090   0.697097 0.718095   0.745497 0.714993

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
    result_fit_rep_matern$summary.std.dev$mode,
    result_fit_rep_matern$summary.range$mode,
    result_fit_rep_matern$summary.nu$mode
  )
)
print(result_df_matern)
```

    ##   parameter true      mean      mode
    ## 1   std_dev  1.0 1.0996757 1.0993732
    ## 2     range  0.1 0.1084833 0.1074197
    ## 3        nu  0.5 0.7190355 0.7149927

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

    ## inlabru version: 2.15.0 
    ## INLA version: 26.08.22 
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
    ##     Pre = 0.121, Running = 17, Post = 0.106, Total = 17.3 
    ## Random effects:
    ##   Name     Model
    ##     field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                            mean    sd 0.025quant 0.5quant
    ## Precision for the Gaussian observations 105.826 8.511     89.341  105.729
    ## Theta1 for field                         -0.065 0.094     -0.234   -0.070
    ## Theta2 for field                          0.809 0.099      0.628    0.805
    ## Theta3 for field                          1.224 0.259      0.766    1.209
    ## Theta4 for field                          0.030 0.071     -0.096    0.026
    ##                                         0.975quant    mode
    ## Precision for the Gaussian observations    122.789 106.128
    ## Theta1 for field                             0.136  -0.095
    ## Theta2 for field                             1.016   0.785
    ## Theta3 for field                             1.780   1.133
    ## Theta4 for field                             0.182   0.006
    ## 
    ## Marginal log-Likelihood:  2.79 
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

    ##                    mean        sd 0.025quant   0.5quant 0.975quant       mode
    ## Theta1.matern -0.064955 0.0944432  -0.233872 -0.0700078   0.136008 -0.0945858
    ## Theta2.matern  0.809048 0.0987540   0.628466  0.8048210   1.016220  0.7846640
    ## Theta3.matern  1.224040 0.2592870   0.766125  1.2086900   1.780060  1.1332500
    ## nu             1.014910 0.0351907   0.952523  1.0122900   1.089720  1.0033200

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

    ##       parameter true      mean       mode
    ## 1 Theta1.matern  0.0 -0.064955 -0.0945858
    ## 2 Theta2.matern  1.0  0.809048  0.7846640
    ## 3 Theta3.matern  1.0  1.224040  1.1332500
    ## 4            nu  0.8  1.014910  1.0033200

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

cv_result <- rSPDE::cross_validation(models, cv_type = "lpo",
                                    true_CV = FALSE, print = FALSE)
```

    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'
    ## Warning in inla.mode.section(file = file.ini, cont.mode, data.dir): ***
    ## control.mode *** option 'restart=TRUE' is set to 'FALSE' since 'fixed=TRUE'

We can now look at the results by printing `cv_result`. Observe that the
best model with respect to each score is displayed in the last row.

``` r

cv_result
```

    ##           Model               mse               mae               dss
    ## 1    stationary 0.133341184740196 0.270557316448001 -1.31508426559492
    ## 2 nonstationary  0.13156561601604 0.269249923003719 -1.33823685372121
    ##            Best     nonstationary     nonstationary     nonstationary
    ##                crps             scrps
    ## 1 0.191833393857515  0.48139511427403
    ## 2 0.189941415871072 0.476060866689055
    ##       nonstationary     nonstationary

The
[`cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.md)
function also has the following useful options:

- `return_scores_folds` option, so that the scores for each fold can be
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

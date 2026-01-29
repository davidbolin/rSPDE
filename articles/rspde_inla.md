# R-INLA implementation of the rational SPDE approach

## Introduction

In this vignette we will present the [`R-INLA`](https://www.r-inla.org)
implementation of the rational SPDE approach. For theoretical details we
refer the reader to the [Rational approximation with the `rSPDE`
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette and to [Bolin, Simas, and Xiong
(2023)](https://doi.org/10.1080/10618600.2023.2231051).

We begin by providing a step-by-step illustration on how to use our
implementation. To this end we will consider a real world data set that
consists of precipitation measurements from the Paraná region in Brazil.

After the initial model fitting, we will show how to change some
parameters of the model. In the end, we will also provide an example in
which we have replicates.

It is important to mention that one can improve the performance by using
the PARDISO solver. Please, go to
<https://www.pardiso-project.org/r-inla/#license> to apply for a
license. Also, use
[`inla.pardiso()`](https://rdrr.io/pkg/INLA/man/pardiso.html) for
instructions on how to enable the PARDISO sparse library.

## Example with real data

To illustrate our implementation of `rSPDE` in
[`R-INLA`](https://www.r-inla.org) we will consider a dataset available
in [`R-INLA`](https://www.r-inla.org). This data has also been used to
illustrate the SPDE approach, see for instance the book [Advanced
Spatial Modeling with Stochastic Partial Differential Equations Using R
and
INLA](https://www.routledge.com/Advanced-Spatial-Modeling-with-Stochastic-Partial-Differential-Equations/Krainski-Gomez-Rubio-Bakka-Lenzi-Castro-Camilo-Simpson-Lindgren-Rue/p/book/9780367570644)
and also the vignette [Spatial Statistics using R-INLA and Gaussian
Markov random
fields](https://sites.stat.washington.edu/peter/591/INLA.html). See also
[Lindgren, Rue, and Lindström
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
$$\Gamma(\mu,\phi):\pi(y) = \frac{1}{\Gamma(\phi)}\left( \frac{\phi}{\mu} \right)^{\phi}y^{\phi - 1}\exp\left( - \frac{\phi y}{\mu} \right).$$
In this parameterization, the distribution has expected value
$E(x) = \mu$ and variance $V(x) = \mu^{2}/(\phi)$, where$1/\phi$ is a
dispersion parameter.

In this example $\mu$ will be modeled using a stochastic model that
includes both covariates and spatial structure, resulting in the latent
Gaussian model for the precipitation measurements $$\begin{aligned}
{y_{i} \mid \mu\left( s_{i} \right),\theta} & {\sim \Gamma\left( \mu\left( s_{i} \right),c\phi \right)} \\
{\log\left( \mu(s) \right)} & {= \eta(s) = \sum\limits_{k}f_{k}\left( c_{k}(s) \right) + u(s)} \\
\theta & {\sim \pi(\theta)}
\end{aligned},$$

where $y_{i}$ denotes the measurement taken at location $s_{i}$,
$c_{k}(s)$ are covariates, $u(s)$ is a mean-zero Gaussian Matérn field,
and $\theta$ is a vector containing all parameters of the model,
including smoothness of the field. That is, by using the `rSPDE` model
we will also be able to estimate the smoothness of the latent field.

### Examining the data

We will be using [`R-INLA`](https://www.r-inla.org). To install
[`R-INLA`](https://www.r-inla.org) go to [R-INLA
Project](https://www.r-inla.org/download-install).

We begin by loading some libraries we need to get the data and build the
plots.

``` r
library(ggplot2)
library(INLA)
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

Let us build plot the precipitation observations using `ggplot`:

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

![](rspde_inla_files/figure-html/plot_precipitations-1.png)

The red line in the figure shows the coast line, and we expect the
distance to the coast to be a good covariate for precipitation. This
covariate is not available, so let us calculate it for each observation
location:

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

![](rspde_inla_files/figure-html/plot_prec_as_func-1.png)

``` r
par(mfrow = c(1, 1))
```

### Creating the rSPDE model

To use the [`R-INLA`](https://www.r-inla.org) implementation of the
`rSPDE` model we need to load the functions:

``` r
library(rSPDE)
```

The `rSPDE`-`INLA` implementation is very reminiscent of
[`R-INLA`](https://www.r-inla.org), so its usage should be
straightforward for [`R-INLA`](https://www.r-inla.org) users. For
instance, to create a `rSPDE` model, one would use
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
in place of
[`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html).
To create an index, one should use
[`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
in place of
[`inla.spde.make.index()`](https://rdrr.io/pkg/INLA/man/inla.spde.make.index.html).
To create the `A` matrix, one should use
[`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md)
in place of
[`inla.spde.make.A()`](https://rdrr.io/pkg/INLA/man/inla.spde.make.A.html),
and so on.

The main differences when comparing the arguments between the
`rSPDE`-`INLA` implementation and the standard SPDE implementation in
[`R-INLA`](https://www.r-inla.org), are the `nu` and `rspde.order`
arguments, which are present in `rSPDE`-`INLA` implementation. We will
see below how use these arguments.

#### Mesh

We can use `fmesher` for creating the mesh. We begin by loading the
`fmesher` package:

``` r
library(fmesher)
```

Let us create a mesh which is based on a non-convex hull to avoid adding
many small triangles outside the domain of interest:

``` r
prdomain <- fm_nonconvex_hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)

plot(prmesh, asp = 1, main = "")
lines(PRborder, col = 3)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")
```

![](rspde_inla_files/figure-html/mesh_creation-1.png)

#### The observation matrix

We now create the $A$ matrix, that connects the mesh to the observation
locations and then create the `rSPDE` model.

For this task, as we mentioned earlier, we need to use an
`rSPDE`specific function, whose name is very reminiscent to
[`R-INLA`](https://www.r-inla.org)’s standard SPDE approach, namely
[`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md)
(in place of [`R-INLA`](https://www.r-inla.org)’s
[`inla.spde.make.A()`](https://rdrr.io/pkg/INLA/man/inla.spde.make.A.html)).
The reason for the need of this specific function is that the size of
the $A$ matrix depends on the order of the rational approximation. The
details can be found in the introduction of the [Rational approximation
with the `rSPDE`
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette.

The default order is 2 for our covariance-based rational approximation.
As mentioned in the introduction of the [Rational approximation with the
`rSPDE`
package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
vignette, an approximation of order 2 in the covariance-based rational
approximation has approximately the same computational cost as the
operator-based rational approximation of order 1.

Recall that the latent process $u$ is a solution of
$$\left( \kappa^{2}I - \Delta \right)^{\alpha/2}(\tau u) = \mathcal{W},$$
where $\alpha = \nu + d/2$. We want to estimate all three parameters
$\tau,\kappa$ and $\nu$, which is the default option of  
the `rSPDE`-`INLA` implementation. However, there is also an option to
fix the smoothness parameter $\nu$ to some predefined value and only
estimate $\tau$ and $\kappa$. This will be discussed later.

In this first example we will assume we want a rational approximation of
order 1. To this end we can use the
[`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md)
function. Since we will assume order 1 and that we want to estimate
smoothness, which are the default options in this function, the required
parameters are simply the mesh and the locations:

``` r
Abar <- rspde.make.A(mesh = prmesh, loc = coords)
```

#### Setting up the rSPDE model

To set up an `rSPDE`model, all we need is the mesh. By default it will
assume that we want to estimate the smoothness parameter $\nu$ and to do
a covariance-based rational approximation of order 2.

Later in this vignette we will also see other options for setting up
`rSPDE` models such as keeping the smoothness parameter fixed and/or
increasing the order of the covariance-based rational approximation.

Therefore, to set up a model all we have to do is use the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function:

``` r
rspde_model <- rspde.matern(mesh = prmesh)
```

Note that this function is very reminiscent of
[`R-INLA`](https://www.r-inla.org)’s
[`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html)
function. This is a pattern we have tried to keep consistent in the
package: All the `rSPDE` versions of some
[`R-INLA`](https://www.r-inla.org) function will either replace `inla`
or `inla.spde` or `inla.spde2` by `rspde`.

#### The `inla.stack`

Since the covariates are already evaluated at the observation locations,
we only want to apply the $A$ matrix to the spatial effect and not the
fixed effects. We can use the
[`inla.stack()`](https://rdrr.io/pkg/INLA/man/inla.stack.html) function.

The difference, however, is that we need to use the function
[`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
(in place of the standard
[`inla.spde.make.index()`](https://rdrr.io/pkg/INLA/man/inla.spde.make.index.html))
to create the index.

If one is using the default options, that is, to estimate the smoothness
parameter $\nu$ and to do a rational approximation of order 2, the usage
of
[`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
is identical to the usage of
[`inla.spde.make.index()`](https://rdrr.io/pkg/INLA/man/inla.spde.make.index.html):

``` r
mesh.index <- rspde.make.index(name = "field", mesh = prmesh)
```

We can then create the stack in a standard manner:

``` r
stk.dat <- inla.stack(
  data = list(y = Y), A = list(Abar, 1), tag = "est",
  effects = list(
    c(
      mesh.index
    ),
    list(
      seaDist = inla.group(seaDist),
      Intercept = 1
    )
  )
)
```

Here the observation matrix $A$ is applied to the spatial effect and the
intercept while an identity observation matrix, denoted by $1$, is
applied to the covariates. This means the covariates are unaffected by
the observation matrix.

The observation matrices in $A = list(Abar,1)$ are used to link the
corresponding elements in the effects-list to the observations. Thus in
our model the latent spatial field `mesh.index` and the intercept are
linked to the log-expectation of the observations, i.e. $\eta(s)$,
through the $A$-matrix. The covariates, on the other hand, are linked
directly to $\eta(s)$. The `stk.dat` object defined above implies the
following principal linkage between model components and observations
$$\eta(s) \sim Ax(s) + A{\mspace{6mu}\text{Intercept}} + \text{seaDist}.$$$\eta(s)$
will then be used in the observation-likelihood,
$$y_{i} \mid \eta\left( s_{i} \right),\theta \sim \Gamma\left( \exp\left( \eta\left( s_{i} \right) \right),c\phi \right).$$

### Model fitting

We will build a model using the distance to the sea $x_{i}$ as a
covariate through an improper CAR(1) model with
$\beta_{ij} = 1(i \sim j)$, which [`R-INLA`](https://www.r-inla.org)
calls a random walk of order 1.

``` r
f.s <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model)
```

Here `-1` is added to remove R’s implicit intercept, which is replaced
by the explicit `+Intercept` from when we created the stack.

To fit the model we proceed as in the standard SPDE approach and we
simply call [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

``` r
rspde_fit <- inla(f.s,
  family = "Gamma", data = inla.stack.data(stk.dat),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
```

### INLA results

We can look at some summaries of the posterior distributions for the
parameters, for example the fixed effects (i.e. the intercept) and the
hyper-parameters (i.e. dispersion in the gamma likelihood, the precision
of the RW1, and the parameters of the spatial field):

``` r
summary(rspde_fit)
```

    ## Time used:
    ##     Pre = 0.384, Running = 5.76, Post = 0.0342, Total = 6.18 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.942 0.042       1.86    1.942      2.023 1.942   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.43    1.039     12.484
    ## Precision for seaDist                          7626.99 4390.004   2376.846
    ## Theta1 for field                                 -3.99    3.060    -11.165
    ## Theta2 for field                                  1.95    0.635      0.983
    ## Theta3 for field                                  2.23    2.648     -1.512
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations    14.40   1.66e+01   14.342
    ## Precision for seaDist                           6593.47   1.90e+04 4965.490
    ## Theta1 for field                                  -3.54   3.31e-01   -1.142
    ## Theta2 for field                                   1.87   3.41e+00    1.467
    ## Theta3 for field                                   1.85   8.44e+00   -0.186
    ## 
    ## Marginal log-Likelihood:  -1254.87 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

Let $\theta_{1} = \text{Theta1}$, $\theta_{2} = \text{Theta2}$ and
$\theta_{3} = \text{Theta3}$. In terms of the SPDE
$$\left( \kappa^{2}I - \Delta \right)^{\alpha/2}(\tau u) = \mathcal{W},$$
where $\alpha = \nu + d/2$, we have that
$$\tau = \exp\left( \theta_{1} \right),\quad\kappa = \exp\left( \theta_{2} \right),$$
and by default
$$\nu = 2(\frac{\exp\left( \theta_{3} \right)}{1 + \exp\left( \theta_{3} \right)}).$$
The number 2 comes from the upper bound for $\nu$, which will be
discussed later in this vignette.

In general, we have
$$\nu = \nu_{UB}(\frac{\exp\left( \theta_{3} \right)}{1 + \exp\left( \theta_{3} \right)}),$$
where $\nu_{UB}$ is the value of the upper bound for the smoothness
parameter $\nu$.

Another choice for prior for $\nu$ is a truncated lognormal distribution
and will also be discussed later in this vignette.

### `rSPDE`-`INLA` results

We can obtain outputs with respect to parameters in the original scale
by using the function
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.md):

``` r
result_fit <- rspde.result(rspde_fit, "field", rspde_model)
```

    ## Warning in rspde.result(rspde_fit, "field", rspde_model): the mean or mode of
    ## nu is very close to nu.upper.bound, please consider increasing nu.upper.bound,
    ## and refitting the model.

``` r
summary(result_fit)
```

    ##           mean       sd  0.025quant  0.5quant 0.975quant        mode
    ## tau   0.205752 0.438083 1.54764e-05 0.0321271    1.42160 4.10688e-08
    ## kappa 8.795130 7.643860 2.69096e+00 6.3483300   29.68280 3.85637e+00
    ## nu    1.488500 0.530635 3.56472e-01 1.7063900    1.99952 1.99999e+00

As mentioned above, when we create the model object with
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md),
we must choose an upper bound for `nu` by using the argument
`nu.upper.bound`. If such an argument is not passed, the default value
of `2` will be used. However, if the mean or mode of `nu` is too close
to `nu.upper.bound`, a warning will be given suggesting the user to
increase `nu.upper.bound` and refit the data.

To create plots of the posterior marginal densities, we can use the
[`gg_df()`](https://davidbolin.github.io/rSPDE/reference/gg_df.md)
function, which creates `ggplot2`-friendly data frames. The following
figure shows the posterior marginal densities of the three parameters
using the
[`gg_df()`](https://davidbolin.github.io/rSPDE/reference/gg_df.md)
function.

``` r
posterior_df_fit <- gg_df(result_fit)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](rspde_inla_files/figure-html/plot_post-1.png)

This function is reminiscent to the
[`inla.spde.result()`](https://rdrr.io/pkg/INLA/man/inla.spde.result.html)
function with the main difference that it has the
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
implemented.

We can also obtain the results for the `matern` parameterization by
setting the `parameterization` argument to `matern`:

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

    ##             mean       sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 1.524550 9.133690 -0.0314183 0.505780   9.237280 0.320613
    ## range   0.459179 0.234540  0.0750641 0.445388   0.975396 0.423208
    ## nu      1.488500 0.530635  0.3564720 1.706390   1.999520 1.999990

In a similar manner, we can obtain posterior plots on the `matern`
parameterization:

``` r
posterior_df_fit_matern <- gg_df(result_fit_matern)

ggplot(posterior_df_fit_matern) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](rspde_inla_files/figure-html/plot_post_matern-1.png)

### Predictions

Let us now obtain predictions (i.e. do kriging) of the expected
precipitation on a dense grid in the region.

We begin by creating the grid in which we want to do the predictions. To
this end, we can use the
[`rspde.mesh.projector()`](https://davidbolin.github.io/rSPDE/reference/rspde.mesh.project.md)
function. This function has the same arguments as the function
[`inla.mesh.projector()`](https://rdrr.io/pkg/INLA/man/inla.mesh.project.html),
with the only difference being that the `rSPDE` version also has an
argument `nu` and an argument `rspde.order`. Thus, we proceed in the
same fashion as we would in [`R-INLA`](https://www.r-inla.org)’s
standard SPDE implementation:

``` r
nxy <- c(150, 100)
projgrid <- rspde.mesh.projector(prmesh,
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

![](rspde_inla_files/figure-html/plot_prd-1.png)

Now, there are a few ways we could calculate the kriging prediction. The
simplest way is to evaluate the mean of all individual random effects in
the linear predictor and then to calculate the exponential of their sum
(since $\mu(s) = \exp\left( \eta(s) \right)$ ). A more accurate way is
to calculate the prediction jointly with the estimation, which
unfortunately is quite computationally expensive if we do prediction on
a fine grid. However, in this illustration, we proceed with this option
to show how one can do it.

To this end, first, link the prediction coordinates to the mesh nodes
through an $A$ matrix

``` r
A.prd <- projgrid$proj$A[xy.in, ]
```

Since we are using distance to the sea as a covariate, we also have to
calculate this covariate for the prediction locations.

``` r
seaDist.prd <- apply(spDists(coord.prd,
  PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)
```

We now make a stack for the prediction locations. We have no data at the
prediction locations, so we set `y= NA`. We then join this stack with
the estimation stack.

``` r
ef.prd <- list(
  c(mesh.index),
  list(
    long = inla.group(coord.prd[
      ,
      1
    ]), lat = inla.group(coord.prd[, 2]),
    seaDist = inla.group(seaDist.prd),
    Intercept = 1
  )
)
stk.prd <- inla.stack(
  data = list(y = NA),
  A = list(A.prd, 1), tag = "prd",
  effects = ef.prd
)
stk.all <- inla.stack(stk.dat, stk.prd)
```

Doing the joint estimation takes a while, and we therefore turn off the
computation of certain things that we are not interested in, such as the
marginals for the random effect. We will also use a simplified
integration strategy (actually only using the posterior mode of the
hyper-parameters) through the command
`control.inla = list(int.strategy = "eb")`, i.e. empirical Bayes.

``` r
rspde_fitprd <- inla(f.s,
  family = "Gamma",
  data = inla.stack.data(stk.all),
  control.predictor = list(
    A = inla.stack.A(stk.all),
    compute = TRUE, link = 1
  ),
  control.compute = list(
    return.marginals = FALSE,
    return.marginals.predictor = FALSE
  ),
  control.inla = list(int.strategy = "eb"),
  num.threads = "1:1"
)
```

We then extract the indices to the prediction nodes and then extract the
mean and the standard deviation of the response:

``` r
id.prd <- inla.stack.index(stk.all, "prd")$data
m.prd <- rspde_fitprd$summary.fitted.values$mean[id.prd]
sd.prd <- rspde_fitprd$summary.fitted.values$sd[id.prd]
```

Finally, we plot the results:

``` r
# Plot the predictions
pred_df <- data.frame(x1 = coord.prd[,1],
                      x2 = coord.prd[,2],
                      mean = m.prd,
                      sd = sd.prd)

ggplot(pred_df, aes(x = x1, y = x2, fill = mean)) +
  geom_raster() +
  scale_fill_viridis()
```

![](rspde_inla_files/figure-html/plot_pred-1.png)

Then, the std. deviations:

``` r
ggplot(pred_df, aes(x = x1, y = x2, fill = sd)) +
  geom_raster() + scale_fill_viridis()
```

![](rspde_inla_files/figure-html/plot_pred_sd-1.png)

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
replicates of a latent spatial field $x(\mathbf{s})$, observed at the
same $m$ locations, $\{\mathbf{s}_{1},\ldots,\mathbf{s}_{m}\}$, for each
replicate. For each $i = 1,\ldots,m,$ we have

$$\begin{aligned}
y_{i} & {= x_{1}\left( \mathbf{s}_{i} \right) + \varepsilon_{i},} \\
\vdots & {= \vdots} \\
y_{i + 29m} & {= x_{30}\left( \mathbf{s}_{i} \right) + \varepsilon_{i + 29m},}
\end{aligned}$$

where $\varepsilon_{1},\ldots,\varepsilon_{30m}$ are iid normally
distributed with mean 0 and standard deviation 0.1.

We use the basis function representation of $x( \cdot )$ to define the
$A$ matrix linking the point locations to the mesh. We also need to
account for the fact that we have 30 replicates at the same locations.
To this end, the $A$ matrix we need can be generated by
[`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
function. The reason being that we are sampling $x( \cdot )$ directly
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

![](rspde_inla_files/figure-html/unnamed-chunk-2-1.png)

We then compute the $A$ matrix, which is needed for simulation, and
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

Notice that for the simulated data, we should use the $A$ matrix from
[`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
function instead of the
[`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md)
function.

We will now simulate a latent process with standard deviation
$\sigma = 1$ and range $0.1$. We will use $\nu = 0.5$ so that the model
has an exponential covariance function. To this end we create a model
object with the
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
`operator_information` object. We then obtain the simulated data $y$ by
connecting with the $A$ matrix and adding the gaussian noise.

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

    ## Warning: Removed 7648 rows containing missing values or values outside the scale range
    ## (`geom_raster()`).

![](rspde_inla_files/figure-html/unnamed-chunk-6-1.png)

### Fitting the R-INLA rSPDE model

Let us then use the rational SPDE approach to fit the data.

We begin by creating the $A$ matrix and index with replicates, and the
`inla.stack` object. It is important to notice that since we have
replicates we should provide the `index` and `repl` arguments for
[`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md)
function, and also the argument `n.repl` in
[`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
function. They behave identically as in their
[`R-INLA`](https://www.r-inla.org)’s counterparts, namely,
[`inla.spde.make.A()`](https://rdrr.io/pkg/INLA/man/inla.spde.make.A.html)
and `inla.make.index()`.

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
```

We now create the model object.

``` r
rspde_model.rep <- rspde.matern(mesh = mesh_2d, parameterization = "spde") 
```

Finally, we create the formula and fit. It is extremely important not to
forget the `replicate` argument when building the formula as
[`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) function will not
produce warning and might fit some meaningless model.

``` r
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
      list(A = inla.stack.A(st.dat.rep)),
    num.threads = "1:1"
  )
```

We can get the summary:

``` r
summary(rspde_fit.rep)
```

    ## Time used:
    ##     Pre = 0.184, Running = 50, Post = 0.614, Total = 50.8 
    ## Random effects:
    ##   Name     Model
    ##     field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                          mean    sd 0.025quant 0.5quant
    ## Precision for the Gaussian observations 98.14 4.874      89.23    97.90
    ## Theta1 for field                        -8.02 4.832     -19.42    -7.26
    ## Theta2 for field                         3.94 0.785       2.85     3.82
    ## Theta3 for field                         2.05 2.692      -1.66     1.63
    ##                                         0.975quant   mode
    ## Precision for the Gaussian observations     108.41 97.173
    ## Theta1 for field                             -1.36 -3.865
    ## Theta2 for field                              5.79  3.266
    ## Theta3 for field                              8.40 -0.263
    ## 
    ## Marginal log-Likelihood:  -4596.60 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

and the summary in the user’s scale:

``` r
result_fit_rep <- rspde.result(rspde_fit.rep, "field", rspde_model.rep)
```

    ## Warning in rspde.result(rspde_fit.rep, "field", rspde_model.rep): the mean or
    ## mode of nu is very close to nu.upper.bound, please consider increasing
    ## nu.upper.bound, and refitting the model.

``` r
summary(result_fit_rep)
```

    ##             mean         sd  0.025quant    0.5quant 0.975quant        mode
    ## tau    0.0308873  0.0999215 3.14880e-13 8.18709e-04   0.268554 3.14880e-13
    ## kappa 74.2683000 94.1564000 1.72450e+01 4.42538e+01 319.828000 2.20023e+01
    ## nu     1.4342000  0.5627370 3.13424e-01 1.64775e+00   1.999510 1.99999e+00

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

    ##   parameter        true        mean         mode
    ## 1       tau  0.08920621  0.03088728 3.148803e-13
    ## 2     kappa 20.00000000 74.26831357 2.200230e+01
    ## 3        nu  0.50000000  1.43420329 1.999988e+00

## An example with a non-stationary model

It is also possible to consider models in which $\sigma$ (std.
deviation) and $\rho$ (range parameter) are non-stationary. One can also
use the parameterization in terms of the SPDE parameters $\kappa$ and
$\tau$.

An example of such a model is given in the vignette [inlabru
implementation of the rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inlabru.md).

## Further options of the `rSPDE`-`INLA` implementation

We will now discuss some of the arguments that were introduced in our
[`R-INLA`](https://www.r-inla.org) implementation of the rational
approximation that are not present in
[`R-INLA`](https://www.r-inla.org)’s standard SPDE implementation.

In each case we will provide an illustrative example.

### Changing the upper bound for the smoothness parameter

When we fit a
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
model we need to provide an upper bound for the smoothness parameter
$\nu$. The reason for that is that the sparsity of the precision matrix
should be kept fixed during [`R-INLA`](https://www.r-inla.org)’s
estimation and the higher the value of $\nu$ the denser the precision
matrix gets.

This means that the higher the value of $\nu$, the higher the
computational cost to fit the model. Therefore, ideally, want to choose
an upper bound for $\nu$ as small as possible.

To change the value of the upper bound for the smoothness parameter, we
must change the argument `nu.upper.bound`. The default value for
`nu.upper.bound` is 2. Other common choices for `nu.upper.bound` are 1,
3 and 4.

It is clear from the discussion above that the smaller the value of
`nu.upper.bound` the faster the estimation procedure will be.

However, if we choose a value of `nu.upper.bound` which is too low, the
“correct” value of $\nu$ might not belong to the interval
$\left( 0,\nu_{UB} \right)$, where $\nu_{UB}$ is the value of
`nu.upper.bound`. Hence, one might be forced to increase
`nu.upper.bound` and estimate again, which, obviously will increase the
computational cost as we will need to do more than one estimation.

Let us illustrate by considering the same model we considered above for
the precipitation in Paraná region in Brazil and consider
`nu.upper.bound` equal to 4, which is generally a good choice for
`nu.upper.bound`.

We simply use the function
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
with the argument `nu.upper.bound` set to 4:

``` r
rspde_model_2 <- rspde.matern(mesh = prmesh, nu.upper.bound = 4)
```

Since we are considering the default `rspde.order`, the $A$ matrix and
the mesh index objects are the same as the previous ones. Let us then
update the formula and fit the model:

``` r
f.s.2 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_2)

rspde_fit_2 <- inla(f.s.2,
  family = "Gamma", data = inla.stack.data(stk.dat),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
```

Let us see the summary of the fit:

``` r
summary(rspde_fit_2)
```

    ## Time used:
    ##     Pre = 0.173, Running = 29.2, Post = 0.0275, Total = 29.4 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.051      1.841    1.941      2.041 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                     mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations    14.415     1.20      12.17
    ## Precision for seaDist                          19028.097 38117.60     761.28
    ## Theta1 for field                                   0.707     4.79      -8.58
    ## Theta2 for field                                   0.396     5.64     -10.88
    ## Theta3 for field                                  -3.358     9.42     -22.19
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations   14.372      16.90   14.301
    ## Precision for seaDist                          8536.754  103468.99 1877.130
    ## Theta1 for field                                  0.658      10.29    0.448
    ## Theta2 for field                                  0.454      11.33    0.701
    ## Theta3 for field                                 -3.261      14.90   -2.850
    ## 
    ## Marginal log-Likelihood:  -1252.09 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

Let us compare with the cost from the previous fit, with the default
value of `nu.upper.bound` of 2:

``` r
# nu.upper.bound = 2
rspde_fit$cpu.used
```

    ##        Pre    Running       Post      Total 
    ## 0.38407564 5.76345634 0.03418279 6.18171477

``` r
# nu.upper.bound = 4
rspde_fit_2$cpu.used
```

    ##         Pre     Running        Post       Total 
    ##  0.17318273 29.15342498  0.02749324 29.35410094

We can see that the fit for `nu.upper.bound` equal to 2 was considerably
faster.

Finally, let us get the result results for the field and see the
estimate of $\nu$:

``` r
result_fit_2 <- rspde.result(rspde_fit_2, "field", rspde_model_2)
summary(result_fit_2)
```

    ##              mean          sd  0.025quant 0.5quant  0.975quant        mode
    ## tau   20196.30000 3.78019e+05 2.02290e-04 1.911370 27607.70000 4.00281e-07
    ## kappa 98151.40000 2.04325e+06 2.08296e-05 1.601210 79317.00000 8.31032e-09
    ## nu        1.46665 1.75488e+00 2.33066e-15 0.150398     3.99999 2.33066e-15

### Changing the order of the rational approximation

To change the order of the rational approximation all we have to do is
set the argument `rspde.order` to the desired value. The current
available possibilities are `1,2,3`,…, up to `8`.

The higher the order of the rational approximation, the more accurate
the results will be, however, the higher the computational cost will be.

The default `rspde.order` of 1 is generally a good choice, fast, and
reasonably accurate. See the vignette [Rational approximation with the
rSPDE package](https://davidbolin.github.io/rSPDE/articles/rspde_cov.md)
for further details on the order of the rational approximation and some
comparison with the Matérn covariance.

Let us fit the above model with the covariance-based rational
approximation of order `3`. Since we are changing the order of the
rational approximation, that is, we are changing the `rspde.order`
argument, we need to recompute the $A$ matrix and the mesh index.
Therefore, we proceed as follows:

- We build a new model:

``` r
rspde_model_order_3 <- rspde.matern(mesh = prmesh, 
  rspde.order = 3
)
```

- We create a new $A$ matrix:

``` r
Abar_3 <- rspde.make.A(mesh = prmesh, loc = coords, rspde.order = 3)
```

- We create a new index:

``` r
mesh.index.3 <- rspde.make.index(
  name = "field", mesh = prmesh,
  rspde.order = 3
)
```

Now the remaining steps are the same as before:

``` r
stk.dat.3 <- inla.stack(
  data = list(y = Y), A = list(Abar_3, 1), tag = "est",
  effects = list(
    c(
      mesh.index.3
    ),
    list(
      long = inla.group(coords[, 1]),
      lat = inla.group(coords[, 2]),
      seaDist = inla.group(seaDist),
      Intercept = 1
    )
  )
)

f.s.3 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_order_3)

rspde_fit_order_3 <- inla(f.s.3,
  family = "Gamma", data = inla.stack.data(stk.dat.3),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat.3), compute = TRUE),
  num.threads = "1:1"
)
```

Let us see the summary:

``` r
summary(rspde_fit_order_3)
```

    ## Time used:
    ##     Pre = 0.179, Running = 19.1, Post = 0.0344, Total = 19.3 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.041       1.86    1.941      2.022 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.46    1.040     12.515
    ## Precision for seaDist                          7501.98 4203.960   2378.419
    ## Theta1 for field                                 -1.79    1.281     -4.599
    ## Theta2 for field                                  1.60    0.405      0.872
    ## Theta3 for field                                  0.29    1.043     -1.489
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations   14.426   1.66e+01   14.361
    ## Precision for seaDist                          6534.777   1.83e+04 4978.996
    ## Theta1 for field                                 -1.696   3.83e-01   -1.215
    ## Theta2 for field                                  1.582   2.46e+00    1.485
    ## Theta3 for field                                  0.214   2.57e+00   -0.166
    ## 
    ## Marginal log-Likelihood:  -1255.48 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

We can see in the above summary that the computational cost
significantly increased. Let us compare the cost of having
`rspde.order=3` with the cost of having `rspde.order=1`:

``` r
# order = 1
rspde_fit$cpu.used
```

    ##        Pre    Running       Post      Total 
    ## 0.38407564 5.76345634 0.03418279 6.18171477

``` r
# order = 3
rspde_fit_order_3$cpu.used
```

    ##         Pre     Running        Post       Total 
    ##  0.17929173 19.13547301  0.03442478 19.34918952

One can check the order of the rational approximation by using the
[`rational.order()`](https://davidbolin.github.io/rSPDE/reference/rational.order.md)
function. It also allows another way to change the order of the rational
order, by using the corresponding `rational.order<-()` function.

The
[`rational.order()`](https://davidbolin.github.io/rSPDE/reference/rational.order.md)
and `rational.order<-()` functions can be applied to the `inla.rspde`
object, to the `A` matrix and to the `index` objects.

Here to check the models:

``` r
rational.order(rspde_model)
```

    ## [1] 1

``` r
rational.order(rspde_model_order_3)
```

    ## [1] 3

Here to check the `A` matrices:

``` r
rational.order(Abar)
```

    ## [1] 1

``` r
rational.order(Abar_3)
```

    ## [1] 3

Here to check the indexes:

``` r
rational.order(mesh.index)
```

    ## [1] 1

``` r
rational.order(mesh.index.3)
```

    ## [1] 3

Let us now change the order of the `rspde_model` object to be 2:

``` r
rational.order(rspde_model) <- 2
rational.order(Abar) <- 2
rational.order(mesh.index) <- 2
```

Let us fit this new model:

``` r
 f.s.2 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model)

stk.dat.2 <- inla.stack(
  data = list(y = Y), A = list(Abar, 1), tag = "est",
  effects = list(
    c(
      mesh.index
    ),
    list(
      long = inla.group(coords[, 1]),
      lat = inla.group(coords[, 2]),
      seaDist = inla.group(seaDist),
      Intercept = 1
    )
  )
)

rspde_fit_order_2 <- inla(f.s.2,
  family = "Gamma", data = inla.stack.data(stk.dat.2),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat.2), compute = TRUE),
  num.threads = "1:1"
)
```

Here is the summary:

``` r
summary(rspde_fit_order_2)
```

    ## Time used:
    ##     Pre = 0.187, Running = 10.1, Post = 0.0311, Total = 10.3 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.041      1.861    1.941      2.022 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.47    1.042      12.51
    ## Precision for seaDist                          7307.10 4132.278    2160.07
    ## Theta1 for field                                 -3.07    2.077      -7.85
    ## Theta2 for field                                  1.83    0.494       1.03
    ## Theta3 for field                                  1.29    1.643      -1.21
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations    14.43   1.66e+01   14.377
    ## Precision for seaDist                           6382.41   1.79e+04 4814.421
    ## Theta1 for field                                  -2.80   6.80e-02   -1.457
    ## Theta2 for field                                   1.79   2.94e+00    1.543
    ## Theta3 for field                                   1.08   5.07e+00    0.025
    ## 
    ## Marginal log-Likelihood:  -1255.21 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

### Estimating models with fixed smoothness

We can fix the smoothness, say $\nu$, of the model by providing a
non-NULL positive value for `nu`.

When the smoothness, $\nu$, is fixed, we can have two possibilities:

- $\alpha = \nu + d/2$ is integer;

- $\alpha = \nu + d/2$ is not integer.

The first case, i.e., when $\alpha$ is integer, has less computational
cost. Furthermore, the $A$ matrix is different than the $A$ matrix for
the non-integer $\alpha$.

The $A$ matrix is the same for all values of $\nu$ such that $\alpha$ is
integer. So, the $A$ matrix for these cases only need to be computed
once. The same holds for the `index` obtained from the
[`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
function.

In the second case the $A$ matrix only depends on the order of the
rational approximation and not on $\nu$. Therefore, if the matrix $A$
has already been computed for some `rspde.order`, then the $A$ matrix
will be same for all the values of $\nu$ such that $\alpha$ is
non-integer for that `rspde.order`. The same holds for the `index`
obtained from the
[`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
function.

If $\nu$ is fixed, we have that the parameters returned by
[`R-INLA`](https://www.r-inla.org) are \$\$\kappa =
\exp(\theta_1)\quad\hbox{and}\quad\tau = \exp(\theta_2).\$\$ We will now
provide illustrations for both scenarios. It is also noteworthy that
just as for the case in which we estimate $\nu$, we can also change the
order of the rational approximation by changing the value of
`rspde.order`. For both illustrations with fixed $\nu$ below, we will
only consider the order of the rational approximation of 1, that is, the
default order.

#### Estimating models with fixed smoothness and non-integer $\alpha$

Recall that:
$$\nu = \nu_{UB}(\frac{\exp\left( \theta_{3} \right)}{1 + \exp\left( \theta_{3} \right)}).$$
Thus, to illustrate, let us consider a fixed $\nu$ given by the mean of
$\nu$ obtained from the first model we considered in this vignette,
namely, the fit given by `rspde_fit`, which is approximately
$\nu = 1.21$.

Notice that for this $\nu$, the value of $\alpha$ is non-integer, so we
can use the $A$ matrix and the index of the first fitted model, which is
also of order 2.

Therefore, all we have to do is build a new model in which we set `nu`
to `1.21`:

``` r
rspde_model_fix <- rspde.matern(mesh = prmesh, nu = 1.21)
```

Let us now fit the model:

``` r
f.s.fix <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_fix)

rspde_fix <- inla(f.s.fix,
  family = "Gamma", data = inla.stack.data(stk.dat),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
```

Here we have the summary:

``` r
summary(rspde_fix)
```

    ## Time used:
    ##     Pre = 0.174, Running = 4.02, Post = 0.0256, Total = 4.22 
    ## Fixed effects:
    ##            mean   sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.04      1.863    1.941       2.02 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.44    1.042      12.49
    ## Precision for seaDist                          7462.98 4213.256    2402.15
    ## Theta1 for field                                 -1.89    0.359      -2.60
    ## Theta2 for field                                  1.60    0.268       1.08
    ##                                                0.5quant 0.975quant    mode
    ## Precision-parameter for the Gamma observations    14.41      16.59   14.35
    ## Precision for seaDist                           6476.38   18387.24 4923.22
    ## Theta1 for field                                  -1.88      -1.18   -1.88
    ## Theta2 for field                                   1.60       2.13    1.60
    ## 
    ## Marginal log-Likelihood:  -1256.21 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

Now, the summary in the original scale:

``` r
result_fix <- rspde.result(rspde_fix, "field", rspde_model_fix)
summary(result_fix)
```

    ##           mean        sd 0.025quant 0.5quant 0.975quant     mode
    ## tau   0.161772 0.0591347  0.0749454  0.15214   0.304701 0.134069
    ## kappa 5.143310 1.3943700  2.9531600  4.95768   8.396250 4.604790

#### Estimating models with fixed smoothness and integer $\alpha$

Since we are in dimension $d = 2$, and $\nu > 0$, the smallest value of
$\nu$ that makes $\alpha = \nu + 1$ an integer is $\nu = 1$. This value
is also close to the estimated mean of the first model we fitted and
enclosed by the posterior marginal density of $\nu$ for the first fit.
Therefore, let us fit the model with $\nu = 1$.

To this end we need to compute a new $A$ matrix:

``` r
Abar.int <- rspde.make.A(
  mesh = prmesh, loc = coords,
  nu = 1
)
```

a new index:

``` r
mesh.index.int <- rspde.make.index(
  name = "field", mesh = prmesh,
  nu = 1
)
```

create a new model (remember to set `nu=1`):

``` r
rspde_model_fix_int1 <- rspde.matern(mesh = prmesh,
  nu = 1)
```

The remaining is standard:

``` r
stk.dat.int <- inla.stack(
  data = list(y = Y), A = list(Abar.int, 1), tag = "est",
  effects = list(
    c(
      mesh.index.int
    ),
    list(
      long = inla.group(coords[, 1]),
      lat = inla.group(coords[, 2]),
      seaDist = inla.group(seaDist),
      Intercept = 1
    )
  )
)

f.s.fix.int.1 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_fix_int1)

rspde_fix_int_1 <- inla(f.s.fix.int.1,
  family = "Gamma",
  data = inla.stack.data(stk.dat.int), verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(
    A = inla.stack.A(stk.dat.int),
    compute = TRUE
  ),
  num.threads = "1:1"
)
```

Let us check the summary:

``` r
summary(rspde_fix_int_1)
```

    ## Time used:
    ##     Pre = 0.179, Running = 1.05, Post = 0.0327, Total = 1.26 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.942 0.041       1.86    1.942      2.023 1.942   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.43    1.041     12.488
    ## Precision for seaDist                          7624.79 4460.558   2402.584
    ## Theta1 for field                                 -1.40    0.327     -2.052
    ## Theta2 for field                                  1.52    0.292      0.951
    ##                                                0.5quant 0.975quant    mode
    ## Precision-parameter for the Gamma observations    14.39     16.583   14.32
    ## Precision for seaDist                           6549.97  19265.727 4903.23
    ## Theta1 for field                                  -1.40     -0.765   -1.39
    ## Theta2 for field                                   1.52      2.100    1.51
    ## 
    ## Marginal log-Likelihood:  -1255.98 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

and check the summary in the user’s scale:

``` r
rspde_result_int <- rspde.result(rspde_fix_int_1, "field", rspde_model_fix_int1)
summary(rspde_result_int)
```

    ##           mean        sd 0.025quant 0.5quant 0.975quant     mode
    ## tau   0.259929 0.0857263   0.129291 0.247525   0.462895 0.223491
    ## kappa 4.764460 1.4170800   2.600690 4.554970   8.124940 4.150690

### Changing the priors

We begin by recalling that the fitted `rSPDE` model in
[`R-INLA`](https://www.r-inla.org) contains the parameters
$\text{Theta1}$, $\text{Theta2}$ and $\text{Theta3}$. Let, again,
$\theta_{1} = \text{Theta1}$, $\theta_{2} = \text{Theta2}$ and
$\theta_{3} = \text{Theta3}$. In terms of the SPDE
$$\left( \kappa^{2}I - \Delta \right)^{\alpha/2}(\tau u) = \mathcal{W},$$
where $\alpha = \nu + d/2$.

We also have the range parameter $\rho = \frac{\sqrt{8\nu}}{\kappa}$ and
the standard deviation
$\sigma = \sqrt{\frac{\Gamma(\nu)}{\tau^{2}\kappa^{2\nu}(4\pi)^{d/2}\Gamma(\nu + d/2)}}$.

#### Changing the priors of $\tau$ and $\kappa$

We begin by dealing with $\tau$ and $\kappa$.

We have that
$$\tau = \exp\left( \theta_{1} \right),\quad\kappa = \exp\left( \theta_{2} \right).$$
The
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function assumes a lognormal prior distribution for $\tau$ and $\kappa$.
This prior distribution is obtained by assuming that $\theta_{1}$ and
$\theta_{2}$ follow normal distributions. By default we assume
$\theta_{1}$ and $\theta_{2}$ to be independent and to follow normal
distributions
$\theta_{1} \sim N\left( \log\left( \tau_{0} \right),10 \right)$ and
$\theta_{2} \sim N\left( \log\left( \kappa_{0} \right),10 \right)$.

$\kappa_{0}$ is suitably defined in terms of the mesh and $\tau_{0}$ is
defined in terms of $\kappa_{0}$ and on the prior of the smoothness
parameter.

If one wants to define the prior
$$\theta_{1} \sim N\left( \text{mean\_theta\_1},\text{sd\_theta\_1} \right),$$
one can simply set the argument
`prior.tau = list(meanlog=mean_theta_1, sdlog=sd_theta_1)`. Analogously,
to define the prior
$$\theta_{2} \sim N\left( \text{mean\_theta\_2},\text{sd\_theta\_2} \right),$$
one can set the argument
`prior.kappa = list(meanlog=mean_theta_2, sdlog=sd_theta_2)`.

It is important to mention that, by default, the initial values of
$\tau$ and $\kappa$ are $\exp\left( \text{mean\_theta\_1} \right)$ and
$\exp\left( \text{mean\_theta\_2} \right)$, respectively. So, if the
user does not change these parameters, and also does not change the
initial values, the initial values of $\tau$ and $\kappa$ will be,
respectively, $\tau_{0}$ and $\kappa_{0}$.

If one sets `prior.tau = list(meanlog=mean_theta_1)`, the prior for
$\theta_{1}$ will be
$$\theta_{1} \sim N\left( \text{mean\_theta\_1},1 \right),$$ whereas, if
one sets, `prior.tau = list(sdlog=sd_theta_1)`, the prior will be
$$\theta_{1} \sim N\left( \log\left( \tau_{0} \right),\text{sd\_theta\_1} \right).$$
Analogously, if one sets `prior.kappa = list(meanlog=mean_theta_2)`, the
prior for $\theta_{2}$ will be
$$\theta_{2} \sim N\left( \text{mean\_theta\_2},1 \right),$$ whereas, if
one sets, `prior.kappa = list(sdlog=sd_theta_2)`, the prior will be
$$\theta_{2} \sim N\left( \log\left( \kappa_{0} \right),\text{sd\_theta\_2} \right).$$

#### Changing the priors of $\rho$ (range) and $\sigma$ (std. dev.)

Let us now consider the priors for the range, $\rho$, and std.
deviation, $\sigma$. This parameterization is used with the argument
`parameterization = "matern"`, which is the default.

In this case, we have that
$$\sigma = \exp\left( \theta_{1} \right),\quad\rho = \exp\left( \theta_{2} \right).$$
We have two options for prior. By default, which is the option
`prior.theta.param = "theta"`, the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function assumes a lognormal prior distribution for $\sigma$ and $\rho$.
This prior distribution is obtained by assuming that $\theta_{1}$ and
$\theta_{2}$ follow normal distributions. By default we assume
$\theta_{1}$ and $\theta_{2}$ to be independent and to follow normal
distributions
$\theta_{1} \sim N\left( \log\left( \sigma_{0} \right),10 \right)$ and
$\theta_{2} \sim N\left( \log\left( \rho_{0} \right),10 \right)$.

$\rho_{0}$ is suitably defined in terms of the mesh and $\sigma_{0}$ is
defined in terms of $\rho_{0}$ and on the prior of the smoothness
parameter.

Similarly to the priors of $\tau$ and $\kappa$, if one wants to define
the prior
$$\theta_{1} \sim N\left( \text{mean\_theta\_1},\text{sd\_theta\_1} \right),$$
one can simply set the argument
`prior.tau = list(meanlog=mean_theta_1, sdlog=sd_theta_1)`. Analogously,
to define the prior
$$\theta_{2} \sim N\left( \text{mean\_theta\_2},\text{sd\_theta\_2} \right),$$
one can set the argument
`prior.kappa = list(meanlog=mean_theta_2, sdlog=sd_theta_2)`.

Another option is to set `prior.theta.param = "spde"`. In this case, a
change of variables is performed. So, we assume a lognormal prior on
$\tau$ and $\kappa$. Then, by the relations
$\rho = \frac{\sqrt{8\nu}}{\kappa}$ and
$\sigma = \sqrt{\frac{\Gamma(\nu)}{\tau^{2}\kappa^{2\nu}(4\pi)^{d/2}\Gamma(\nu + d/2)}}$,
we obtain a prior for $\rho$ and $\sigma$.

#### Changing the prior of $\nu$

Finally, let us consider the smoothness parameter $\nu$.

By default, we assume that $\nu$ follows a beta distribution on the
interval $\left( 0,\nu_{UB} \right)$, where $\nu_{UB}$ is the upper
bound for $\nu$, with mean $\nu_{0} = \min\{ 1,\nu_{UB}/2\}$ and
variance
$\frac{\nu_{0}\left( \nu_{UB} - \nu_{0} \right)}{1 + \phi_{0}}$, and we
call $\phi_{0}$ the precision parameter, whose default value is
$$\phi_{0} = \max\{\frac{\nu_{UB}}{\nu_{0}},\frac{\nu_{UB}}{\nu_{UB} - \nu_{0}}\} + \phi_{inc}.$$
The parameter $\phi_{inc}$ is an increment to ensure that the prior beta
density has boundary values equal to zero (where the boundary values are
defined either by continuity or by limits). The default value of
$\phi_{inc}$ is 1. The value of $\phi_{inc}$ can be changed by changing
the argument `nu.prec.inc` in the
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
function. The higher the value of $\phi_{inc}$ (that is, the value of
`nu.prec.inc`) the more informative the prior distribution becomes.

Let us denote a beta distribution with support on
$\left( 0,\nu_{UB} \right)$, mean $\mu$ and precision parameter $\phi$
by $\mathcal{B}_{\nu_{UB}}(\mu,\phi)$.

If we want $\nu$ to have a prior
$$\nu \sim \mathcal{B}_{\nu_{UB}}\left( \text{nu\_1},\text{prec\_1} \right),$$
one simply needs to set `prior.nu = list(mean=nu_1, prec=prec_1)`. If
one sets `prior.nu = list(mean=nu_1)`, then $\nu$ will have prior
$$\nu \sim \mathcal{B}_{\nu_{UB}}\left( \text{nu\_1},\phi_{1} \right),$$
where
$$\phi_{1} = \max\{\frac{\nu_{UB}}{\text{nu\_1}},\frac{\nu_{UB}}{\nu_{UB} - \text{nu\_1}}\} + \text{nu.prec.inc}.$$

Of one sets `prior.nu = list(prec=prec_1)`, then $\nu$ will have prior
$$\nu \sim \mathcal{B}_{\nu_{UB}}\left( \nu_{0},\text{prec\_1} \right).$$
It is also noteworthy that we have that, in terms of
[`R-INLA`](https://www.r-inla.org)’s parameters,

$$\nu = \nu_{UB}(\frac{\exp\left( \theta_{3} \right)}{1 + \exp\left( \theta_{3} \right)}).$$

It is important to mention that, by default, if a beta prior
distribution is chosen for the smoothness parameter $\nu$, then the
initial value of $\nu$ is the mean of the prior beta distribution. So,
if the user does not change this parameter, and also does not change the
initial value, the initial value of $\nu$ will be
$\min\{ 1,\nu_{UB}/2\}$.

We also assume that, in terms of [`R-INLA`](https://www.r-inla.org)’s
parameters,
$$\nu = \nu_{UB}(\frac{\exp\left( \theta_{3} \right)}{1 + \exp\left( \theta_{3} \right)}).$$

We can have another possibility of prior distribution for $\nu$, namely,
truncated lognormal distribution. The truncated lognormal distribution
is defined in the following sense. We assume that $\log(\nu)$ has prior
distribution given by a [truncated normal
distribution](https://en.wikipedia.org/wiki/Truncated_normal_distribution#One_sided_truncation_(of_upper_tail))
with support $\left( - \infty,\log\left( \nu_{UB} \right) \right)$,
where $\nu_{UB}$ is the upper bound for $\nu$, with location parameter
$\mu_{0} = \log\left( \nu_{0} \right) = \log(\min\{ 1,\nu_{UB}/2\})$ and
scale parameter $\sigma_{0} = 1$. More precisely, let
$\Phi( \cdot ;\mu,\sigma)$ stand for the cumulative distribution
function (CDF) of a normal distribution with mean $\mu$ and standard
deviation $\sigma$. Then, $\log(\nu)$ has cumulative distribution
function given by
$$F_{\log{(\nu)}}(x) = \frac{\Phi\left( x;\mu_{0},\sigma_{0} \right)}{\Phi\left( \nu_{UB} \right)},\quad x \leq \nu_{UB},$$
and $F_{\log{(\nu)}}(x) = 1$ if $x > \nu_{UB}$. We will call $\mu_{0}$
and $\sigma_{0}$ the log-location and log-scale parameters of $\nu$,
respectively, and we say that $\log(\nu)$ follows a truncated normal
distribution with location parameter $\mu_{0}$ and scale parameter
$\sigma_{0}$.

To change the prior distribution of $\nu$ to the truncated lognormal
distribution, we need to set the argument `prior.nu.dist="lognormal"`.

To change these parameters in the prior distribution to, say, `log_nu_1`
and `log_sigma_1`, one can simply set
`prior.nu = list(loglocation=log_nu_1, logscale=sigma_1)`.

If one sets `prior.nu = list(loglocation=log_nu_1)`, the prior for
$\theta_{3}$ will be a truncated normal normal distribution with
location parameter `log_nu_1` and scale parameter `1`. Analogously, if
one sets, `prior.nu = list(logscale=sigma_1)`, the prior for
$\theta_{3}$ will be a truncated normal distribution with location
parameter $\log\left( \nu_{0} \right) = \log(\min\{ 1,\nu_{UB}/2\})$ and
scale parameter `sigma_1`.

It is important to mention that, by default, if a truncated lognormal
prior distribution is chosen for the smoothness parameter $\nu$, then
the initial value of $\nu$ is the exponential of the log-location
parameter of $\nu$. So, if the user does not change this parameter, and
also does not change the initial value, the initial value of $\nu$ will
be $\min\{ 1,\nu_{UB}/2\}$.

Let us consider an example with the same dataset used in the first model
of this vignette where we change the prior distribution of $\nu$ from
beta to lognormal.

``` r
rspde_model_beta <- rspde.matern(mesh = prmesh, prior.nu.dist = "lognormal")
```

Since we did not change `rspde.order` and are not fixing $\nu$, we can
use the same $A$ matrix and same index from the first example.

Therefore, all we have to do is update the formula and fit the model:

``` r
f.s.beta <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_beta)

rspde_fit_beta <- inla(f.s.beta,
  family = "Gamma", data = inla.stack.data(stk.dat),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
```

We have the summary:

``` r
summary(rspde_fit_beta)
```

    ## Time used:
    ##     Pre = 0.172, Running = 6.56, Post = 0.0262, Total = 6.76 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.942 0.041       1.86    1.942      2.023 1.942   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                    mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.433    1.038     12.473
    ## Precision for seaDist                          7693.167 4499.817   2404.657
    ## Theta1 for field                                 -1.783    1.362     -4.800
    ## Theta2 for field                                  1.570    0.402      0.849
    ## Theta3 for field                                  0.316    1.157     -1.610
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations   14.404   1.66e+01   14.365
    ## Precision for seaDist                          6613.108   1.94e+04 4949.561
    ## Theta1 for field                                 -1.665   4.77e-01   -1.079
    ## Theta2 for field                                  1.550   2.42e+00    1.448
    ## Theta3 for field                                  0.217   2.87e+00   -0.271
    ## 
    ## Marginal log-Likelihood:  -1255.65 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

Also, we can have the summary in the user’s scale:

``` r
result_fit_beta <- rspde.result(rspde_fit_beta, "field", rspde_model_beta)
summary(result_fit_beta)
```

    ##          mean       sd 0.025quant 0.5quant 0.975quant      mode
    ## tau   0.35168 0.454859  0.0085236 0.196327    1.59893 0.0139569
    ## kappa 5.22068 2.305360  2.3495500 4.673070   11.19560 3.7982100
    ## nu    1.10806 0.448041  0.3350500 1.092740    1.89001 0.7850190

and the plot of the posterior marginal densities

``` r
posterior_df_fit_beta <- gg_df(result_fit_beta)

ggplot(posterior_df_fit_beta) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](rspde_inla_files/figure-html/unnamed-chunk-32-1.png)

### Changing the starting values

The starting values to be used by [`R-INLA`](https://www.r-inla.org)’s
optimization algorithm can be changed by setting the arguments
`start.ltau`, `start.lkappa` and `start.nu`.

- `start.ltau` will be the initial value for $\log(\tau)$, that is, the
  logarithm of $\tau$.

- `start.lkappa` will be the inital value for $\log(\kappa)$, that is,
  the logarithm of $\kappa$.

- `start.nu` will be the initial value for $\nu$. Notice that here the
  initial value is *not* on the log scale.

One can change the initial value of one or more parameters.

For instance, let us consider the example with precipitation data,
`rspde.order=3`, but change the initial values to the ones close to the
fitted value when considering the default `rspde.order` (which is 1):

``` r
rspde_model_order_3_start <- rspde.matern(mesh = prmesh, rspde.order = 3,
  nu.upper.bound = 2,
  start.lkappa = result_fit$summary.log.kappa$mean,
  start.ltau = result_fit$summary.log.tau$mean,
  start.nu = min(result_fit$summary.nu$mean, 2 - 1e-5)
)
```

Since we already computed the $A$ matrix and the index for
`rspde.order=3`, all we have to do is to update the formula and fit:

``` r
f.s.3.start <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_order_3_start)

rspde_fit_order_3_start <- inla(f.s.3.start,
  family = "Gamma",
  data = inla.stack.data(stk.dat.3),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(
    A = inla.stack.A(stk.dat.3),
    compute = TRUE
  ),
  num.threads = "1:1"
)
```

We have the summary:

``` r
summary(rspde_fit_order_3_start)
```

    ## Time used:
    ##     Pre = 0.174, Running = 15.9, Post = 0.0336, Total = 16.2 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.038      1.866    1.941      2.015 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.47    1.043     12.496
    ## Precision for seaDist                          7265.31 3877.218   2339.349
    ## Theta1 for field                                 -3.14    0.809     -4.700
    ## Theta2 for field                                  1.82    0.254      1.316
    ## Theta3 for field                                  1.59    1.014     -0.439
    ##                                                0.5quant 0.975quant    mode
    ## Precision-parameter for the Gamma observations    14.44      16.60   14.40
    ## Precision for seaDist                           6423.08   17126.24 4994.32
    ## Theta1 for field                                  -3.15      -1.52   -3.19
    ## Theta2 for field                                   1.83       2.31    1.84
    ## Theta3 for field                                   1.60       3.55    1.65
    ## 
    ## Marginal log-Likelihood:  -1256.40 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

### Changing the type of the rational approximation

We have three rational approximations available. The BRASIL algorithm
[Hofreither (2021)](https://doi.org/10.1007/s11075-020-01042-0), and two
“versions” of the Clenshaw-Lord Chebyshev-Pade algorithm, one with lower
bound zero and another with the lower bound given in [Bolin, Simas, and
Xiong (2023)](https://doi.org/10.1080/10618600.2023.2231051).

The type of rational approximation can be chosen by setting the
`type.rational.approx` argument in the `rspde.matern` function. The
BRASIL algorithm corresponds to the choice `brasil`, the Clenshaw-Lord
Chebyshev pade with zero lower bound and non-zero lower bounds are
given, respectively, by the choices `chebfun` and `chebfunLB`.

Let us fit a model assigning a `brasil` rational approximation. We will
consider a model with the order of the rational approximation being 1:

``` r
rspde_model_brasil <- rspde.matern(prmesh, 
              type.rational.approx = "brasil")

f.s.brasil <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_brasil)

rspde_fit_order_1_brasil <- inla(f.s.brasil,
  family = "Gamma", data = inla.stack.data(stk.dat),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
  num.threads = "1:1"
)
```

Let us get the summary:

``` r
summary(rspde_fit_order_1_brasil)
```

    ## Time used:
    ##     Pre = 0.188, Running = 5.43, Post = 0.0268, Total = 5.65 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.942 0.042       1.86    1.942      2.023 1.942   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                    mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.434    1.038     12.490
    ## Precision for seaDist                          7609.809 4339.961   2390.076
    ## Theta1 for field                                 -2.379    2.055     -7.053
    ## Theta2 for field                                  1.700    0.534      0.812
    ## Theta3 for field                                  0.883    1.830     -1.989
    ##                                                0.5quant 0.975quant    mode
    ## Precision-parameter for the Gamma observations   14.400   1.66e+01   14.34
    ## Precision for seaDist                          6594.899   1.88e+04 4987.12
    ## Theta1 for field                                 -2.142   8.39e-01   -1.00
    ## Theta2 for field                                  1.655   2.88e+00    1.43
    ## Theta3 for field                                  0.674   5.04e+00   -0.33
    ## 
    ## Marginal log-Likelihood:  -1255.05 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

Finally, similarly to the order of the rational approximation, one can
check the order with the
[`rational.type()`](https://davidbolin.github.io/rSPDE/reference/rational.type.md)
function, and assign a new type with the `rational.type<-()` function.

``` r
rational.type(rspde_model)
```

    ## [1] "brasil"

``` r
rational.type(rspde_model_brasil)
```

    ## [1] "brasil"

Let us change the type of the rational approximation on the model with
rational approximation of order 3:

``` r
rational.type(rspde_model_order_3) <- "brasil"

f.s.3 <- y ~ -1 + Intercept + f(seaDist, model = "rw1") +
  f(field, model = rspde_model_order_3)

rspde_fit_order_3_brasil <- inla(f.s.3,
  family = "Gamma", data = inla.stack.data(stk.dat.3),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb"),
  control.predictor = list(A = inla.stack.A(stk.dat.3), compute = TRUE),
  num.threads = "1:1"
)
```

Let us get the summary:

``` r
summary(rspde_fit_order_3_brasil)
```

    ## Time used:
    ##     Pre = 0.18, Running = 17.5, Post = 0.0341, Total = 17.7 
    ## Fixed effects:
    ##            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 1.941 0.041      1.861    1.941      2.022 1.941   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     seaDist RW1 model
    ##    field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                                   mean       sd 0.025quant
    ## Precision-parameter for the Gamma observations   14.46    1.040     12.516
    ## Precision for seaDist                          7472.06 4173.502   2377.984
    ## Theta1 for field                                 -4.41    3.521    -12.619
    ## Theta2 for field                                  2.17    0.846      0.899
    ## Theta3 for field                                  2.42    2.866     -1.738
    ##                                                0.5quant 0.975quant     mode
    ## Precision-parameter for the Gamma observations    14.43   1.66e+01   14.361
    ## Precision for seaDist                           6513.73   1.82e+04 4970.466
    ## Theta1 for field                                  -3.92   6.87e-01   -1.438
    ## Theta2 for field                                   2.07   4.13e+00    1.521
    ## Theta3 for field                                   2.02   9.10e+00    0.012
    ## 
    ## Marginal log-Likelihood:  -1254.71 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

## References

Bolin, David, Alexandre B. Simas, and Zhen Xiong. 2023.
“Covariance-Based Rational Approximations of Fractional SPDEs for
Computationally Efficient Bayesian Inference.” *Journal of Computational
and Graphical Statistics*.
<https://doi.org/10.1080/10618600.2022.2139648>.

Hofreither, Clemens. 2021. “An Algorithm for Best Rational Approximation
Based on Barycentric Rational Interpolation.” *Numerical Algorithms* 88
(1): 365–88.

Lindgren, Finn, Håvard Rue, and Johan Lindström. 2011. “An Explicit Link
Between Gaussian Fields and Gaussian Markov Random Fields: The
Stochastic Partial Differential Equation Approach.” *Journal of the
Royal Statistical Society. Series B. Statistical Methodology* 73 (4):
423–98.

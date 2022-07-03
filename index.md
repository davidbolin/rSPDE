# The `rSPDE` Package <img src="man/figures/logo.png" align="right" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/rSPDE)](https://cran.r-project.org/package=rSPDE)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)

`rSPDE` is an R package used for computing rational approximations of fractional SPDEs. These rational approximations can be used for computatially efficient statistical inference.

Basic statistical operations such as likelihood evaluations and kriging predictions using the fractional approximations are also implemented. The package also contains an interface to [R-INLA][ref4].

# Introduction #

Several popular Gaussian random field models can be represented as solutions to stochastic partial differential equations (SPDEs) of the form 


\\[
L^{\\beta}(\\tau u) = \\mathcal{W}.
\\]

Here \\(\\mathcal{W}\\) is a Gaussian white noise, \\(L\\) is a second-order differential operator, the fractional power \\(\\beta&gt;0\\) determines the smoothness of \\(u\\), and \\(\\tau&gt;0\\) scales the variance of \\(u\\). The simplest example is a model on \\(\\mathbb {R}^d\\) with \\(L = \\kappa^2 - \\Delta\\), which results in a Gaussian random field \\(u\\) with a Matérn covariance function

\\[
C(h) = \\dfrac{ \\sigma^2 }{ 2 ^ {\\nu-1} \\Gamma (\\nu)} (\\kappa h) ^ {\\nu} K_{\\nu} (\\kappa h).
\\]

If \\(2 \\beta\\) is an integer and if the domain \\(\\mathcal {D}\\) where the model is defined is bounded, then \\(u\\) can be approximated by a Gaussian Markov random field (GMRF) \\(\\mathbf { \\mathrm{u} }\\) via a finite element method (FEM) for the SPDE. Specifically, the approximation can be written as

\\[
u_h (s) = \\sum _ { i=1 } ^ n u_i \\varphi_i (s).
\\]

Here \\(\\{\\varphi_i\\}\\) are piecewise linear basis functions defined by some triangulation of \\(\\mathcal {D}\\) and the vector of weights \\( \\mathbf { \\mathrm { u } } = (u_1,\\ldots,u_n)^\\top \\) is normally distributed, \\(N(\\mathbf { \\mathrm{u} }, \\tilde{ \\mathbf { \\mathrm{Q} } }^{-1})\\), where \\(\\tilde{ \\mathbf{ \\mathrm{Q} } }\\) is sparse. See [An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach][ref8] for further details. 

The `rSPDE` package provides corresponding computationally efficient approximations for the case when \\(\\beta\\) is a general fractional power. The main idea is to combine the FEM approximation with a rational approximation of the fractional operator.
As a result, one can easily do inference and prediction using fractional SPDE models such as 

\\[
( \\kappa^2-\\Delta )^\\beta u = \\mathcal{ W }.
\\]

In particular, it allows for bayesian inference of all model parameters, including the fractional parameter \\(\\beta\\).

For illustration purposes, the package contains a simple FEM implementation for models on R. See the 
[An introduction to the rSPDE package][ref2] vignette for an introduction to the package. The [Rational approximation with the rSPDE package][ref6] and [Operator-based rational approximation ][ref5] vignettes provide
introductions to how to create and fit `rSPDE` models. For an introduction to the [`R-INLA`](https://www.r-inla.org) implementation
of the `rSPDE` models see the [R-INLA implementation of the rational SPDE approach][ref3]. The [`rSPDE` documentation][ref7] contains descriptions and examples of the functions in the `rSPDE` package.

# Installation instructions #
The latest CRAN release of the package can be installed directly from CRAN with `install.packages("rSPDE")`.
The latest stable version (which is sometimes slightly more recent than the CRAN version), can be installed by using the command
```r
remotes::install_github("davidbolin/rspde", ref = "stable")
```
in R. The development version can be installed using the command
```r
remotes::install_github("davidbolin/rspde", ref = "devel")
```

If you want to install the package using the `remotes::install_github`-method on Windows, you first need to install `Rtools` and add the paths to `Rtools` and `gcc` to the Windows `PATH` environment variable. This can be done for the current R session only using the commands
```r
rtools = "C:\\Rtools\\bin"
gcc = "C:\\Rtools\\gcc-4.6.3\\bin"
Sys.setenv(PATH = paste(c(gcc, rtools, Sys.getenv("PATH")), collapse = ";"))
```
where the variables `rtools` and `gcc` need to be changed if `Rtool`s is not installed directly on `C:`,
and `gcc`'s version might need to be changed depending on the version of `Rtools`.

# Example

We will illustrate the `rSPDE` package with a kriging example. 

The data consist of precipitation measurements from the Paraná region in Brazil
and were provided by the Brazilian National Water Agency. The data were collected
at 616 gauge stations in Paraná state, south of Brazil, for each day in 2011.
We will not analyse the full spatio-temporal data set, but instead look at the 
total precipitation in January

For further details on the dataset and on the commands we refer the reader
to the [rSPDE-INLA Vignette][ref3].

```r
library(rSPDE)
library(ggplot2)
library(INLA)
library(fields)
library(splancs)
library(gridExtra)
library(lattice)

#Load the data
data(PRprec)
data(PRborder)

#Get the precipitation in January
Y <- rowMeans(PRprec[, 3 + 1:31])

#Treat the data and plot
ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])
alt <- PRprec$Altitude[ind]

ggplot() + geom_point(aes(x = coords[, 1], y = coords[, 2], colour = Y),
                      size = 2,
                      alpha = 1) + 
  scale_colour_gradientn(colours = tim.colors(100)) + 
  geom_path(aes(x = PRborder[, 1], y = PRborder[, 2])) + 
  geom_path(aes(x = PRborder[1034:1078, 1], 
                y = PRborder[1034:1078, 2]), colour = "red")
```

<img src="articles/rspde_inla_files/figure-html/plot_precipitations-1.png">

```r
#Get distance from the sea
seaDist <- apply(spDists(coords, PRborder[1034:1078, ], longlat = TRUE), 1, 
                 min)
                 
#Create the mesh
prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)
plot(prmesh, asp = 1, main = "")
lines(PRborder, col = 3)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")
```

<img src="articles/rspde_inla_files/figure-html/mesh_creation-1.png">

```r
#Create the observation matrix
Abar <- rspde.make.A(mesh = prmesh, loc = coords)

#Create the rspde model object
rspde_model <- rspde.matern(mesh = prmesh)

#Create the index and inla.stack object
mesh.index <- rspde.make.index(name = "field", mesh = prmesh)
stk.dat <- inla.stack(
  data = list(y = Y), A = list(Abar, 1), tag = "est", 
  effects = list(c(mesh.index, 
                   list(Intercept = 1)), 
                 list(long = inla.group(coords[, 1]), 
                      lat = inla.group(coords[,2]),
                      seaDist = inla.group(seaDist))))
                      
#Create the formula object and fit the model
f.s <- y ~ -1 + Intercept +  f(seaDist, model = "rw1") + 
  f(field, model = rspde_model) 
  
rspde_fit <- inla(f.s, family = "Gamma", data = inla.stack.data(stk.dat), 
            verbose = FALSE, 
            control.inla=list(int.strategy='eb'),
            control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE))
            
summary(rspde_fit)
#>
#> 
#> Call:
#>    c("inla(formula = f.s, family = \"Gamma\", data = inla.stack.data(stk.dat), ",
#>" verbose = FALSE, control.predictor = list(A = inla.stack.A(stk.dat), ", " 
#>compute = TRUE), control.inla = list(int.strategy = \"eb\"))" ) 
#> Time used:
#>     Pre = 4.4, Running = 28.1, Post = 0.106, Total = 32.7 
#> Fixed effects:
#>            mean    sd 0.025quant 0.5quant 0.975quant  mode kld
#> Intercept 0.648 0.019      0.611    0.648      0.686 0.648   0
#> 
#> Random effects:
#>   Name	  Model
#>     seaDist RW1 model
#>    field RGeneric2
#> 
#> Model hyperparameters:
#>                                                    mean       sd 0.025quant 0.5quant 0.975quant     mode
#> Precision parameter for the Gamma observations   13.200    0.876     11.548   13.177     14.990   13.139
#> Precision for seaDist                          9378.627 6990.227   2651.630 7372.387  27837.501 5025.442
#> Theta1 for field                                 -1.147    0.311     -1.788   -1.157     -0.521   -1.174
#> Theta2 for field                                  1.089    0.116      0.846    1.090      1.315    1.099
#> Theta3 for field                                 -0.899    0.179     -1.293   -0.898     -0.531   -0.884
#> 
#> Marginal log-Likelihood:  -1261.74 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

#Get the summary on the user's scale
result_fit <- rspde.result(rspde_fit, "field", rspde_model)
summary(result_fit)
#>
#>          mean       sd 0.025quant 0.5quant 0.975quant     mode
#>tau   0.332668 0.106630   0.176185 0.315581   0.595938 0.282519
#>kappa 2.988480 0.344316   2.355250 2.979790   3.711410 2.962190
#>nu    1.161300 0.145758   0.884741 1.160160   1.457410 1.159710

#Plot the posterior densities
par(mfrow=c(1,3))
plot(result_fit)
```

<img src="articles/rspde_inla_files/figure-html/plot_post-1.png">

```r
#Create a grid to predict
nxy <- c(150, 100)
projgrid <- rspde.mesh.projector(prmesh, xlim = range(PRborder[, 1]), 
ylim = range(PRborder[,2]), dims = nxy)
xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[, 1], PRborder[, 2]))
coord.prd <- projgrid$lattice$loc[xy.in, ]

#Compute A matrix and seaDist at predict locations and build the stack
A.prd <- projgrid$proj$A[xy.in, ]
seaDist.prd <- apply(spDists(coord.prd, 
    PRborder[1034:1078, ], longlat = TRUE), 1, min)
ef.prd = list(c(mesh.index, list(Intercept = 1)), 
    list(long = inla.group(coord.prd[, 
    1]), lat = inla.group(coord.prd[, 2]), 
    seaDist = inla.group(seaDist.prd)))
stk.prd <- inla.stack(data = list(y = NA), 
    A = list(A.prd, 1), tag = "prd", 
    effects = ef.prd)
stk.all <- inla.stack(stk.dat, stk.prd)

rspde_fitprd <- inla(f.s, family = "Gamma", 
             data = inla.stack.data(stk.all), 
             control.predictor = list(A = inla.stack.A(stk.all),
             compute = TRUE, link = 1))

id.prd <- inla.stack.index(stk.all, "prd")$data
sd.prd <- m.prd <- matrix(NA, nxy[1], nxy[2])
m.prd[xy.in] <- rspde_fitprd$summary.fitted.values$mean[id.prd]
sd.prd[xy.in] <- rspde_fitprd$summary.fitted.values$sd[id.prd]

#Plot the predictions
grid.arrange(levelplot(m.prd, col.regions = tim.colors(99), 
             xlab = "", ylab = "", main = "mean", 
                       scales = list(draw = FALSE)), 
             levelplot(sd.prd, col.regions = topo.colors(99), 
             xlab = "", ylab = "", scales = list(draw = FALSE), 
                       main = "standard deviation"))
```
<img src="articles/rspde_inla_files/figure-html/plot_pred-1.png">

[ref]: https://www.tandfonline.com/doi/full/10.1080/10618600.2019.1665537  "The rational SPDE approach for Gaussian random fields with general smoothness"
[ref2]: https://davidbolin.github.io/rSPDE//articles/rSPDE.html "An introduction to the rSPDE package"
[ref3]: https://davidbolin.github.io/rSPDE//articles/rspde_inla.html "INLA Vignette"
[ref4]: https://www.r-inla.org "INLA homepage"
[ref5]: https://davidbolin.github.io/rSPDE//articles/rspde_base.html
[ref6]: https://davidbolin.github.io/rSPDE//articles/rspde_cov.html
[ref7]: https://davidbolin.github.io/rSPDE/reference/index.html "`rSPDE` documentation."
[ref8]: https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2011.00777.x "An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach"
################################################################################
## Regenerate the shipped mesh-free weighted-L2 coefficient tables.
##
## The mesh-free covariance tables do not depend on the mesh or on kappa, so
## they are the same for every model and are stored in R/sysdata.rda rather
## than rebuilt at set-up. Building the largest of them takes about four
## seconds; loading it costs the few milliseconds of its quadrature.
##
## Run this whenever anything that changes the coefficients changes: the fit
## itself (wl2_fit, wl2_resjac, wl2_nnls, wl2_lm), the quadrature
## (wl2_weyl_cells), the sweep in wl2_coefficient_table, or the convergence
## tolerance wl2_tol. tests/testthat/test.rational.wl2.R checks that a fresh
## build still reproduces what is stored; if that test starts failing, this
## script is what to run.
##
##     Rscript data-raw/wl2_tables.R
##
## It rewrites R/sysdata.rda in place, preserving everything else in it.
################################################################################

stopifnot(file.exists("DESCRIPTION"), file.exists("R/sysdata.rda"))
pkgload::load_all(".", quiet = TRUE)

## The configurations worth storing: the mesh-free covariance tables, which are
## what matern.operators() uses unless the user asks for a spectral interval.
## The operator tables are not stored: they always come with an x_min, and a
## mesh-free table would only serve as a starting value, which is worth a few
## tenths of a second.
CONFIG <- list(
  d = 1:3,
  m = 1:6,
  m_alpha = 0:2
)

build_all <- function(verbose = TRUE) {
  out <- list()
  for (d in CONFIG$d) {
    for (m_alpha in CONFIG$m_alpha) {
      for (m in CONFIG$m) {
        ## alpha = nu + d / 2 with nu > 0, and the fit needs a finite trace, so
        ## some (d, floor(alpha)) pairs have no valid alpha at all.
        if (length(rSPDE:::wl2_alpha_grid(m_alpha, d, "covariance", 0.01)) == 0) {
          next
        }
        rSPDE:::wl2_clear_cache()
        t0 <- proc.time()[3]
        tab <- rSPDE:::wl2_coefficient_table(
          d = d, m = m, m_alpha = m_alpha, x_min = NULL, type = "covariance",
          cache = FALSE, shipped = FALSE
        )
        dt <- proc.time()[3] - t0
        ## Only the columns are stored; the attributes, the quadrature included,
        ## follow from the configuration and are reattached on lookup.
        attributes(tab) <- attributes(tab)[c("names", "row.names", "class")]
        key <- paste("covariance", d, m, m_alpha, sep = "_")
        out[[key]] <- tab
        if (verbose) {
          cat(sprintf(
            "  %-22s %3d alpha points  %6.2f s  max rel_err %.3e\n",
            key, nrow(tab), dt, max(tab$rel_err)
          ))
          utils::flush.console()
        }
      }
    }
  }
  out
}

cat("Building the mesh-free weighted-L2 tables.\n")
t_all <- proc.time()[3]
wl2_meshfree_tables <- build_all()
cat(sprintf(
  "\n%d tables, %.0f s, %.0f KB\n", length(wl2_meshfree_tables),
  proc.time()[3] - t_all, length(serialize(wl2_meshfree_tables, NULL)) / 1024
))

## Rewrite sysdata.rda, keeping every other object in it.
e <- new.env()
load("R/sysdata.rda", envir = e)
assign("wl2_meshfree_tables", wl2_meshfree_tables, envir = e)
save(
  list = ls(e), envir = e, file = "R/sysdata.rda",
  compress = "xz", version = 2
)
cat(sprintf(
  "wrote R/sysdata.rda (%d objects, %.1f MB)\n",
  length(ls(e)), file.size("R/sysdata.rda") / 1024^2
))

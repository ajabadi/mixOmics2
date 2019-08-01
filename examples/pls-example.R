#' \dontrun{
## successful: FALSE

library(mixOmics.data)

## ---------------- X and Y matrices using linnerud dataset
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls1 <- pls(X=X, Y=Y, mode = "classic")
plotVar(linn.pls)

## ---------------- formula method for numerics
## 'formula' argument should be explicitly mentioned for correct method dispatch
linn.pls2 <- pls(formula = Y ~ X, mode = "classic")
## exclude calls and see if all outputs  are identical
identical(linn.pls1[-1], linn.pls2[-1])
#> TRUE
## ---------------- MultiAssayExperiment using linnerud.mae
## 'data' argument should be explicitly mentioned for correct method dispatch
linn.pls3 <- pls(X='exercise', Y='physiological', mode = "classic", data = linnerud.mae)
identical(linn.pls1[-1], linn.pls3[-1])
#> TRUE

## ---------------- MultiAssayExperiment and formula
linn.pls4 <- pls(formula = physiological ~ exercise, data = linnerud.mae, mode = "classic",)
identical(linn.pls1[-1], linn.pls4[-1])
#> TRUE

## ---------------- colData of MultiAssayExperiment and formula
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.pls <- pls(formula = Y~X, ncomp = 3)

plotVar(toxicity.pls, cutoff = 0.8)
plotLoadings(toxicity.pls, ndisplay = 10)

#' }

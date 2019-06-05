library(mixOmics.data)

X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "classic")

plotVar(linn.pls)

\dontrun{
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.pls <- pls(X, Y, ncomp = 3)
}

## ---- example with formula

linn.pls.form <- pls(formula = Y~X)

## --- example with MultiAssayExperiment:


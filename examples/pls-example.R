X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "classic")

\dontrun{
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.pls <- pls(X, Y, ncomp = 3)
}

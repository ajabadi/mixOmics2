context("pls")


##TODO works
## pls works for matrix-matrix/numeric
## pls fails when xy and formula given - args_conflict
## pls fails when numeric xy and MAE data - args_conflict
## pls fails in formula_mae with non-NULL xy - args_conflict
## pls fails when formula/single X is not a valid assay (symbol or char) from data  - invalid_X - but works otherwise
##
## Do we wanna accept data.frames as well?
## pls works for

library(MultiAssayExperiment)

test_that("pls works as expected",{
  X=assay(miniACC,1) ## RNASeq2GeneNorm
  Y=assay(miniACC,2) ## gistict
  Y_colData = "days_to_death"
  ##---- "xy"
  pls.res.xy <- pls(X =X, Y=Y)
  pls.res.formula <- pls(formula = Y~X)
  pls.res.formula.mae <- pls(formula = gistict ~ RNASeq2GeneNorm, data = miniACC)
  pls.res.xy.mae <- pls(X='gistict' , Y='RNASeq2GeneNorm', data = miniACC)
  ##TODO
  # pls.res.xySym.mae <- pls(X=gistict , Y=RNASeq2GeneNorm, data = miniACC)

  expect_identical(pls.res.xy ,pls.res.formula,pls.res.formula.mae, pls.res.xy.mae )
})
#
# test_that("pls fails when expectattions do not match arguments",{
#
#   ##---- "xy"
#   expect_condition(pls(X =X, Y=Y,formula = Y~X, data = NULL ), class = "args_conflict")
#   expect_condition(pls(X =X, Y=Y,formula = NULL, data = miniACC ), class = "args_conflict")
#   ##---- "formula"
#   expect_condition(pls(X=NULL, Y=NULL,formula = Y~X, data = miniACC ), class = "args_conflict")
#   ##---- 'formula_mae'
#   expect_condition(pls(X=NULL, Y=Y,formula = RNASeq2GeneNorm ~ gistict, data = miniACC ), class = "args_conflict")
#   ##---- 'xy_mae'
#   expect_condition(pls(X=X, Y=Y,formula = RNASeq2GeneNorm ~ gistict, data = miniACC ), class = "args_conflict")
#   expect_condition(pls(X=gistictz, Y=RNASeq2GeneNorm,formula = NULL, data = miniACC ), class = "invalid_X")
#   expect_condition(pls(X=gistict, Y=invalidY,formula = NULL, data = miniACC ), class = "invalid_Y")
# })


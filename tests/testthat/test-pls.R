context("pls")


##TODO

## pls gives inv_xy when Y coldata is not numeric or facor (may be characters should be accepted too?), but works when it is a factor with a warning
## character Y assay with common values for samples is converted to a named numeric and returned in pls result for reference
## ----------------------------------------------------- test data, this is the only section that needs input - you can replicate using new test data

## ------ pls works with numeric ~ matrix
test_that("pls produces identical 'mixo_pls' classes for designated valid signatures when Y is a column data",{
  ## suppress NZV warnings
  suppressMessages({
    pls.res.xy <-          pls(X =Xm_Yc, Y=Ycn )
    pls.res.formula <-     pls(formula =Ycn ~ Xm_Yc)
    pls.res.formula.mae <- pls(formula = f_Yc, data = mae_data)
    pls.res.xy.mae <-      pls(X=X , Y=Yc , data = mae_data)
    })

  expect_identical(pls.res.xy[-1] ,pls.res.formula[-1],pls.res.formula.mae[-1], pls.res.xy.mae[-1] )
  expect_true(class( pls.res.xy)=="mixo_pls")
})

## ------ pls works with matrix ~ matrix
test_that("pls produces identical 'mixo_pls' classes for designated valid signatures when Y is an assay",{
  ## suppress NZV warnings
  suppressMessages({
    pls.res.xy <-          pls(X =Xm_Ya, Y=Yam)
    pls.res.formula <-     pls(formula =Yam ~ Xm_Ya)
    pls.res.formula.mae <- pls(formula = f_Ya, data = mae_data)
    pls.res.xy.mae <-      pls(X=X , Y=Ya , data = mae_data)
  })

  expect_identical(pls.res.xy[-1] ,pls.res.formula[-1],pls.res.formula.mae[-1], pls.res.xy.mae[-1] )
  expect_true(class( pls.res.xy)=="mixo_pls")
})

## ------ correct error with invalid signature combination for X,Y, formula, data
test_that("pls fails with invalid signature and produces appropriate error",{

  expect_error(pls(X=Xm_Ya, Y=Ycn,formula = Y~X ))
  expect_condition(pls(X=Y~Z ), class = "inv_signature")
  expect_condition(pls(X=NULL, Y=Yam,formula = RNASeq2GeneNorm ~ gistict, data = mae_data ), class = "inv_signature")
})

## ------ correct error with invalid assays
test_that("pls fails with invalid assay and produces appropriate error",{

  ##---- "xy"
  expect_condition(pls(formula = Y~X, data = mae_data ), class = "inv_xy")
  expect_condition(pls(X = "invalidX", Y="invalidY", data = mae_data ), class = "inv_xy")

  expect_condition(pls(X=Xm_Ya, Y=Yam,data = mae_data),  class = "inv_xy")
  ##---- "formula"
  expect_condition(pls(formula = Y~X, data = mae_data ), class = "inv_xy")
  ##---- 'formula_mae'
  expect_condition(pls(formula = wrong_LHS ~ gistict, data = mae_data ), class = "inv_xy")
})

## ------ correct error with invalid formula format
test_that("pls fails with invalid formula formats and produces expected errors",{
  expect_condition(pls(formula = Y~X+Y), class = "inv_formula")
  expect_condition(pls(formula = Y+U~X), class = "inv_formula")
})

## ------ correct error with invalid formula elements
test_that("pls fails with invalid formula formats and produces expected errors",{
  expect_condition(pls(formula = Y~U), class = "simpleError")
})

## ------ correct error with non-numeric/factor Y coldata
test_that("pls fails with invalid Y",{
  expect_condition(pls(X=X , Y=Y_inv , data = mae_data), class = "inv_xy")
  expect_warning(pls(X=X , Y=Y_char , data = mae_data), class = "char_Y")
})

## ------ for unnamed and class numeric first two arguments, warning is shown and passed as X and Y for backward
## compatibility
test_that("pls creates a warning for matrix-like formula and data and passes them as X and Y",{
  expect_warning(pls(Xm_Ya, Yam), class = "arg_change")
})

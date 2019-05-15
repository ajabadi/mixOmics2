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
## ------- test data
mae_data <- miniACC
X_index <- 1 ## RNASeq2GeneNorm
Y_index_assay <- 2 ## gistict
Y_index_coldata <- 2 ## "years_to_birth" - must be numeric
## -------

X_assay_name <- names(assays(mae_data))[X_index]
X_assay_mat <- assay(mae_data, X_assay_name)
Y_assay_name <- names(assays(mae_data))[Y_index_assay]
Y_assay_matrix <- assay(mae_data,Y_assay_name)
Y_colData_name <- names(colData(mae_data))[Y_index_coldata]
Y_colData_num <-  colData(mae_data)[,Y_colData_name]

formula_assay<- as.formula(paste0(Y_assay_name, " ~ ", X_assay_name))
formula_coldata <- as.formula(paste0(Y_colData_name, " ~ ", X_assay_name))
## ------

test_that("pls produces identical 'mixo_pls' classes for designated valid signatures",{
  ## suppress NZV warnings
  pls.res.xy <-          suppressWarnings(pls(X =X_assay_mat, Y=Y_assay_matrix))
  pls.res.formula <-     suppressWarnings(pls(formula = Y_assay_matrix~X_assay_mat))
  pls.res.formula.mae <- suppressWarnings(pls(formula = formula_assay, data = mae_data))
  pls.res.xy.mae <-      suppressWarnings(pls(X=X_assay_name , Y=Y_assay_name, data = mae_data))

  expect_identical(pls.res.xy[-1] ,pls.res.formula[-1],pls.res.formula.mae[-1], pls.res.xy.mae[-1] )
  expect_true(class( pls.res.xy)=="mixo_pls")
})
## pls_methods_wrapper
test_that("pls works if Y is an assay or colData",{
  ## suppress NZV warnings
  pls.y.assay <-          suppressWarnings(pls(X =X, Y=Y))
  pls.y.oldata<-     suppressWarnings(pls(formula = Y~X))
  pls.res.formula.mae <- suppressWarnings(pls(formula = gistict ~ RNASeq2GeneNorm, data = mae_data))
  pls.res.xy.mae <-      suppressWarnings(pls(X=X_assay_name , Y=Y_assay_name, data = mae_data))

  expect_identical(pls.res.xy[-1] ,pls.res.formula[-1],pls.res.formula.mae[-1], pls.res.xy.mae[-1] )
  expect_true(class( pls.res.xy)=="mixo_pls")
})

test_that("pls fails with invalid signature and produces appropriate error",{

  expect_condition(pls(X=X, Y=Y,formula = Y~X ), class = "invalid_signature")
  expect_condition(pls(X=Y~Z ), class = "invalid_signature")
})

test_that("pls fails with invalid assay is chosen and produces appropriate error",{

  ##---- "xy"
  expect_condition(pls(formula = Y~X, data = mae_data ), class = "invalid_assay")
  expect_condition(pls(X = "invalidX", Y="invalidY", data = mae_data ), class = "invalid_assay")

  expect_condition(pls(X=X, Y=Y,data = mae_data), class = "invalid_signature")
  ##---- "formula"
  expect_condition(pls(formula = Y~X, data = mae_data ), class = "invalid_assay")
  ##---- 'formula_mae'
  expect_condition(pls(X=NULL, Y=Y,formula = RNASeq2GeneNorm ~ gistict, data = mae_data ), class = "args_conflict")
})

test_that("pls fails with invalid formula formats and produces expected errors",{

  ##---- "xy"
  expect_condition(pls(formula = Y~X+Y), class = "invalid_formula")

})

##TODO
## invalid formula length

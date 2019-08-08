context("pca") ## context helps devtools::test() to categorise tests
test_that('pca results for nutrimouse are the same when using either list or MAE class (with and without quotes OR using assay index for Y)',{
  pca.res1 <- pca(X = nutrimouse$lipid)
  pca.res2 <- pca(X = "lipid", data = nutrimouse.mae)
  ## ignore call slot and expect identical
  expect_identical(pca.res1[-1], pca.res2[-1])
  assayName <- "lipid"
  pca.res3 <- pca(X = assayName , data = nutrimouse.mae)
  expect_identical(pca.res1[-1], pca.res3[-1])
})

##TODO tests for:
## NULL assay

test_that('pca returns error for invalid assay',{
  ## expect error
  expect_error(pca(data = nutrimouse.mae, X = "lipidz"), class = "inv_assay")
})

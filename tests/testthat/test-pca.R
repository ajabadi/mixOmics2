context("pca") ## context helps devtools::test() to categorise tests
test_that('pca results for nutrimouse are the same when using either list or MAE class (with and without quotes OR using assay index for Y)',{
  pca.res1 <- pca(X = nutrimouse$lipid)
  pca.res2 <- pca(data = nutrimouse.mae, X = "lipid")
  ## ignore call slot and expect identical
  expect_identical(pca.res1[-1], pca.res2[-1])
  assay <- 'lipid'
  data <- nutrimouse.mae
  pca.res3 <- pca(data = data, X = assay)
  expect_identical(pca.res1[-1], pca.res3[-1])
})

##TODO tests for:
## NULL assay

test_that('pca returns error for invalid assay',{
  ## expect error
  expect_error(pca(data = nutrimouse.mae, X = "lipidz"), class = "inv_assay")
  expect_error(pca(data = nutrimouse, X = "lipid"), class = "inv_data")
  expect_error(pca(X = "lipid"), class = "inv_data")

})

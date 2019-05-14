context("pca") ## context helps devtools::test() to categorise tests
test_that('pca results for nutrimouse are the same when using either list or MAE class (with and without quotes OR using assay index for Y)',{
  pca.res1 <- pca(X = nutrimouse$lipid)
  pca.res2 <- pca(X = nutrimouse.mae, Assay = "lipid")
  pca.res3 <- pca(X = nutrimouse.mae, Assay = lipid)
  pca.res3 <- pca(X = nutrimouse.mae, Assay = 1)
  ## ignore call slot and expect identical
  expect_identical(pca.res1[-1], pca.res2[-1], pca.res3[-1])
})

##TODO tests for:
## NULL assay

test_that('pca returns error for invalid assay',{
  ## expect error
  expect_error(pca(X = nutrimouse.mae, Assay = "lipidz"), class = "invalid_assay")
  expect_error(pca(X = nutrimouse.mae, Assay = 800), class = "invalid_assay")
})

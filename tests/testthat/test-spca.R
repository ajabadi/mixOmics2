context("spca")

test_that('spca results for nutrimouse are the same when using either list or MAE class (with and without quotes for Y)',{
  spca.res1 <- spca(X = nutrimouse$gene, keepX=c(5,5))
  spca.res2 <- spca(X = nutrimouse.mae, assay = "gene", keepX=c(5,5))
  spca.res3 <- spca(X = nutrimouse.mae, assay = gene, keepX=c(5,5))
  ## ignore call slot and expect identical
  expect_identical(spca.res1[-1], spca.res2[-1], spca.res3[-1])
})

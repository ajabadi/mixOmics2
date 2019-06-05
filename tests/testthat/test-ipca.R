context("ipca")

test_that('ipca results for nutrimouse are the same when using either list or MAE class (with and without quotes for Y)',{
  ipca.res1 <- ipca(X = nutrimouse$gene)
  ipca.res2 <- ipca(X = nutrimouse.mae, assay = "gene")
  ipca.res3 <- ipca(X = nutrimouse.mae, assay = gene)
  ## ignore call slot and expect identical
  expect_identical(ipca.res1[-1], ipca.res2[-1], ipca.res3[-1])
})

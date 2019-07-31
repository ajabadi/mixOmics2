context("sipca")

test_that('sipca results for nutrimouse are the same when using either list or MAE class (with and without quotes for Y)',{
  sipca.res1 <- sipca(liver.toxicity$gene, ncomp = 3, mode="deflation", keepX=c(50,50,50))
  sipca.res2 <- sipca(liver.toxicity.mae, assay='gene', ncomp = 3, mode="deflation", keepX=c(50,50,50))
  ## ignore call slot and expect identical
  expect_identical(sipca.res1[-1], sipca.res2[-1])
})

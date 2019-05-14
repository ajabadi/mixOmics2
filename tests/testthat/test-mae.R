context("test-mae")

# ####################### ipca - calls function from ica.def.par.R
# #### run test
# test_that('ipca results for nutrimouse are the same when using either list or MAE class (with and without quotes for Y)',{
#   ipca.res1 <- ipca(X = liver.toxicity$gene)
#   ipca.res2 <- ipca(X = liver.toxicity.mae, Assay = "gene")
#   ipca.res3 <- ipca(X = liver.toxicity.mae, Assay = gene)
#   ## ignore call slot and expect identical
#   expect_identical(ipca.res1[-1], ipca.res2[-1], ipca.res3[-1])
# })
#
# ####################### sipca - calls function from ica.def.par.R
# #### run test
# test_that('sipca results for nutrimouse are the same when using either list or MAE class (with and without quotes for Y)',{
#   keepX.test <- rep(5,3)
#   sipca.res1 <- sipca(X = liver.toxicity$gene, keepX=keepX.test)
#   sipca.res2 <- sipca(X = liver.toxicity.mae, Assay = "gene", keepX=keepX.test)
#   sipca.res3 <- sipca(X = liver.toxicity.mae, Assay = gene, keepX=keepX.test)
#   ## ignore call slot and expect identical
#   expect_identical(sipca.res1[-1], sipca.res2[-1], sipca.res3[-1])
# })
#
# ####################### plsda
## TODO
##  both Y and formula creates error of class args_conflict
# #### run test
# test_that('plsda results for nutrimouse are the same when using either matrix or MAE class with formula',{
#   plsda.res1 <- plsda(nutrimouse$gene,nutrimouse$diet)
#   plsda.res2 <- plsda(X = nutrimouse.mae, formula=diet~gene)
#   plsda.res3 <- plsda(X=nutrimouse.mae,Assay=gene, phenotype=diet)
#   ## ignore call slot and expect identical
#   expect_identical(plsda.res1[-1], plsda.res1[-1], plsda.res3[-1])
# })

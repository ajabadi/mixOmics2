#' \dontrun{
## successful: TRUE

library(mixOmics.data)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
spca.rat

## variable representation
plotVar(spca.rat, cex=2)

plotVar(spca.rat,style="3d")
# }

# example with MultiAssayExperiment class
# --------------------------------

spca.rat <- spca(liver.toxicity.mae, assay='gene', ncomp = 3, keepX = rep(50, 3))
spca.rat

# \dontrun{


## samples representation
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3],
          group = as.numeric(liver.toxicity$treatment[, 3]))
plotIndiv(spca.rat, cex = 0.01,
          col = as.numeric(liver.toxicity$treatment[, 3]),style="3d")


# example with multilevel decomposition and CLR log ratio transformation
# ----------------

data("diverse.16S")
pca.res = pca(X = diverse.16S$data.TSS, ncomp = 5,
              logratio = 'CLR', multilevel = diverse.16S$sample)
plot(pca.res)
plotIndiv(pca.res, ind.names = FALSE, group = diverse.16S$bodysite, title = '16S diverse data',
          legend=TRUE)

  #' }

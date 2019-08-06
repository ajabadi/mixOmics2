####################################################################################
## ---------- internal
.internal <-

#' Title
#'
#' detailed description

## ----------------------------------- Parameters


## ----------------------------------- Value


## ----------------------------------- Misc
#' @author Florian Rohart, Kim-Anh Lê Cao, Ignacio González, Al J Abadi
#' @seealso \code{\link{nipals}}, \code{\link{prcomp}}, \code{\link{biplot}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}} and http://www.mixOmics.org
#' for more details.
#' @references On log ratio transformations: Filzmoser, P., Hron, K., Reimann,
#' C.: Principal component analysis for compositional data with outliers.
#' Environmetrics 20(6), 621-632 (2009) Lê Cao K.-A., Costello ME, Lakis VA,
#' Bartolo, F,Chua XY, Brazeilles R, Rondeau P. MixMC: Multivariate insights
#' into Microbial Communities. PLoS ONE, 11(8): e0160169 (2016). On multilevel
#' decomposition: Westerhuis, J.A., van Velzen, E.J., Hoefsloot, H.C., Smilde,
#' A.K.: Multivariate paired data analysis: multilevel plsda versus oplsda.
#' Metabolomics 6(1), 119-128 (2010) Liquet, B., Lê Cao, K.-A., Hocini, H.,
#' Thiebaut, R.: A novel approach for biomarker selection and the integration
#' of repeated measures experiments from two assays. BMC bioinformatics 13(1),
#' 325 (2012)
#' @keywords algebra

## ----------------------------------- Examples
#' @example examples/pca-example.R
####################################################################################
## ---------- Generic

#' @foo.. aguments passed to the generic.
#' @usage \S4method{pca}{ANY}(X, ncomp = 2, center = TRUE, scale = FALSE, max.iter = 500,
#' tol = 1e-09, logratio = c('none','CLR','ILR'), ilr.offset = 0.001,
#' V = NULL, multilevel = NULL)
#' @export
setGeneric('pca', function (X, ncomp=2,...) standardGeneric('pca'))

####################################################################################
## ---------- Methods

## ----------------------------------- ANY
#' @export
setMethod('pca', 'ANY', .pca)

## ----------------------------------- MultiAssayExperiment
#' @importFrom SummarizedExperiment assay assays
#' @rdname pca
#' @export
setMethod('pca', 'MultiAssayExperiment', function(X, ncomp=2,..., assay=NULL){

})

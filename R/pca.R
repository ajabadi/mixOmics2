#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Leigh Coonan, Student, University of Quuensland, Australia
#   Fangzhou Yao, Student, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Al J Abadi, Melbourne Integartive Genomics, The University of Melbourne, Australia
#
# created: 2009
# last modified: 2019
#
# Copyright (C) 2009
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################

#' Principal Components Analysis
#'
#' Performs a principal components analysis on the given data matrix that can
#' contain missing values. If data are complete 'pca' uses Singular Value
#' Decomposition, if there are some missing values, it uses the NIPALS
#' algorithm.
#' The calculation is done either by a singular value decomposition of the
#' (possibly centered and scaled) data matrix, if the data is complete or by
#' using the NIPALS algorithm if there is data missing. Unlike
#' \code{\link{princomp}}, the print method for these objects prints the
#' results in a nice format and the \code{plot} method produces a bar plot of
#' the percentage of variance explaned by the principal components (PCs).
#'
#' When using NIPALS (missing values), we make the assumption that the first
#' (\code{min(ncol(X),} \code{nrow(X)}) principal components will account for
#' 100 \% of the explained variance.
#'
#' Note that \code{scale= TRUE} cannot be used if there are zero or constant
#' (for \code{center = TRUE}) variables.
#'
#' Components are omitted if their standard deviations are less than or equal
#' to \code{comp.tol} times the standard deviation of the first component. With
#' the default null setting, no components are omitted. Other settings for
#' \code{comp.tol} could be \code{comp.tol = sqrt(.Machine$double.eps)}, which
#' would omit essentially constant components, or \code{comp.tol = 0}.
#'
#' According to Filzmoser et al., a ILR log ratio transformation is more
#' appropriate for PCA with compositional data. Both CLR and ILR are valid.
#'
#' Logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#'
#' Logratio can only be applied if the data do not contain any 0 value (for
#' count data, we thus advise the normalise raw data with a 1 offset). For ILR
#' transformation and additional offset might be needed.
#'
## --------------------------------------------------------------------------------------- arguments
#'
#' @param X a numeric matrix (or data frame) which provides the data for the
#' principal components analysis. It can contain missing values. Alternatively, a \code{MultiAssayExperiment} object.
#' @param Assay if \code{X} is a \code{MultiAssayExperiment} object, name of an assay from X (string or symbol).
#' @param ncomp integer, if data is complete \code{ncomp} decides the number of
#' components and associated eigenvalues to display from the \code{pcasvd}
#' algorithm and if the data has missing values, \code{ncomp} gives the number
#' of components to keep to perform the reconstitution of the data using the
#' NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X),
#' ncol(X))}
#' @param center a logical value indicating whether the variables should be
#' shifted to be zero centered. Alternately, a vector of length equal the
#' number of columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. The default is
#' \code{FALSE} for consistency with \code{prcomp} function, but in general
#' scaling is advisable. Alternatively, a vector of length equal the number of
#' columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param max.iter integer, the maximum number of iterations in the NIPALS
#' algorithm.
#' @param tol a positive real, the tolerance used in the NIPALS algorithm.
#' @param logratio one of ('none','CLR','ILR'). Specifies the log ratio
#' transformation to deal with compositional values that may arise from
#' specific normalisation in sequencing data. Default to 'none'
#' @param ilr.offset When logratio is set to 'ILR', an offset must be input to
#' avoid infinite value after the logratio transform, default to 0.001.
#' @param V Matrix used in the logratio transformation id provided.
#' @param multilevel sample information for multilevel decomposition for repeated measurements.
#' @param ... arguments passed to methods.
## --------------------------------------------------------------------------------------- value
#' @return \code{pca} returns a list with class \code{"pca"} and
#' \code{"prcomp"} containing the following components: \item{ncomp}{the number
#' of principal components used.} \item{sdev}{the eigenvalues of the
#' covariance/correlation matrix, though the calculation is actually done with
#' the singular values of the data matrix or by using NIPALS.}
#' \item{rotation}{the matrix of variable loadings (i.e., a matrix whose
#' columns contain the eigenvectors).} \item{loadings}{same as 'rotation' to
#' keep the mixOmics spirit} \item{x}{the value of the rotated data (the
#' centred (and scaled if requested) data multiplied by the rotation/loadings
#' matrix), also called the principal components.} \item{variates}{same as 'x'
#' to keep the mixOmics spirit} \item{center, scale}{the centering and scaling
#' used, or \code{FALSE}.} \item{explained_variance}{explained variance from
#' the multivariate model, used for plotIndiv}
## ---------------------------------------------------------------------------------------
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
## --------------------------------------------------------------------------------------- examples
#' @example examples/pca-example.R

#############################################################
## generic function
#############################################################
#' @usage \S4method{pca}{ANY}(X, ncomp = 2, center = TRUE, scale = FALSE, max.iter = 500, tol = 1e-09, logratio = c('none','CLR','ILR'), ilr.offset = 0.001, V = NULL, multilevel = NULL)
## arguemnts must be copied from internal to both @usage and setGeneric plus the '...' in generic so the methods can add arguments - if we only include X, RStudio won't suggest the rest automatically for autofill


## --------------- the importFrom section for automation of NAMESPACE which should ideally be distrbuted to corresponding function files that import them
#' @import MASS lattice igraph ggplot2 corpcor parallel RColorBrewer
#' @importFrom grDevices as.graphicsAnnot chull col2rgb colorRamp colorRampPalette colors dev.cur dev.new dev.off dev.prev dev.set devAskNewPage graphics.off gray gray.colors heat.colors rgb jpeg pdf tiff x11 adjustcolor rainbow
#' @importFrom graphics abline arrows axis barplot box image layout legend lines locator mtext par plot plot.default points polygon rect segments strheight strwidth symbols text title Axis boxplot rasterImage matplot
#' @importFrom stats as.dendrogram as.dist coefficients cor cov dist hclust lm lsfit median na.omit order.dendrogram predict quantile reorder var sd pnorm aggregate t.test
#' @importFrom utils setTxtProgressBar txtProgressBar packageDescription relist download.file
#' @importFrom ellipse ellipse
#' @importFrom methods hasArg is
#' @importFrom dplyr group_by mutate summarise arrange row_number filter n
#' @importFrom tidyr gather
#' @importFrom reshape2 melt dcast
#' @importFrom rARPACK svds
#' @importFrom gridExtra grid.arrange

#' @export
setGeneric("pca", def = function( X, Assay=NULL, ncomp = 2, center = TRUE, scale = FALSE, max.iter = 500, tol = 1e-09,
                                  logratio = 'none', ilr.offset = 0.001, V = NULL, multilevel = NULL,...) standardGeneric("pca"))

#############################################################
## internal function
#############################################################
.pca <- function( X,
                  ncomp = 2,
                  center = TRUE,
                  scale = FALSE,
                  max.iter = 500,
                  tol = 1e-09,
                  logratio = 'none',# one of ('none','CLR','ILR')
                  ilr.offset = 0.001,
                  V = NULL,
                  multilevel = NULL) {
  #-- checking general input parameters --------------------------------------#
  #---------------------------------------------------------------------------#

  #-- check that the user did not enter extra arguments
  arg.call = match.call()
  user.arg = names(arg.call)[-1]

  err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                 error = function(e) e)

  if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)

  #-- X matrix
  if (is.data.frame(X))
    X = as.matrix(X)

  if (!is.matrix(X) || is.character(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)

  if (any(apply(X, 1, is.infinite)))
    stop("infinite values in 'X'.", call. = FALSE)

  #-- put a names on the rows and columns of X --#
  X.names = colnames(X)
  if (is.null(X.names))
    X.names = paste("V", 1:ncol(X), sep = "")

  ind.names = rownames(X)
  if (is.null(ind.names))
    ind.names = 1:nrow(X)

  #-- ncomp
  if (is.null(ncomp))
    ncomp = min(nrow(X),ncol(X))

  ncomp = round(ncomp)

  if ( !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)

  if (ncomp > min(ncol(X), nrow(X)))
    stop("use smaller 'ncomp'", call. = FALSE)

  #-- log.ratio
  choices = c('CLR', 'ILR','none')
  logratio = choices[pmatch(logratio, choices)]

  if (any(is.na(logratio)) || length(logratio) > 1)
    stop("'logratio' should be one of 'CLR' ,'ILR'or 'none'.", call. = FALSE)

  if (logratio != "none" && any(X < 0))
    stop("'X' contains negative values, you can not log-transform your data")


  #-- cheking center and scale
  if (!is.logical(center))
  {
    if (!is.numeric(center) || (length(center) != ncol(X)))
      stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
           call. = FALSE)
  }

  if (!is.logical(scale))
  {
    if (!is.numeric(scale) || (length(scale) != ncol(X)))
      stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
           call. = FALSE)
  }

  #-- max.iter
  if (is.null(max.iter) || !is.numeric(max.iter) || max.iter < 1 || !is.finite(max.iter))
    stop("invalid value for 'max.iter'.", call. = FALSE)

  max.iter = round(max.iter)

  #-- tol
  if (is.null(tol) || !is.numeric(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)

  #-- end checking --#
  #------------------#


  #-----------------------------#
  #-- logratio transformation --#

  if (is.null(V) & logratio == "ILR") # back-transformation to clr-space, will be used later to recalculate loadings etc
    V = clr.backtransfo(X)

  X = logratio.transfo(X = X, logratio = logratio, offset = if(logratio == "ILR") {ilr.offset} else {0})

  #as X may have changed
  if (ncomp > min(ncol(X), nrow(X)))
    stop("use smaller 'ncomp'", call. = FALSE)

  #-- logratio transformation --#
  #-----------------------------#

  #---------------------------------------------------------------------------#
  #-- multilevel approach ----------------------------------------------------#

  if (!is.null(multilevel))
  {
    # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
    multilevel = data.frame(multilevel)

    if ((nrow(X) != nrow(multilevel)))
      stop("unequal number of rows in 'X' and 'multilevel'.")

    if (ncol(multilevel) != 1)
      stop("'multilevel' should have a single column for the repeated measurements.")

    multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements

    Xw = withinVariation(X, design = multilevel)
    X = Xw
  }
  #-- multilevel approach ----------------------------------------------------#
  #---------------------------------------------------------------------------#

  X = scale(X, center = center, scale = scale)
  cen = attr(X, "scaled:center")
  sc = attr(X, "scaled:scale")

  if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance.",
         call. = FALSE)

  is.na.X = is.na(X)
  na.X = FALSE
  if (any(is.na.X)) na.X = TRUE
  NA.X = any(is.na.X)

  cl = match.call()
  cl[[1]] = as.name('pca')
  result = list(call = cl, X = X, ncomp = ncomp,NA.X = NA.X,
                center = if (is.null(cen)) {FALSE} else {cen},
                scale = if (is.null(sc)) {FALSE} else {sc},
                names = list(X = X.names, sample = ind.names))


  #-- pca approach -----------------------------------------------------------#
  #---------------------------------------------------------------------------#

  if(logratio == 'CLR' | logratio=='none')
  {
    #-- if there are missing values use NIPALS agorithm
    if (any(is.na.X))
    {
      res = nipals(X, ncomp = ncomp, reconst = TRUE, max.iter = max.iter, tol = tol)
      result$sdev = res$eig / sqrt(max(1, nrow(X) - 1))
      names(result$sdev) = paste("PC", 1:length(result$sdev), sep = "")
      result$rotation = res$p
      dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
      X[is.na.X] = res$rec[is.na.X]
      result$x = X %*% res$p
      dimnames(result$x) = list(ind.names, paste("PC", 1:ncol(result$x), sep = ""))
    } else {
      #-- if data is complete use singular value decomposition

      #-- borrowed from 'prcomp' function
      res = svd(X, nu = 0)

      result$sdev = res$d[1:ncomp] / sqrt(max(1, nrow(X) - 1))
      result$rotation = res$v[, 1:ncomp, drop = FALSE]
      result$x = X %*% res$v[, 1:ncomp, drop = FALSE]
    }
  } else {
    # if 'ILR', transform data and then back transform in clr space (from RobCompositions package)
    # data have been transformed above
    res = svd(X, nu = max(1, nrow(X) - 1))
    if (ncomp < ncol(X))
    {
      result$sdev = res$d[1:ncomp] / sqrt(max(1, nrow(X) - 1))  # Note: what differs with RobCompo is that they use: cumsum(eigen(cov(X))$values)/sum(eigen(cov(X))$values)
      # calculate loadings using back transformation to clr-space
      result$rotation = V %*% res$v[, 1:ncomp, drop = FALSE]
      # extract component score from the svd, multiply matrix by vector using diag, NB: this differ from our mixOmics PCA calculations
      # NB: this differ also from Filmoser paper, but ok from their code: scores are unchanged
      result$x = res$u[, 1:ncomp, drop = FALSE] %*% diag(res$d[1:ncomp, drop = FALSE])
    } else {
      result$sdev = res$d / sqrt(max(1, nrow(X) - 1))
      result$rotation = V %*% res$v
      result$x = res$u%*% diag(res$d)
    }
  }

  names(result$sdev) = paste("PC", 1:length(result$sdev), sep = "")
  dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
  dimnames(result$x) = list(ind.names, paste("PC", 1:ncol(result$x), sep = ""))

  result$var.tot=sum(X^2 / max(1, nrow(X) - 1))# same as all res$d, or variance after nipals replacement of the missing values

  # to be similar to other methods, add loadings and variates as outputs
  result$loadings = list(X=result$rotation)
  result$variates = list(X=result$x)

  # output multilevel if needed
  if(!is.null(multilevel))
    result=c(result, list(Xw = Xw, design = multilevel))

  class(result) = c("pca","prcomp")
  if(!is.null(multilevel))
    class(result)=c("mlpca",class(result))

  #calcul explained variance
  result$explained_variance = result$sdev^2 / result$var.tot
  result$cum.var = cumsum(result$explained_variance)

  return(invisible(result))
}


#############################################################
## S4 method definitions.
#############################################################

## ---- X: ANY (input handlers will expect matrix or data.frame) ----
#' @export
setMethod("pca", "ANY", function(X, Assay=NULL,...) .pca(X,...))

## ---- X: MultiAssayExperiment ----
#' @rdname pca
#' @importFrom SummarizedExperiment assay
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
setMethod("pca","MultiAssayExperiment", function(X, Assay,...) {
  try_res <- tryCatch(Assay, error = function(e) e)
  if("simpleError" %in% class(try_res)){
    Assay <- as.character(substitute(Assay)) ## internal_mae2dm will check if it is valid
  }
  ## get all inputs so you can refer to provided names
  X <- internal_mae2dm(X = X, Assay = Assay)
  .pca(X=X,...) } )

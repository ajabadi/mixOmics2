#############################################################################################################
# Authors:
#  Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#  Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#  Fangzhou Yao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
#   Al J Abadi, Melbourne Integartive Genomics, The University of Melbourne, Australia

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

#' Independent Principal Component Analysis
#'
#' Performs sparse independent principal component analysis on the given data
#' matrix to enable variable selection.
#'
#' See Details of ipca.
#'
#' Soft thresholding is implemented on the independent loading vectors to
#' obtain sparse loading vectors and enable variable selection.
#'
#'
## --------------------------------------------------------------------------------------- arguments
#'
#' @inheritParams ipca
#' @param keepX numeric vector of length ncomp, the number of variables to keep
## --------------------------------------------------------------------------------------- value
#' @return \code{pca} returns a list with class \code{"ipca"} containing the
#' following components: \item{ncomp}{the number of principal components used.}
#' \item{unmixing}{the unmixing matrix of size (ncomp x ncomp)}
#' \item{mixing}{the mixing matrix of size (ncomp x ncomp} \item{X}{the
#' centered data matrix} \item{x}{the principal components (with sparse
#' independent loadings)} \item{loadings}{the sparse independent loading
#' vectors} \item{kurtosis}{the kurtosis measure of the independent loading
#' vectors}
## ---------------------------------------------------------------------------------------
#' @author Fangzhou Yao and Jeff Coquery.
#' @seealso \code{\link{ipca}}, \code{\link{pca}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}} and http://www.mixOmics.org for more details.
#' @references Yao, F., Coquery, J. and Lê Cao, K.-A. (2011) Principal
#' component analysis with independent loadings: a combination of PCA and ICA.
#' (in preparation)
#'
#' A. Hyvarinen and E. Oja (2000) Independent Component Analysis: Algorithms
#' and Applications, \emph{Neural Networks}, \bold{13(4-5)}:411-430
#'
#' J L Marchini, C Heaton and B D Ripley (2010). fastICA: FastICA Algorithms to
#' perform ICA and Projection Pursuit. R package version 1.1-13.
#' @keywords algebra
## --------------------------------------------------------------------------------------- examples
#' @example examples/sipca-example.R

#############################################################
## S3 generic
#############################################################
#'@rdname sipca
#'@usage {sipca}{default}(X, assay=NULL, ncomp = 3, mode = c("deflation","parallel"), fun = c("logcosh", "exp"), scale = FALSE, max.iter = 200, tol = 1e-04, keepX = rep(50,ncomp), w.init = NULL)
#'@export
sipca <- function(X, assay, ...) UseMethod("sipca")


#############################################################
## internal function
#############################################################
.sipca <-
  function (X, ncomp  = 3, mode = c("deflation","parallel"),
            fun = c("logcosh", "exp"),
            scale = FALSE, max.iter = 200,
            tol = 1e-04, keepX = rep(50,ncomp),
            w.init = NULL)
  {
    dim_x <- dim(X)
    d <- dim_x[dim_x != 1]
    if (length(d) != 2)
      stop("data must be in a matrix form")
    X <- if (length(d) != length(dim_x))
    {matrix(X, d[1], d[2])}
    else {as.matrix(X)}

    alpha <- 1

    mode <- match.arg(mode)
    fun <- match.arg(fun)

    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:ncol(X), sep = "")

    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) ind.names = 1:nrow(X)

    X <- scale(X, scale = FALSE)
    if (scale) {X=scale(X, scale=scale)}
    svd_mat <- svd(X)
    right_sing_vect <- svd_mat$v
    right_sing_vect <- scale(right_sing_vect, center=TRUE, scale=TRUE)
    n <- nrow(t(X))
    p <- ncol(t(X))

    if (ncomp > min(n, p)) {
      message("'ncomp' is too large: reset to ", min(n, p))
      ncomp <- min(n, p)
    }
    if(is.null(w.init))
      w.init <- matrix(1/sqrt(ncomp),ncomp,ncomp)
    else {
      if(!is.matrix(w.init) || length(w.init) != (ncomp^2))
        stop("w.init is not a matrix or is the wrong size")
    }

    X1 <- t(right_sing_vect)[1:ncomp,]

    if (mode == "deflation") {
      unmix_mat <- ica.def(X1, ncomp, tol = tol, fun = fun,
                           alpha = alpha, max.iter = max.iter, verbose = FALSE, w.init = w.init)
    }
    else if (mode == "parallel") {
      unmix_mat <- ica.par(X1, ncomp, tol = tol, fun = fun,
                           alpha = alpha, max.iter = max.iter, verbose = FALSE, w.init = w.init)
    }
    w <- unmix_mat
    independent_mat <- w %*% X1
    #==order independent_mat by kurtosis==#
    kurt <- vector(length=ncomp)
    independent_mat.new <- matrix(nrow = ncomp, ncol = n)
    for(h in 1:ncomp){
      kurt[h] <- (mean(independent_mat[h,]^4)-3*(mean(independent_mat[h,]^2))^2)
    }
    for(i in 1:ncomp){
      independent_mat.new[i,] <- independent_mat[order(kurt,decreasing=TRUE)[i],]
      independent_mat.new[i,] <- independent_mat.new[i,]/as.vector(crossprod(independent_mat.new[i,]))
    }

    #== variable selection==#
    v.sparse=matrix(nrow = ncomp, ncol = n)
    for(i in 1:ncomp){
      nx <- n - keepX[i]
      v.sparse[i,] = ifelse(abs(independent_mat.new[i,]) > abs(independent_mat.new[i,][order(abs(independent_mat.new[i,]))][nx]),
                            (abs(independent_mat.new[i,]) - abs(independent_mat.new[i,][order(abs(independent_mat.new[i,]))][nx])) * sign(independent_mat.new[i,]), 0)
    }
    independent_mat.new = v.sparse


    mix_mat <- t(w) %*% solve(w %*% t(w))

    ipc_mat = matrix(nrow=p, ncol=ncomp)
    ipc_mat = X %*% t(independent_mat.new)
    ##== force orthogonality ==##
    for(h in 1:ncomp){
      if(h==1){ipc_mat[,h]=X %*% (t(independent_mat.new)[,h])}
      if(h>1){ipc_mat[,h]=(lsfit(y=X%*%(t(independent_mat.new)[,h]), ipc_mat[,1:(h-1)],intercept=FALSE)$res)}
      ipc_mat[,h]=ipc_mat[,h]/as.vector(sqrt(crossprod(ipc_mat[,h])))
    }
    ##== force over ==##
    # put rownames of loading vectors
    colnames(independent_mat.new) = colnames(X)

    cl = match.call()
    cl[[1]] = as.name('sipca')

    result = (list(call=cl, X = X, ncomp=ncomp, keepX=keepX, unmixing = t(unmix_mat), mixing = t(mix_mat), loadings = list(X=t(independent_mat.new)), rotation = t(independent_mat.new),
                   kurtosis = kurt[order(kurt,decreasing=TRUE)],names = list(X = X.names, sample = ind.names)))

    result$x = ipc_mat
    result$variates=list(X=ipc_mat)
    dimnames(result$x) = list(ind.names, paste("IPC", 1:ncol(result$rotation), sep = " "))

    class(result) = c("sipca","ipca","pca")

    #calcul explained variance
    explX=explained_variance(X,result$variates$X,ncomp)
    result$explained_variance=explX



    return(invisible(result))
  }

#############################################################
## S3 methods
#############################################################
## --------------------------------------------------------------------------------------- default
#'@export
sipca.default <- function(X, assay=NULL,...) {
  .sipca(X,...)
}

## --------------------------------------------------------------------------------------- MultiAssayExperiment
#'@rdname sipca
#'@importFrom SummarizedExperiment assay
#'@export
sipca.MultiAssayExperiment <-   function(X, assay, ...){
  try_res <- tryCatch(assay, error = function(e) e)
  if("simpleError" %in% class(try_res)){
    assay <- as.character(substitute(assay)) ## internal_mae2dm will check if it is valid
  }
  ## get all inputs so you can refer to provided names
  X <- internal_mae2dm(X = X, Assay = assay)
  .sipca(X=X,...)
}
## --------------------------------------------------------------------------------------- SingleCellExperiment
#'@rdname sipca
#'@export
sipca.SingleCellExperiment <- function(X, assay="logcounts", ...) sipca.MultiAssayExperiment
## --------------------------------------------------------------------------------------- SummarizedExperiment
#'@rdname sipca
#'@export
sipca.SummarizedExperiment <- function(X, assay, ...) sipca.MultiAssayExperiment

###############################################################################
# Authors:
#   Florian Rohart,
#   Benoit Gautier,
#   Kim-Anh Le Cao,
#
# created: 22-04-2015
# last modified: 04-10-2017
#
# Copyright (C) 2015
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
###############################################################################


# ==============================================================================
# block.plsda: perform a horizontal PLS-DA on a combination of datasets,
#   input as a list in X
#   this function is a particular setting of .mintBlock,
#   the formatting of the input is checked in .mintWrapperBlock
# ==============================================================================

# X: a list of data sets (called 'blocks') matching on the same samples.
#   Data in the list should be arranged in samples x variables,
#   with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: a factor or a class vector for the discrete outcome.
# indY: to supply if Y is missing, indicate the position of the outcome in X.
# ncomp: numeric vector of length the number of blocks in \code{X}.
#   The number of components to include in the model for each block
#   (does not necessarily need to take the same value for each block).
#   By default set to 2 per block.
# design: the input design.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid".
#   Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression".
#   Default to "regression"
# scale: boleean. If scale = TRUE, each block is standardized to zero means
#   and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single".
#   Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# all.outputs: calculation of non-essential outputs
#   (e.g. explained variance, loadings.Astar, etc)








#' N-integration with Projection to Latent Structures models (PLS) with
#' Discriminant Analysis
#'
#' Integration of multiple data sets measured on the same samples or
#' observations to classify a discrete outcome, ie. N-integration with
#' Discriminant Analysis. The method is partly based on Generalised Canonical
#' Correlation Analysis.
#'
#' \code{block.plsda} function fits a horizontal integration PLS-DA model with
#' a specified number of components per block). A factor indicating the
#' discrete outcome needs to be provided, either by \code{Y} or by its position
#' \code{indY} in the list of blocks \code{X}.
#'
#' \code{X} can contain missing values. Missing values are handled by being
#' disregarded during the cross product computations in the algorithm
#' \code{block.pls} without having to delete rows with missing data.
#' Alternatively, missing data can be imputed prior using the \code{nipals}
#' function.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and \code{?pls} for more details).
#'
#' Note that our method is partly based on Generalised Canonical Correlation
#' Analysis and differs from the MB-PLS approaches proposed by Kowalski et al.,
#' 1989, J Chemom 3(1) and Westerhuis et al., 1998, J Chemom, 12(5).
#'
#' @param X A list of data sets (called 'blocks') measured on the same samples.
#' Data in the list should be arranged in matrices, samples x variables, with
#' samples order matching in all data sets.
#' @param Y A factor or a class vector indicating the discrete outcome of each
#' sample.
#' @param indY To be supplied if Y is missing, indicates the position of the
#' factor / class vector outcome in the list \code{X}
#' @param ncomp the number of components to include in the model. Default to 2.
#' Applies to all blocks.
#' @param design numeric matrix of size (number of blocks in X) x (number of
#' blocks in X) with values between 0 and 1. Each value indicates the strenght
#' of the relationship to be modelled between two blocks; a value of 0
#' indicates no relationship, 1 is the maximum value. If \code{Y} is provided
#' instead of \code{indY}, the \code{design} matrix is changed to include
#' relationships to \code{Y}.
#' @param scheme Either "horst", "factorial" or "centroid". Default =
#' \code{horst}, see reference.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details. Default = \code{regression}.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances. Default = \code{TRUE}.
#' @param init Mode of initialization use in the algorithm, either by Singular
#' Value Decompostion of the product of each block of X with Y ("svd") or each
#' block independently ("svd.single"). Default = \code{svd}.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default = \code{FALSE}.
#' @param all.outputs boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
#' @return \code{block.plsda} returns an object of class
#' \code{"block.plsda","block.pls"}, a list that contains the following
#' components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates of each block of X.} \item{loadings}{list containing
#' the estimated loadings for the variates.} \item{names}{list containing the
#' names to be used for individuals and variables.} \item{nzv}{list containing
#' the zero- or near-zero predictors information.} \item{iter}{Number of
#' iterations of the algorthm for each component}
#' \item{explained_variance}{Percentage of explained variance for each
#' component and each block}
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.pls}},
#' \code{\link{block.splsda}} and http://www.mixOmics.org for more details.
#' @references On PLSDA:
#'
#' Barker M and Rayens W (2003). Partial least squares for discrimination.
#' \emph{Journal of Chemometrics} \bold{17}(3), 166-173. Perez-Enciso, M. and
#' Tenenhaus, M. (2003). Prediction of clinical outcome with microarray data: a
#' partial least squares discriminant analysis (PLS-DA) approach. \emph{Human
#' Genetics} \bold{112}, 581-592. Nguyen, D. V. and Rocke, D. M. (2002). Tumor
#' classification by partial least squares using microarray gene expression
#' data. \emph{Bioinformatics} \bold{18}, 39-50.
#'
#' On multiple integration with PLS-DA: Gunther O., Shin H., Ng R. T. ,
#' McMaster W. R., McManus B. M. , Keown P. A. , Tebbutt S.J. , Lê Cao K-A. ,
#' (2014) Novel multivariate methods for integration of genomics and proteomics
#' data: Applications in a kidney transplant rejection study, OMICS: A journal
#' of integrative biology, 18(11), 682-95.
#'
#' On multiple integration with sPLS-DA and 4 data blocks:
#'
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO: multi omics integration for biomarker discovery.
#' BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611}
#'
#' mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @examples
#'
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = nutrimouse$diet)
#' # with this design, all blocks are connected
#' design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
#' byrow = TRUE, dimnames = list(names(data), names(data)))
#'
#' res = block.plsda(X = data, indY = 3) # indY indicates where the outcome Y is in the list X
#' plotIndiv(res, ind.names = FALSE, legend = TRUE)
#' plotVar(res)
#'
#' \dontrun{
#' # when Y is provided
#' res2 = block.plsda(list(gene = nutrimouse$gene, lipid = nutrimouse$lipid),
#' Y = nutrimouse$diet, ncomp = 2)
#' plotIndiv(res2)
#' plotVar(res2)
#' }
#'
#' @export block.plsda
block.plsda = function(X,
Y,
indY,
ncomp = 2,
design,
scheme,
mode,
scale = TRUE,
init = "svd",
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
all.outputs = TRUE)

{
    # check inpuy 'Y' and transformation in a dummy matrix
    if (!missing(Y))
    {
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }

        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")

        Y.input = Y
        Y = unmap(Y)
        colnames(Y) = levels(Y.input)
        rownames(Y) = rownames(X[[1]])
    } else if (!missing(indY)) {
        temp = X[[indY]]
        #not called Y to not be an input of the wrapper.sparse.mint.block
        if (is.null(dim(temp)))
        {
            temp = factor(temp)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }

        if (nlevels(temp) == 1)
        stop("'X[[indY]]' should be a factor with more than one level")

        Y.input = temp
        X[[indY]] = unmap(temp)
        colnames(X[[indY]]) = levels(Y.input)
        rownames(X[[indY]]) = rownames(X[[ifelse(indY==1,2,1)]])

    } else if (missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
    }

    # call to '.mintWrapperBlock'
    result = .mintWrapperBlock(X=X, Y=Y, indY=indY, ncomp=ncomp,
        design=design, scheme=scheme, mode=mode, scale=scale,
        init=init, tol=tol, max.iter=max.iter, near.zero.var=near.zero.var,
        all.outputs = all.outputs)

    # calculate weights for each dataset
    weights = .getWeights(result$variates, indY = result$indY)

    # choose the desired output from 'result'
    out=list(call = match.call(),
        X = result$A[-result$indY],
        Y = Y.input,
        ind.mat = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        variates = result$variates,
        loadings = result$loadings,
        crit = result$crit,
        AVE = result$AVE,
        names = result$names,
        init = result$init,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = result$scale,
        design = result$design,
        scheme = result$scheme,
        indY = result$indY,
        weights = weights,
        explained_variance = result$explained_variance)#[-result$indY])

    # give a class
    class(out) = c("block.plsda","block.pls","sgccda","sgcca","DA")
    return(invisible(out))

}




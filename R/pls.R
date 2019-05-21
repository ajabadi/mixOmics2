#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
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


# ========================================================================================================
# pls: perform a PLS
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

#' Partial Least Squares (PLS) Regression
#'
#' Function to perform Partial Least Squares (PLS) regression.
#'
#' \code{pls} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp}
#' components. Multi-response models are fully supported. The \code{X} and
#' \code{Y} datasets can contain missing values.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References). Different modes relate on how the Y matrix is deflated across
#' the iterations of the algorithms - i.e. the different components.
#'
#' - Regression mode: the Y matrix is deflated with respect to the information
#' extracted/modelled from the local regression on X. Here the goal is to
#' predict Y from X (Y and X play an asymmetric role). Consequently the latent
#' variables computed to predict Y from X are different from those computed to
#' predict X from Y.
#'
#' - Canonical mode: the Y matrix is deflated to the information
#' extracted/modelled from the local regression on Y. Here X and Y play a
#' symmetric role and the goal is similar to a Canonical Correlation type of
#' analysis.
#'
#' - Invariant mode: the Y matrix is not deflated
#'
#' - Classic mode: is similar to a regression mode. It gives identical results
#' for the variates and loadings associated to the X data set, but differences
#' for the loadings vectors associated to the Y data set (different
#' normalisations are used). Classic mode is the PLS2 model as defined by
#' Tenenhaus (1998), Chap 9.
#'
#' Note that in all cases the results are the same on the first component as
#' deflation only starts after component 1.
#'
#' The estimation of the missing values can be performed by the reconstitution
#' of the data matrix using the \code{nipals} function. Otherwise, missing
#' values are handled by casewise deletion in the \code{pls} function without
#' having to delete the rows with missing data.
#'
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
## --------------------------------------------------------------------------------------- parameters
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y numeric vector or matrix of responses (for multi-response models).
#' \code{NA}s are allowed.
#' @param ncomp the number of components to include in the model. Default to 2.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Setting this argument to FALSE (when appropriate) will speed up the
#' computations. Default value is FALSE
#' @param logratio one of ('none','CLR'). Default to 'none'
#' @param multilevel Design matrix for repeated measurement analysis, where
#' multlevel decomposition is required. For a one factor decomposition, the
#' repeated measures on each individual, i.e. the individuals ID is input as
#' the first column. For a 2 level factor decomposition then 2nd AND 3rd
#' columns indicate those factors. See examples in \code{?spls}).
#' @param all.outputs boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
## --------------------------------------------------------------------------------------- value
#' @return \code{pls} returns an object of class \code{"pls"}, a list that
#' contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates.} \item{loadings}{list containing the estimated
#' loadings for the \eqn{X} and \eqn{Y} variates.} \item{names}{list containing
#' the names to be used for individuals and variables.} \item{tol}{the
#' tolerance used in the iterative algorithm, used for subsequent S3 methods}
#' \item{iter}{Number of iterations of the algorthm for each component}
#' \item{max.iter}{the maximum number of iterations, used for subsequent S3
#' methods} \item{nzv}{list containing the zero- or near-zero predictors
#' information.} \item{scale}{whether scaling was applied per predictor.}
#' \item{logratio}{whether log ratio transformation for relative proportion
#' data was applied, and if so, which type of transformation.}
#' \item{explained_variance}{amount of variance explained per component (note
#' that contrary to PCA, this amount may not decrease as the aim of the method
#' is not to maximise the variance, but the covariance between data sets).}
#' \item{input.X}{numeric matrix of predictors in X that was input, before any
#' saling / logratio / multilevel transformation.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.}
#' \item{defl.matrix}{residual matrices X for each dimension.}
## ---------------------------------------------------------------------------------------
#' @author Sébastien Déjean and Ignacio González and Kim-Anh Lê Cao.
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}} and
#' http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#'
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#'
#' Abdi H (2010). Partial least squares regression and projection on latent
#' structure regression (PLS Regression). \emph{Wiley Interdisciplinary
#' Reviews: Computational Statistics}, 2(1), 97-106.
#' @keywords regression multivariate
## --------------------------------------------------------------------------------------- examples
#' @example examples/pls-example.R
#'

#############################################################
## generic function
#############################################################
#' @usage \S4method{pls}{ANY}(X, Y, ncomp = 2, scale = TRUE, mode = c("regression", "canonical", "invariant", "classic"), tol = 1e-06, max.iter = 100, near.zero.var = FALSE, logratio = c("none", "CLR"), multilevel = NULL, all.outputs = TRUE)
#' @usage \S4method{pls}{ANY}(formula=Y~X, ncomp = 2, scale = TRUE, mode = c("regression", "canonical", "invariant", "classic"), tol = 1e-06, max.iter = 100, near.zero.var = FALSE, logratio = c("none", "CLR"), multilevel = NULL, all.outputs = TRUE)
## arguemnts must be copied from internal to both @usage and setGeneric plus the '...' in generic so the methods can add arguments - if we only include X, RStudio won't suggest the rest automatically for autofill
#' @export
setGeneric("pls", def = function(formula, data, X, Y, ncomp = 2, scale = TRUE,
                                 mode = c("regression", "canonical", "invariant", "classic"),
                                 tol = 1e-06, max.iter = 100, near.zero.var = FALSE, logratio = "none",
                                 multilevel = NULL, all.outputs = TRUE,...) standardGeneric("pls"))

#############################################################
## internal function
#############################################################
#' @importFrom matrixStats colSds
#' @importFrom matrixStats colVars
.pls = function(X,
Y,
ncomp = 2,
scale = TRUE,
mode = c("regression", "canonical", "invariant", "classic"),
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
logratio = "none",
multilevel = NULL,
all.outputs = TRUE)
{

    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(X = X, Y = Y,ncomp = ncomp, scale = scale, near.zero.var = near.zero.var, mode = mode,
        max.iter = max.iter, tol = tol, logratio = logratio, multilevel = multilevel, DA = FALSE, all.outputs=all.outputs)

    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        variates = result$variates,
        loadings = result$loadings,
        loadings.star = result$loadings.star,
        names = result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        logratio = logratio,
        explained_variance = result$explained_variance,
        input.X = result$input.X,
        mat.c = result$mat.c#,
        #defl.matrix = result$defl.matrix
        )

    class(out) = c("mixo_pls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlpls",class(out))
    }

    return(invisible(out))

}

#############################################################
## S4 method definitions.
#############################################################

# # ## ----------- wrong signature
setMethod("pls", signature("ANY"), definition = function(formula, data,...){
    stop_custom("inv_signature", "incorrect input format to 'spls'. Read ?spls for supported signatures")
})

## ----------- default - both matrices (Y can be numeric) - rest missing
setMethod("pls", signature(formula="missing", data="missing"), definition = function(formula=NULL, data=NULL,...){
    mc <- as.list(match.call()[-1])
    mc
    # if(!all(c("X", "Y") %in% names(mc))) stop_custom("inv_xy", " 'X' and'Y' muast be both provided")mc[c[1,2]] <- mc[c("X","Y")] ## ensure X is the first one
    # mc$X <- eval.parent(mc[[1]])
    # mc$Y <- as.matrix(eval.parent(mc[[2]]))
    # do.call(.pls, mc)
})

## ----------- formula = Y_numeric ~ X_matrix
#' @rdname pls
#' @param formula formula of form \code{dependent_data/group_membership (Y) ~ independent_data (X)}.
#' @export
setMethod("pls", signature(formula="formula", data="missing"), definition = function(formula, data=NULL,...){
    mc <- as.list(match.call()[-1])

    fterms <- as.list(eval.parent(mc$formula))[-1]
    if(any(sapply(as.list(mc$formula), length)!=1)) stop_custom("inv_formula", "formula must be of form: Y~X")
    mc$X <- eval.parent(fterms[[2]])
    mc$Y <- as.matrix(eval.parent(fterms[[1]]))
    mc$formula <- NULL

    do.call(.pls, mc)
})

## ----------- if formula=assay ~ phenotype/assay and data=MAE is provided
#' @rdname pls
#' @param data A \code{MultiAssayExperiment} dataset with at least 2 assays or one assay and one numeric \code{colData}.
#' @export
setMethod("pls", signature(formula="formula", data="MultiAssayExperiment"), definition = function(formula, data,...){
    mc <- as.list(match.call()[-1])

    mc$data <- eval.parent(mc$data)
    mc$formula <- eval.parent(mc$formula)
    fterms <- as.list(mc$formula)[-1]
    ## if the length of each of the formula elements (LHS and RHS) is anything other than 1, stop
    if(any(sapply(fterms, length)!=1)) stop_custom("inv_formula", "formula must be of form: Y~X")
    ## so that formula can also be stored in a variable - makes testing easier as well
    mc[c("Y","X")] <- as.character(fterms)
    mc <- names2mat(mc = mc)
    mc$formula <- mc$data <-  NULL

    do.call(.pls, mc)
})

## ----------- if X=X_assay, Y=Y_assay/Y_colData and data=MAE is provided
#' @export
#' @rdname pls
setMethod("pls", signature(formula="missing", data="MultiAssayExperiment"), definition = function(formula=NULL, data,...){
    mc <- as.list(match.call()[-1])
    mc$data <- eval.parent(data)
    ## ensure there's X and Y in arguments
    if(!all(c("X", "Y") %in% names(mc))) stop_custom("inv_xy", " assay names 'X' and/or 'Y' not provided")
    ## if they are not class character, they can be objects or names
    if((class(try(mc$X)) != "character") | (class(try(mc$Y)) !="character")){
        if(class(try(eval.parent(mc$X)))=="character"){ ## X stored in an abject
            mc$X <- eval.parent(mc$X)
        }
        if(class(try(eval.parent(mc$Y)))=="character"){
            mc$Y <- eval.parent(mc$Y)
        }
        ## in case they are names
        mc$X <- as.character(mc$X)
        mc$Y <- as.character(mc$Y)
    }
    mc <- names2mat(mc = mc)
    mc$data <-  NULL
    do.call(.pls, mc)
})

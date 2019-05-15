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
setGeneric("pls", def = function(X=NULL, Y=NULL, formula=NULL, data=NULL, ncomp = 2, scale = TRUE,
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
## get an expectation for the methods explicit arguments and the arguments, creates numeric X and Y and runs internal pls

pls_methods_helper <- function(Expect, formals,...)({
    formals$Expect <- Expect
    XY.list <- do.call(check_generic_args,args = formals, envir = parent.frame(n=2))
    X <- XY.list$X
    Y <- XY.list$Y
    .pls(X,Y,...)
})


## ----------- default - both matrices (Y can be numeric)
setMethod("pls", signature("ANY"), definition = function(X,Y, formula, data,...){
    stop_custom("args_conflict", "incorrect input format to 'spls'. Read the documenation for valid inputs for each method")
})

## ----------- default - both matrices (Y can be numeric)
setMethod("pls", signature(X="matrix", Y="matrix"), definition = function(X,Y, formula=NULL, data=NULL,...){
    mc <- as.list(match.call(expand.dots = FALSE)[-1])
    do.call(pls_methods_helper, args = list(Expect="xy", formals=mc,...))
})


## ----------- formula = Y_numeric ~ X_matrix
#' @export
#' @rdname pls
setMethod("pls", signature(formula="formula"), definition = function(X=NULL,Y=NULL, formula, data=NULL,...){
    args.default <- c(X=NULL, y=NULL, formula=NULL, data=NULL)
    mc <- match.call(expand.dots = FALSE)[-1]
    mc <- match(c("X","Y",))
    args.g <- args.default
    args.g[names(args.list)] <- args.list
    args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
    do.call(pls_methods_helper, args = list(Expect="formula", formals=args.g,...))
})

## ----------- if formula=assay ~ phenotype/assay and data=MAE is provided
#' @export
#' @rdname pls
setMethod("pls", signature(formula="formula", data="MultiAssayExperiment"), definition = function(X=NULL,Y=NULL, formula, data,...){
    ## making sure there are no missing arguments
    args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
    do.call(pls_methods_helper, args = list(Expect="formula_mae", formals=args.g,...))
})

## ----------- if X=X_assay, Y=Y_assay/Y_colData and data=MAE is provided
#' @export
#' @rdname pls
setMethod("pls", signature(X="character", Y="character", data="MultiAssayExperiment"), definition = function(X,Y, formula=NULL, data,...){
    args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
    do.call(pls_methods_helper, args = list(Expect="xy_mae", formals=args.g,...))
})

#############################################################
## S4 method definitions.
#############################################################

#' ## ----------- default - both matrices (Y can be numeric)
#' setMethod("pls", signature("ANY"), definition = function(X=NULL,Y=NULL, formula=NULL, data=NULL,...){
#'     stop_custom("args_conflict", "incorrect input format to 'spls'. Read the documenation for valid inputs for each method")
#' })
#'
#' ## ----------- default - both matrices (Y can be numeric)
#' setMethod("pls", signature(X="matrix", Y="matrix", formula="NULL", data="NULL"), definition = function(X,Y, formula=NULL, data=NULL,...){
#'
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     do.call(pls_methods_helper, args = list(Expect="xy", formals=args.g,...))
#' })
#'
#'
#' ## ----------- formula = Y_numeric ~ X_matrix
#' #' @export
#' #' @rdname pls
#' setMethod("pls", signature(formula="formula"), definition = function(X=NULL,Y=NULL, formula, data=NULL,...){
#'     args.default <- c(X=NULL, y=NULL, formula=NULL, data=NULL)
#'     mc <- match.call(expand.dots = FALSE)[-1]
#'     mc <- match(c("X","Y",))
#'     args.g <- args.default
#'     args.g[names(args.list)] <- args.list
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     do.call(pls_methods_helper, args = list(Expect="formula", formals=args.g,...))
#' })
#'
#' ## ----------- if formula=assay ~ phenotype/assay and data=MAE is provided
#' #' @export
#' #' @rdname pls
#' setMethod("pls", signature(formula="formula", data="MultiAssayExperiment"), definition = function(X=NULL,Y=NULL, formula, data,...){
#'     ## making sure there are no missing arguments
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     do.call(pls_methods_helper, args = list(Expect="formula_mae", formals=args.g,...))
#' })
#'
#' ## ----------- if X=X_assay, Y=Y_assay/Y_colData and data=MAE is provided
#' #' @export
#' #' @rdname pls
#' setMethod("pls", signature(X="character", Y="character", data="MultiAssayExperiment"), definition = function(X,Y, formula=NULL, data,...){
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     do.call(pls_methods_helper, args = list(Expect="xy_mae", formals=args.g,...))
#' })


#'
#' #############################################################
#' ## S4 method definitions.
#' #############################################################
#' ## ----------- default - both matrices (Y can be numeric)
#' setMethod("pls", signature(X="matrix", Y="matrix"), definition = function(X,Y, formula, data,...){
#'     ## making sure there are no missing arguments
#'     # args.default <- c(X=NULL, y=NULL, formula=NULL, data=NULL)
#'     # args.list <- as.list(match.call(expand.dots = FALSE)[-1])
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     ## add 'Expect' formal for internal - 'xy' since the call matches these two only but we
#'     ## want to make sure other unused arguments are NULL so the user is not confused
#'     args.g$Expect="xy"
#'     args.g
#'     # args.list
#'     XY.list <- do.call(check_generic_args,args = args.g, envir = parent.frame())
#'     XY.list
#'     # X <- XY.list$X
#'     # Y <- XY.list$Y
#'     # .pls(X,Y,...)
#' })
#'
#'
#' ## ----------- formula = Y_numeric ~ X_matrix
#' #' @export
#' #' @rdname pls
#' setMethod("pls", signature(formula="formula", data="NULL"), definition = function(X,Y, formula, data,...){
#'     message("formula")
#'     ## making sure there are no missing arguments
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     ## add 'Expect' formal for internal - 'xy' since the call matches these two only but we
#'     ## want to make sure other unused arguments are NULL so the user is not confused
#'     args.g$Expect="formula"
#'     XY.list <- do.call(check_generic_args,args = args.g, envir = parent.frame())
#'     X <- XY.list$X
#'     Y <- XY.list$Y
#'     .pls(X,Y,...)
#' })
#'
#' ## ----------- if formula=assay ~ phenotype/assay and data=MAE is provided
#' #' @export
#' #' @rdname pls
#' setMethod("pls", signature(formula="formula", data="MultiAssayExperiment"), definition = function(X,Y, formula, data,...){
#'     ## making sure there are no missing arguments
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     ## add 'Expect' formal for internal - 'formula' since the call matches these two only but we
#'     ## want to make sure other unused arguments are NULL so the user is not confused
#'     args.g$Expect="formula_mae"
#'     XY.list <- do.call(check_generic_args,args = args.g, envir = parent.frame())
#'     X <- XY.list$X
#'     Y <- XY.list$Y
#'     .pls(X,Y,...)
#' })
#'
#' ## ----------- if X=X_assay, Y=Y_assay/Y_colData and data=MAE is provided
#' #' @export
#' #' @rdname pls
#' setMethod("pls", signature(formula="NULL", data="MultiAssayExperiment"), definition = function(X,Y, formula, data,...){
#'     ## making sure there are no missing arguments
#'     args.g <- mget(names(formals()),sys.frame(sys.nframe()))[c("X", "Y","formula", "data")]
#'     # args.g[c("X","Y")] <- as.character(substitute(args.g[c("X","Y")]))
#'     ## add 'Expect' formal for internal - 'formula' since the call matches these two only but we
#'     ## want to make sure other unused arguments are NULL so the user is not confused
#'     args.g$Expect="xy_mae"
#'     XY.list <- do.call(check_generic_args,args = args.g, envir = parent.frame())
#'     X <- XY.list$X
#'     Y <- XY.list$Y
#'     .pls(X,Y,...)
#' })

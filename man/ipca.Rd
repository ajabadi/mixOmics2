% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ipca.R
\docType{methods}
\name{ipca}
\alias{ipca}
\alias{ipca,MultiAssayExperiment-method}
\title{Independent Principal Component Analysis}
\usage{
\S4method{ipca}{ANY}X, ncomp = 2,  mode = c("deflation","parallel"),
fun = c("logcosh", "exp", "kur"), scale = FALSE, w.init = NULL,
max.iter = 200, tol = 1e-04)

\S4method{ipca}{MultiAssayExperiment}(X, ncomp = 2, ..., assay = NULL)
}
\arguments{
\item{X}{A numeric matrix (or data frame) which provides the data for the
principal components analysis. It can contain missing values.
Alternatively, a \code{MultiAssayExperiment} object.}

\item{ncomp}{Integer, if data is complete \code{ncomp} decides the number of
components and associated eigenvalues to display from the \code{pcasvd}
algorithm and if the data has missing values, \code{ncomp} gives the number
of components to keep to perform the reconstitution of the data using the
NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X),
ncol(X))}.}

\item{...}{aguments passed to the generic.}

\item{assay}{Name or index of an assay from \code{X}.}

\item{mode}{character string. What type of algorithm to use when estimating
the unmixing matrix, choose one of \code{"deflation"}, \code{"parallel"}.}

\item{fun}{the function used in approximation to neg-entropy in the FastICA
algorithm. Default set to \code{logcosh}, see details of FastICA.}

\item{w.init}{initial un-mixing matrix (unlike FastICA, this matrix is fixed
here).}
}
\value{
\code{ipca} returns a list with class \code{"ipca"} containing the
following components: \item{ncomp}{the number of independent principal
components used.} \item{unmixing}{the unmixing matrix of size (ncomp x
ncomp)} \item{mixing}{the mixing matrix of size (ncomp x ncomp)}
\item{X}{the centered data matrix} \item{x}{the indepenent principal
components} \item{loadings}{the independent loading vectors}
\item{kurtosis}{the kurtosis measure of the independent loading vectors}
}
\description{
Performs independent principal component analysis on the given data matrix,
a combination of Principal Component Analysis and Independent Component
Analysis.

In PCA, the loading vectors indicate the importance of the variables in the
principal components. In large biological data sets, the loading vectors
should only assign large weights to important variables (genes, metabolites
...). That means the distribution of any loading vector should be
super-Gaussian: most of the weights are very close to zero while only a few
have large (absolute) values.

However, due to the existence of noise, the distribution of any loading
vector is distorted and tends toward a Gaussian distribtion according to the
Central Limit Theroem. By maximizing the non-Gaussianity of the loading
vectors using FastICA, we obtain more noiseless loading vectors. We then
project the original data matrix on these noiseless loading vectors, to
obtain independent principal components, which should be also more noiseless
and be able to better cluster the samples according to the biological
treatment (note, IPCA is an unsupervised approach).

\bold{Algorithm} 1. The original data matrix is centered.

2. PCA is used to reduce dimension and generate the loading vectors.

3. ICA (FastICA) is implemented on the loading vectors to generate
independent loading vectors.

4. The centered data matrix is projected on the independent loading vectors
to obtain the independent principal components.
}
\examples{
#' \dontrun{
## successful: TRUE

library(mixOmics.data)

# implement IPCA on a microarray dataset
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")
ipca.res

# samples representation
plotIndiv(ipca.res, ind.names = as.character(liver.toxicity$treatment[, 4]),
          group = as.numeric(as.factor(liver.toxicity$treatment[, 4])))

# example with MultiAssayExperiment class
# --------------------------------

ipca.res <- ipca(liver.toxicity.mae, assay='gene', ncomp = 3, mode="deflation")
ipca.res



  plotIndiv(ipca.res, cex = 1,
            col = as.numeric(as.factor(liver.toxicity$treatment[, 4])),style="3d")

# variables representation with cutoff
plotVar(ipca.res, cex = 1, cutoff = 0.5)


  ## 3d
  plotVar(ipca.res, rad.in = 0.5, cex = 0.5,style="3d", cutoff = 0.8)

  #' }
}
\references{
Yao, F., Coquery, J. and Lê Cao, K.-A. (2011) Principal
component analysis with independent loadings: a combination of PCA and ICA.
(in preparation)

A. Hyvarinen and E. Oja (2000) Independent Component Analysis: Algorithms
and Applications, \emph{Neural Networks}, \bold{13(4-5)}:411-430

J L Marchini, C Heaton and B D Ripley (2010). fastICA: FastICA Algorithms to
perform ICA and Projection Pursuit. R package version 1.1-13.
}
\seealso{
\code{\link{sipca}}, \code{\link{pca}}, \code{\link{plotIndiv}},
\code{\link{plotVar}}, and http://www.mixOmics.org for more details.
}
\author{
Fangzhou Yao and Jeff Coquery.
}
\keyword{algebra}

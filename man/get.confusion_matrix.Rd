% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.confusion_matrix.R
\name{get.confusion_matrix}
\alias{get.confusion_matrix}
\alias{get.BER}
\title{Create confusion table and calculate the Balanced Error Rate}
\usage{
get.confusion_matrix(truth, all.levels, predicted)
}
\arguments{
\item{truth}{A factor vector indicating the true classes of the samples
(typically \code{Y} from the training set).}

\item{all.levels}{Levels of the 'truth' factor. Optional parameter if there
are some missing levels in \code{truth} compared to the fitted predicted
model}

\item{predicted}{Vector of predicted classes (typically the prediction from
the test set). Can contain NA.}

\item{confusion}{result from a \code{get.confusion_matrix} to calculate the
Balanced Error Rate}
}
\value{
\code{get.confusion_matrix} returns a confusion matrix.

\code{get.BER} returns the BER from a confusion matrix
}
\description{
Create confusion table between a vector of true classes and a vector of
predicted classes, calculate the Balanced Error rate
}
\details{
BER is appropriate in case of an unbalanced number of samples per class as
it calculates the average proportion of wrongly classified samples in each
class, weighted by the number of samples in each class. BER is less biased
towards majority classes during the performance assessment.
}
\examples{


# Example
# -----------------------------------

X <- liver.toxicity$gene
Y <- as.factor(liver.toxicity$treatment[, 4])

## if training is perfomed on 4/5th of the original data
samp <- sample(1:5, nrow(X), replace = TRUE)
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
test.predict <- predict(plsda.train, X[test, ], dist = "max.dist")
Prediction <- test.predict$class$max.dist[, 2]

# the confusion table compares the real subtypes with the predicted subtypes for a 2 component model
confusion.mat = get.confusion_matrix(truth = Y[test],
predicted = Prediction)

get.BER(confusion.mat)

}
\references{
mixOmics article:

Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
feature selection and multiple data integration. PLoS Comput Biol 13(11):
e1005752
}
\seealso{
\code{\link{predict}}.
}
\author{
Florian Rohart
}

\name{qtQDA}
\alias{qtQDA}
\alias{discriminant}
\alias{cervical}
\alias{qtQDA.resampling}

\title{Classify RNA-Seq samples}

\description{
Given training set and training labels, classifies the test samples into existing categories.
}

\usage{
qtQDA(training, training.labels, test, prior = NULL,
      dispersion = c("common", "trended", "tagwise"), num.genes = NULL)
qtQDA.resampling(training, training.labels, prior=NULL,
      dispersion=c("common", "trended", "tagwise"), num.genes=NULL,
      resampling=c("bootstrap","cross.validation"), nfold=7, nbs=10)
discriminant(training, training.labels, test, prior = NULL,
      dispersion = c("common", "trended", "tagwise"))
}

\arguments{
\item{training}{a \code{DGEList} object or a matrix containing counts for training samples.}
\item{training.labels}{a character vector of class labels for each sample. 
The labels show which category each sample comes from.}
\item{test}{a \code{DGEList} object, a vector or a matrix containing counts for test sample(s) to be classified.
Number of rows is assumed to be the same as for \code{training}.}
\item{prior}{a numeric vector of prior probabilities for each class. 
Default value is \code{NULL}. If \code{NULL}, prior probabilities will be computed from training samples.}
\item{dispersion}{a string to specify which dispersion to be used. The values can be "common", "trended" or "tagwise". Default value is "tagwise".}
\item{num.genes}{a numeric vector containing number of features to be selected. If more than one value is provided, values will be sorted in descending order and the output will be printed out accordingly. If \code{NULL}, classification will be executed using all of the features.}
\item{resampling}{a character value indicating the resampling method to be used. Values can be "bootstrap", "cross.validation".}
\item{nfold}{a numeric value indicating number of folds to be constructed. Default value is 7.}
\item{nbs}{a numeric value indicating number of bootstrap to be run. Default value is 10.}
}

\details{
This functions are used to classify the \code{test} samples into existing categories using training data.
The data here are assumed to be marginally negative binomial, but dependent.

\code{discriminant} function first estimates the parameters of underlying model from \code{training} data using the statistical methodology implemented in \code{R} package \code{edgeR}, then quantile transforms the data values as proposed in Kochan et al (2019) using those parameter estimations. 
To incorporate the dependence, \code{discriminant} then computes the covariance matrices for each class separately and regularises them with \code{cov.shrink} function implemented in \code{R} package \code{corpcor}. 
Covariance matrix skrinkage approach was developed by Schafer and Strimmer (2005).
Lastly, \code{discriminant} performs the quadratic discriminant analysis performed on the transformed data and estimates the test sample labels.

\code{qtQDA} and \code{qtQDA.resampling} utilise \code{\link{glmLRTqtQDA}} to obtain 	
a data frame where the rows are same as for \code{training} but sorted in descending order by likelihood ratio statistics (LR). Then applying quantile transformation and Gaussian quadratic discriminant analysis using \code{discriminant} function, \code{qtQDA} and \code{qtQDA.resampling} produce the estimated test set labels.

Dispersion parameters are selected from one of the three sophisticated approaches produced from \code{estimateDisp} function in \code{edgeR}.
}

\value{
\code{discriminant} produces a character vector with the test set labels.

\code{qtQDA} produces a list object containing one component:
\item{classes}{a data frame with the test set labels.}

\code{qtQDA.resampling} produces a list object -when \code{resampling} is \code{cross.validation}- containing two components:
\item{class estimations for last fold}{a data frame with the test set labels for the samples selected by last fold.
If \code{resampling} is \code{bootstrap} the name of this list object is 'class estimations for last bootstrap'.}
\item{mean estimated error rates}{a numeric vector  of average error rates estimated for each fold (or bootstrap).}

For both \code{qtQDA} and \code{qtQDA.resampling}, if \code{num.genes} is \code{NULL} test set labels (for \code{qtQDA.resampling} test sets are reselected for each cross validation or bootstrap) are computed by training the algorithm with the entire list of genes; and if \code{num.genes} is of length 2 or more, then a data frame with the column names corresponding to the number of genes provided.
}

\references{
Schafer, J and Strimmer, K (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. 
\emph{Statistical Applications in Genetics and Molecular Biology}, 4(1), pp. 32.
\url{https://doi.org/10.2202/1544-6115.1175}

Chen, Y, Lun, ATL, and Smyth, GK (2014). Differential expression analysis of complex RNA-seq experiments using edgeR. \emph{Statistical Analysis of Next Generation Sequence Data, Somnath Datta and Daniel S. Nettleton (eds), Springer, New York}, pages 51-74. 
\url{http://www.statsci.org/smyth/pubs/edgeRChapterPreprint.pdf}
}

\examples{
data("cervical")

### In this example test sets are included in the training sets to control the algorithm. We expect to see minor or no error in classification.
training <- cervical
training.labels <- c(rep("N", 29), rep("T", 29))
test <- cervical[, c(1:5,54:58)]
test.labels <-c(rep("N", 5), rep("T", 5))
qtQDA(training=training, test=test, training.labels=training.labels, 
num.genes=NULL, prior=NULL, dispersion="tagwise")
qtQDA(training=training, test=test, training.labels=training.labels, 
num.genes=c(20, 50, 100, 200, 300, 500, 714), prior=NULL, dispersion="tagwise")
discriminant(training, training.labels, test, prior=NULL, dispersion="tagwise")
}




















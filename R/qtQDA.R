#' A classification function with qtQDA
#'
#' This function takes  a matrix of genewise counts accross samples, compute posterior probability that it came from each of a list of possible NB-MVN distributions and classifies new samples into existing groups
#' @param training training data where rows represent genes and columns represent samples
#' @param train.labels class labels of samples in the training set
#' @param test test data where rows represent genes and columns represent samples
#' @return numeric vector of posterior probabilities for each group
#' @export
#' @examples
#' y <- matrix(rnbinom(n=1000*200, prob=1/2, size=3),1000,200)
#' colnames(y) <- c(paste0("T",1:100),paste0("N",1:100))
#' y[1:500,101:200] <- y[1:500,101:200]+100
#' training <- y[, trainindex <- sample(c(1:79,110:200),170, replace=FALSE)]
#' test <- y[, testindex <- sample(c(80:109),30, replace=FALSE)]
#' train.labels <- c(rep("T",100), rep("N",100))[trainindex]
#' estimated.class <- qtQDA(training=training, train.labels=train.labels, test=test)
#' estimated.class
#' actual.class <- testindex
#' sum(estimated.class==actual.class)
#' error.rate <- sum(estimated.class==actual.class)/length(actual.class)
#' error.rate

library(corpcor)
library(limma)
library(edgeR)
source("./R/discriminant.R")

qtQDA <- function(training, train.labels, test){
C.QTBC <- list()
covariance.QTBC <- list()
sigma.QTBC <- list()
training <- DGEList(training)
test <- DGEList(test)
type.train <- unique(train.labels)
m <- ncol(test$counts)
prior.prob <- sapply( 1:length(type.train), function(q) {sum(train.labels==type.train[q])/ length(train.labels) } )
X <- model.matrix(~0+factor(train.labels))
training <- estimateDisp(training, X)
fit <- glmQLFit(training, X)
zscoreddata <- edgeR::zscoreNBinom(training$counts, size = 1/training$tagwise.dispersion, mu = fit$fitted.values)
C.QTBC <- lapply( 1:length(type.train), function(q){zscoreddata[ ,grepl( type.train[q], colnames( zscoreddata))] })
covariance.QTBC <- lapply( 1:length(type.train), function(q) { cov(t(C.QTBC [[q]])) } )
sigma.QTBC <- lapply( 1:length(type.train), function(q) { cov.shrink(covariance.QTBC [[q]]) })          
decision.QTBC <- rep(0, m)
discriminant.score <- matrix(0, nrow = nrow(t(test$counts)), ncol = 2)
for(j in 1:m){
        Y.new <- test[,j]$counts
        discriminant.score[j,] <- discriminant(Y.new, prior.prob, fit$coefficients, training$tagwise.dispersion, sigma.QTBC)
}
type <- type.train[apply(discriminant.score, 1, which.max)]
return(type)
}

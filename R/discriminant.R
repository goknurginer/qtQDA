discriminant <- function(training, training.labels, test, prior=NULL, dispersion=c("common", "trended", "tagwise"))
# Given prior probabilities, for a vector of genewise counts,
# compute posterior probability that came from each of a list
# of possible NB-MVN distributions.
# Created 5 November 2018. Last modified 31 July 2019.
{
    # Assign default values
    training <- DGEList(training, group = training.labels, genes = rownames(training))
    test <- DGEList(test)
    training.levels <- unique(training.labels[order(training.labels)])
    nclass <- length(training.levels)
    discriminant.score <- NULL
    LogLik <- rep(0, nclass)
    
    # Check NULL arguments
    if (is.null(prior)) prior <- sapply(1:length(training.levels), function(x) {sum(training.labels == training.levels[x]) / length(training.labels)})
    
    # Estimate dispersion and fit glmQL model
    X <- model.matrix(~0 + factor(training.labels))
    training <- estimateDisp(training, X)
    fit.qtQDA <- glmQLFit(training, X)
    coefficients <- fit.qtQDA$coefficients
    
    # Assign dispersion based on selected option in dispersion argument
    if (dispersion == "common") disp <- training$common.dispersion
    else if (dispersion == "trended") disp <- training$trended.dispersion
    else disp <- training$tagwise.dispersion
    
    # Quantile transform data values as given in the paper
    zscoreddata <- edgeR::zscoreNBinom(training$counts, size=1/disp, mu=fit.qtQDA$fitted.values)
    zscoreddata.by.class <- lapply( 1:nclass, function(q){zscoreddata[, grepl(training.levels[q], colnames(zscoreddata))] })
    
    # Compute covariance matrix for each class
    cov.by.class <- lapply(1:nclass, function(q) {cov(t(zscoreddata.by.class[[q]]))})

    # Reguralize covariance matrix using cov.shrik from corpcor package
    options(warn = -1)
    Sigma <- lapply(1:nclass, function(q) {cov.shrink(cov.by.class[[q]], verbose=FALSE)})
    options(warn = 0)

    # Compute discriminant scores for each test sample and assign them into appropriate class
    for (i in 1:ncol(test)) {
        y <- test[, i]$counts
        for (s in 1:nclass) {
            mean <- exp(log(sum(y)) + coefficients)
            z <- edgeR::zscoreNBinom(y, size=1/disp, mu=mean[,s])
            R <- chol(Sigma[[s]])
            z <- backsolve(R, z, transpose=TRUE)
            LogLik[s] <- -0.5*sum(z^2)
        }
        posterior <- log(prior) + LogLik
        discriminant.score <- rbind(discriminant.score, posterior)
    }
    rownames(discriminant.score) <- colnames(test)
    type <- as.character(training.levels[apply(discriminant.score, 1, which.max)])
    return(type)
}

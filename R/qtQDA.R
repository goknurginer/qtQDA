# classification with qtQDA

classify.qtQDA <- function(training, train.labels, test)
        
        #  training           : training data where rows represent genes and columns represent samples 
        #  train.labels       : class labels of samples in the training set
        #  test               : test data where rows represent genes and columns represent samples
        
        
        #  Necla Kochan and Goknur Giner 
        #  18 May 2019
{
        # Default values  
        C.QTBC <- list()
        covariance.QTBC <- list()
        sigma.QTBC <- list()
        
        training <- DGEList(training)
        test <- DGEList(test)
        
        type.train <- unique(train.labels)
        m <- ncol(test$counts)
        prior.prob <- sapply( 1:length(type.train), function(q) {sum(train.labels==type.train[q])/ length(train.labels) } )
        
        # construct design matrix for proposed classifier
        X <- model.matrix(~0+factor(train.labels))
        training <- estimateDisp(training, X)
        fit <- glmQLFit(training, X)
        
        # ZScoreNBinom transformation for QTBC and QTNBC
        zscoreddata <- edgeR::zscoreNBinom(training$counts, size = 1/training$tagwise.dispersion, mu = fit$fitted.values)
        
        # data sets displaying each subtype 
        C.QTBC <- lapply( 1:length(type.train), function(q){zscoreddata[ ,grepl( type.train[q], colnames( zscoreddata))] })
        
        #covariance matrix for each class
        covariance.QTBC <- lapply( 1:length(type.train), function(q) { cov(t(C.QTBC [[q]])) } )
        
        # shrinking the covariance matrix for each class
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

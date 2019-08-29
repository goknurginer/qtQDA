qtQDA.resampling <- function(training, training.labels, prior=NULL,  dispersion=c("common", "trended", "tagwise"), num.genes=NULL, resampling=c("bootstrap","cross.validation"), nfold=7, nbs=10)
# Given training set and training labels for a test set,
# classify the test set into existing categories.
# Created 9 June 2019. Last modified 22 August 2019.
{

# Assign default values
training <- DGEList(training, group = training.labels, genes = rownames(training))
training.levels <- unique(training.labels[order(training.labels)])
nclass <- length(training.levels)
estimates <- NULL
results <- NULL
mean.er.cv <- NULL
mean.er.bs <- NULL
colnames.results <- NULL

# Check NULL arguments
if(is.null(training.labels)) stop("training.labels can not be NULL") 
if(is.null(prior)) prior <- sapply( 1:length(training.levels), function(x) {sum(training.labels==training.levels[x])/ length(training.labels)})
if(is.null(num.genes)) num.genes <- nrow(training)
num.genes <- num.genes[order(num.genes, decreasing = TRUE)]

# Estimate dispersion, fit glmLRT and save features that are sorted by LR
toptable <- glmLRTqtQDA(training=training, training.labels=training.labels)

# Apply classification for given number of genes
for (i in 1:length(num.genes)){
  message(num.genes[i], " genes selected and classes are being estimated")
  colnames.results <- c(colnames.results, paste0("N", num.genes[i]))
  selected.features <- rownames(toptable[1:num.genes[i],])
  training <- training[rownames(training$counts) %in% selected.features, ]

  if (resampling=="bootstrap") {
    er.bs.all <- NULL
    for (h in 1:nbs){
      # Define indices
      a <- table(training.labels)[1]
      b <- table(training.labels)[2]
      ptest <- 0.3 # proportion of test sets
      index <- c(sample(a), sample((a+1):(a+b)))
      test.index <- c(index[1:round(ptest*a)], index[(a+1):(a+1+round(ptest*b))])
      training.index <- index[-test.index]
      
      # Select training and test sets and convert them into DGEList objects
      training.bs <- DGEList(training$counts[, training.index], group = training$samples[colnames(training$counts[, training.index]),]$group, genes = rownames(training$counts[, training.index]))
      training.labels.bs <- training.bs$samples$group
      training.levels.bs <- unique(training.labels.bs[order(training.labels.bs)])
      test.bs <- DGEList(training$counts[, test.index], group = training$samples[colnames(training$counts[, test.index]),]$group, genes = rownames(training$counts[, test.index]))
      test.labels <- as.character(test.bs$samples$group)    
      
      # Define bs values
      nclass.bs <- ncol(test.bs$counts)
      prior <- sapply( 1:length(training.levels.bs), function(x) {sum(training.labels.bs==training.levels.bs[x])/ length(training.labels.bs)})
      
      # Estimate classes for each test set
      type <- discriminant(training=training.bs, training.labels=training.labels.bs, test=test.bs, prior=prior, dispersion = dispersion)
      
      # Compute number of errors
      er.bs <- sum(!type==test.labels)/nrow(test.bs$samples)
      er.bs.all <- c(er.bs.all, er.bs) 
    }
    results <- cbind(results, type)
    
    # Compute mean number of error across nbs bootstrap
    mean.er.bs <-c(mean.er.bs, mean(er.bs.all))
  }
  else {
    er.cv.all <- NULL
    folds <- cvFolds(NROW(training$samples), K=nfold)
    for (h in 1:nfold) {
      # Select training and test sets and convert them into DGEList objects
      training.cv <- DGEList(training$counts[, folds$subsets[folds$which!=h]], group = training$samples[colnames(training$counts[, folds$subsets[folds$which!=h]]),]$group, genes = rownames(training$counts[, folds$subsets[folds$which!=h]]))
      training.labels.cv <- training.cv$samples$group
      training.levels.cv <- unique(training.labels.cv)
      test.cv <- DGEList(training$counts[, folds$subsets[folds$which == h]], group = training$samples[colnames(training$counts[, folds$subsets[folds$which == h]]),]$group, genes = rownames(training$counts[, folds$subsets[folds$which == h]]))
      test.labels <- as.character(test.cv$samples$group)  
      
      # Define cv values
      nclass.cv <- ncol(test.cv$counts)
      prior <- sapply( 1:length(training.levels.cv), function(x) {sum(training.labels.cv==training.levels.cv[x])/ length(training.labels.cv)})
         
      # Estimate classes for each test set
      type <- discriminant(training=training.cv, training.labels=training.labels.cv,  test=test.cv, prior=prior, dispersion = dispersion)
      er.cv <- sum(!type==test.labels)/nrow(test.cv$samples)
      er.cv.all <- c(er.cv.all, er.cv) 
    }
    results <- cbind(results, type)
    
    # Compute mean number of error across nfold cross validation
    mean.er.cv <-c(mean.er.cv, mean(er.cv.all))
  }
}

if (resampling=="bootstrap") {
  colnames(results) <- colnames.results
  rownames(results) <- test.labels
  estimates$classes <- results
  names(mean.er.bs) <- colnames.results
  estimates$mean.er.bs <- mean.er.bs
  names(estimates) <- c("class estimations for last bootstrap", "mean estimated error rates")
  return(estimates)
}

else {
  colnames(results) <- colnames.results
  rownames(results) <- test.labels
  estimates$classes <- results
  names(mean.er.cv) <- colnames.results
  estimates$mean.er.cv <- mean.er.cv
  names(estimates) <- c("class estimations for last fold", "mean estimated error rates")
  return(estimates)
}  
}

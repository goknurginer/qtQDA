qtQDA <- function(training, training.labels, test, prior=NULL, dispersion=c("common", "trended", "tagwise"), num.genes=NULL)
# Given training set and training labels for a test set,
# classify the test set into existing categories.
# Created 9 June 2019. Last modified 22 August 2019.
{

# Assign default values
training <- DGEList(training, group = training.labels, genes = rownames(training))
test <- DGEList(test)
training.levels <- unique(training.labels[order(training.labels)])
nclass <- ncol(test$counts)
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
  test <- test[rownames(test$counts) %in% selected.features, ]
 
  # Estimate classes for each test set
  type <- discriminant(training=training, training.labels=training.labels, test=test, prior=prior, dispersion = dispersion)
  results <- cbind(results, type)
}
colnames(results) <- colnames.results
rownames(results) <- rownames(test$samples)
estimates$classes <- results
return(estimates)
}

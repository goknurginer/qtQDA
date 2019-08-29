glmLRTqtQDA <- function(training, training.labels)
 # Given training set and training labels
 # outputs the table of differentially expressed genes.
 # Created 10 June 2019. Last modified 30 July 2019.
{
  training <- DGEList(training)
  design <- model.matrix(~training.labels)
  training <- estimateDisp(training, design)
  fit <- glmFit(training, design)
  fit <- glmLRT(fit)
  invisible(topTags(fit, n=Inf)$table)
}

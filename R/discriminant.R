# classification with NBZScore 

discriminant<- function(y, prior.prob, coefficients, dispersion, Sigma)
  #  For a vector of genewise counts, compute posterior probability
  #  that it came from each of a list of possible NB-MVN distributions
  
  #  y          : vector of counts for sample to be classified
  #  prior.prob : vector of prior probabilities for subtypes
  #  mean       : list of NB mean vectors for y, one for each subtype
  #  dispersion : list of NB dispersion vectors for y, one for each subtype
  #  Sigma      : list of correlation matrices, one for each subtype
  
  #  Gordon Smyth
  #  5 Nov 2018
{
  nsubtypes <- length(prior.prob)
  LogLik <- rep(0,nsubtypes)
  for (s in 1:nsubtypes) {
    mean <- exp(log(sum(y)) + coefficients)
    z <- edgeR::zscoreNBinom(y, size=1/dispersion, mu=mean[,s])
    R <- chol(Sigma[[s]])
    z <- backsolve(R,z,transpose=TRUE)
    LogLik[s] <- -0.5*sum(z^2)
  }
  posterior <- log(prior.prob) + LogLik
  posterior
  
}

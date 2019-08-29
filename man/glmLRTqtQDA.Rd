\name{glmLRTqtQDA}
\alias{glmLRTqtQDA}

\title{Fit \code{glmLRT} from \code{edgeR}}
\description{
Given training set and training labels, generates the table of features obtained from \code{glmLRT} in R package \code{edgeR} and sorted by likelihood ratio statistics (LR).
}

\usage{
glmLRTqtQDA(training, training.labels)
}

\arguments{
\item{training}{a \code{DGEList} object or a matrix containing counts for training samples.}
\item{training.labels}{a character vector of class labels for each sample. 
The labels show which category each sample comes from.}
}

\details{
This function estimates the top tags from \code{training} data using the statistical methodologies \code{glmFit} and \code{glmLRT} that implement generalized linear model (glm) methods developed by McCarthy et al (2012) in \code{R} package \code{edgeR}.
}

\value{
a data frame where the rows are same as for \code{training} but sorted in descending order by likelihood ratio statistics (LR).
}

\seealso{
See \code{\link{qtQDA}}  and \code{qtQDA.resampling}.
}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. 
\emph{Nucleic Acids Research}, 40, 4288-4297. 
\url{https://doi.org/10.1093/nar/gks042}
}


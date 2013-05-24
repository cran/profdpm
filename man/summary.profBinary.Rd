\name{summary.profBinary}
\alias{summary.profBinary}
\title{
    Summarize objects of class \code{profBinary}.
}
\description{
    \code{summary.profBinary} is an S3 method to summarize objects of the class \code{profBinary}.
}
\usage{
\method{summary}{profBinary}(object, ...)
}
\arguments{
  \item{object}{an instance of class \code{profBinary}}
  \item{...}{additional arguments (not used)}
}

\details{
    The \code{summary.profBinary} function outputs summary information using the \code{cat} and \code{print} functions. For each unique value of \code{x$clust}, the \code{summary.profBinary} function outputs the number of observation groups assigned to the corresponding cluster. The estimated outcome probabilities and their 95\% credible intervals are also printed for each cluster. The 95\% credible intervals are computed using the marginal posterior distribution, conditional on the estimated data partition. See the package vignette for additional information. 
}
\value{
A list of lists, one for each unique cluster, each with the following elements:
   \item{groups}{The number of observation groups assigned to the corresponding cluster}
   \item{summary}{A data frame containing the estimate and 95\% credible limits for each outcome probability}
}
\author{
Matt Shotwell <matt.shotwell@vanderbilt.edu>
}
\references{
Matthew S. Shotwell (2013). \pkg{profdpm}: An R Package for MAP Estimation in a Class of Conjugate Product Partition Models. \emph{Journal of Statistical Software}, \bold{53(8)}, 1-18. URL http://www.jstatsoft.org/v53/i08/.
}
\seealso{
\code{\link{profBinary}}
}

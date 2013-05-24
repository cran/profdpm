\name{summary.profLinear}
\alias{summary.profLinear}
\title{
    Summarize objects of class \code{profLinear}.
}
\description{
    \code{summary.profLinear} is an S3 method to summarize objects of the class \code{profLinear}.
}
\usage{
\method{summary}{profLinear}(object, ...)
}
\arguments{
  \item{object}{an instance of class \code{profLinear}}
  \item{...}{additional arguments (not used)}
}

\details{
    The \code{summary.profLinear} function outputs summary information using the \code{cat} and \code{print} functions. For each unique value of \code{x$clust}, the \code{summary.profLinear} function outputs the number of observations, and observation groups assigned to the corresponding cluster. The estimated linear coefficients and their approximate 95\% credible intervals are also printed for each cluster. The 95\% credible intervals are computed using a Laplace approximation to the marginal posterior distribution for the linear coefficients, conditional on the estimated data partition. See the package vignette for additional information. 
}
\value{
A list of lists, one for each unique cluster, each with the following elements:
   \item{groups}{The number of observation groups assigned to the corresponding cluster}
   \item{obs}{The number of observations assigned to the corresponding cluster}
   \item{summary}{A data frame containing the estimate and 95\% credible limits for each linear coefficient}
}
\author{
Matt Shotwell <matt.shotwell@vanderbilt.edu>
}
\references{
Matthew S. Shotwell (2013). \pkg{profdpm}: An R Package for MAP Estimation in a Class of Conjugate Product Partition Models. \emph{Journal of Statistical Software}, \bold{53(8)}, 1-18. URL http://www.jstatsoft.org/v53/i08/.
}
\seealso{
\code{\link{profLinear}}
}

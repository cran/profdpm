\name{pci}
\alias{pci}
\title{Partition Comparison Indices}

\description{
  This function computes several partition comparison indices.
}

\usage{
  pci(x1, x2)
}

\arguments{
  \item{x1}{a factor}
  \item{x2}{a factor}
}

\value{
  A named vector with the following elements:
  \item{R}{ Rand index }
  \item{FM}{ Fowlkes and Mallows index }
  \item{W10}{ Wallace 10 index }
  \item{W01}{ Wallace 01 index }
  \item{J}{ Jaccard index }
}

\details{
This function computes indices of similarity between two factors representing the cluster partition of \eqn{n}{n} items. The two vectors must be of the same length. Let \eqn{n_{11}}{n11} be the number of item pairs that occur in the same cluster in both partitions \code{x1} and \code{x2}, \eqn{n_{00}}{n00} the number of item pairs that occur in differenct clusters in both partitions, \eqn{n_{10}}{n10} the number of item pairs that occur in the same cluster in partition \code{x1} but in different clusters in parition \code{x2}, and \eqn{n_{01}}{n01} the number of item pairs that occur in different clusters in partition \code{x1} but in the same cluster in partition \code{x2}. The Rand index is given by \deqn{\frac{n_{11}+n_{00}}{n_{11}+n_{00}+n_{01}+n_{10}}.}{(n11+n00)/(n11+n00+n01+n10).} The Fowlkes and Mallows index is given by \deqn{\frac{n_{11}}{\sqrt{(n_{11}+n_{01})(n_{11}+n_{10})}}.}{n11/sqrt((n11+n10)*(n11+n01)).} The Wallace indices are respectively given by \deqn{\frac{n_{11}}{n_{11}+n_{10}}\qquad\frac{n_{11}}{n_{11}+n_{01}}.}{n11/(n11+n10), and n11/(n11+n01).} The Jaccard index is given by \deqn{\frac{n_{11}}{n_{11}+n_{01}+n_{10}}.}{n11/(n11+n01+n10).}
} 

\references{
Matthew S. Shotwell (2013). \pkg{profdpm}: An R Package for MAP Estimation in a Class of Conjugate Product Partition Models. \emph{Journal of Statistical Software}, \bold{53(8)}, 1-18. URL http://www.jstatsoft.org/v53/i08/.

Rand, W. (1971) Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association 66:846-850

Fowlkes, E. B. and Mallows, C. L. (1983) A Method for Comparing Two Hierarchical Clusterings. Journal of the American Statistical Association 78:553-569

Wallace, D. L. (1983) A Method for Comparing Two Hierarchical Clusterings: comment. Journal of the American Statistical Association 78:569-576
}

\author{Matt Shotwell}

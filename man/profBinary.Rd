\name{profBinary}
\alias{profBinary}
\title{Binary Product Partition Models}

\description{
This function finds the most probable cluster partition in a binary product partition model (\acronym{PPM}). The Dirichlet process mixture of binary models is the default \acronym{PPM}.
}

\usage{
profBinary(formula, data, clust, param, method="agglomerative",
           maxiter=1000, crit=1e-6, verbose=FALSE)
}

\arguments{
\item{formula}{a one-sided formula specifying a set of binary response variables.}
\item{data}{a dataframe where \code{formula} is evaluated}
\item{clust}{optional vector of factors (or coercible to factors) indicating initial clustering among observations.}
\item{param}{optional list containing the any of the named elements \code{alpha}, \code{a0}, and \code{b0} corresponding to the prior parameters of the beta-binary Dirichlet process mixture. The prior parameters of the beta-binary Dirichlet process mixture should all be scalars.}
\item{method}{character string indicating the optimization method to be used. Meaningful values for this string are \code{"stochastic"}, \code{"gibbs"}, \code{"agglomerative"}, and \code{"none"}. 
\itemize{
\item The \code{"stochastic"} method is an iterative stochastic search utilizing `explode' and `merge' operations on the clusters of a partition. At the explode step, a randomly selected subset of observations are redistributed uniformly at random to an existing or new cluster. Each of the exploded observations are then merged with an existing cluster in a sequentially optimal fashion. Optimization involves computing a moving average of the relative change in the marginal posterior distribution over the possible clusters after each iteration. The optimization stopping criterion is the minumim value this quantity can take before stopping the optimization cycle. If the optimization cycle reaches the maximum allowable iterations before meeting the stopping criterion, a warning is issued.

\item The \code{"gibbs"} method implements the Polya urn Gibbs sampler. This method draws samples from the posterior distribution over the cluster partition in a sequential Gibbs fashion. The sample value with greatest posterior mass is returned. See MacEachern (1994) for details.

\item The \code{"agglomerative"} method initially places each observation into seperate clusters. At each iteration, two of the remaining clusters are merged, where the merged clusters are chosen such that the resulting increase in the posterior mass function is maximized. This is repeated until only one cluster remains. The MAP estimate is the cluster partition, among those considered, which maximizes the posterior mass function over the possible cluster partitions. See Ward (1963) for additional details.  

\item The \code{"fast"} method is a modified version of Sequential Update and Greedy Search (SUGS). This SUGS algorithm assigns observations to clusters sequentially. Initally, the first observation forms a singleton cluster. Subsequent observations are assigned to an existing cluster, or form a new singleton cluster by optimizing the associated posterior probabilities, conditional on previous cluster assignments. See Wang and Dunson (2010) for additional details.

\item The \code{"none"} method is typically used in conjunction with the \code{clust} option to specify an initial cluster partition. If \code{"none"} is specified without \code{clust}, a simple algorithm is used to initialize the cluster partition. Otherwise, the cluster partition is initialized using the \code{clust} argument. The posterior statistics are then computed for initialized clusters.
}
}
\item{maxiter}{integer value specifying the maximum number of iterations for the optimization algorithm.}
\item{crit}{numeric scalar constituting a stopping criterion for the \code{"stochastic"} and \code{"gibbs"}  optimization methods.}
\item{verbose}{logical value indicating whether the routine should be verbose in printing.}
}

\value{
An instance of the class \code{profBinary} containing the following objects
\item{y}{the numeric matrix of observations, where rows with missing observations (\code{NA}) are removed}
\item{param}{the list of prior parameters}
\item{clust}{a numeric vector of integers indicating cluster membership for each non-missing observation}
\item{a}{a list of numeric vectors containing the posterior vector \eqn{a}{a} for each cluster}
\item{b}{a list of numeric vectors containing the posterior vector \eqn{b}{b} for each cluster}
\item{logp}{the unnormalized log value of the marginal posterior mass function for the cluster partition evaluated at \code{clust}}
\item{model}{a model frame, resulting from a call to \code{model.frame}}
}

\details{
This function fits a Dirichlet process mixture of binary models (\acronym{DPMBM}) using the profile method. This method will cluster binary observations vectors (rows of \code{y}) into clusters. The cluster partition is estimated by maximizing the marginal posterior distribution over all possible cluster partitions. Each cluster has an associated binary model. The binary model assigns Bernoulli probabilities independently to each binary valued outcome, corresponding to the columns of \code{y}. The prior parameters \code{a0} and \code{b0} assign a beta prior distribution to each outcome probability. Conditional on the estimated cluster partition, each outcome probability is beta distributed \emph{a posteriori}. The function \code{profBinary} returns the associated posterior parameters of the beta destribution for each cluster and outcome probability.  

Missing observations (\code{NA}) are removed automatically and a warning is issued. The return value contains the reduced observation matrix.

}

\seealso{
\code{\link{pci}}
}

\examples{
library(profdpm)

}

\author{Matt Shotwell}
\references{
Ward, J. H. (1963) Heirarchical Grouping to Optimize an Objective Function. Journal of the American Statistical Association 58:236-244 
MacEachern, S. N. (1994) Estimating Normal Means with Conjugate Style Dirichlet Process Prior. Communications in Statistics B 23:727-741 
}


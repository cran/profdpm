\name{profLinear}
\alias{profLinear}
\title{Linear Product Partition Models}

\description{
This function finds the most probable cluster partition in a linear product partition model (\acronym{PPM}). The Dirichlet process mixture of linear models is the default \acronym{PPM}.
}

\usage{
profLinear(y, x, group, clust, param, method="stochastic",
           maxiter=1000, crit=1e-5, verbose=FALSE)
}

\arguments{
\item{y}{numeric vector of observations.}
\item{x}{numeric matrix of covariates.}
\item{group}{optional vector of factors (or coercible to factors) indicating grouping among observations. Observations that are grouped will be clustered together. This is useful if several values form a single longitudinal observation.}
\item{clust}{optional vector of factors (or coercible to factors) indicating initial clustering among observations. Grouped observations (see \code{group}) must have the same \code{clust} value.}
\item{param}{optional list containing the any of the named elements \code{alpha}, \code{a0}, \code{b0}, \code{m0}, and \code{s0} corresponding to the prior parameters of the normal-gamma Dirichlet process mixture. The prior parameters of the normal-gamma Dirichlet process mixture should all be scalars except \code{m0} which should be a vector of length equal to the number of columns in \code{x}.}
\item{method}{character string indicating the optimization method to be used. Meaningful values for this string are \code{"stochastic"}, \code{"gibbs"}, \code{"agglomerative"}, and \code{"none"}. 
\itemize{
\item The \code{"stochastic"} method is still under investigation and may be unsuitable for some datasets. Optimization involves computing a moving average of the relative change in the marginal posterior distribution over the possible clusters after each iteration. The optimization stopping criterion is the minumim value this quantity can take before stopping the optimization cycle. If the optimization cycle reaches the maximum allowable iterations before meeting the stopping criterion, a warning is issued. Even if the optimization criterion is met, the resulting clusters may not represent the global maximum under the marginal posterior distribution over the possible clusters.

\item The \code{"gibbs"} method implements the Polya urn Gibbs sampler. This method draws samples from the posterior distribution over the cluster partition in a sequential Gibbs fashion. The sample value with greatest posterior mass is returned. See MacEachern(1994) for details.

\item The \code{"agglomerative"} method initially places each observation into seperate clusters. At each iteration, two of the remaining clusters are merged, where the merged clusters are chosen such that the resulting increase in the posterior mass function is maximized. This is repeated until only one cluster remains. The MAP estimate is the cluster partition, among those considered, which maximizes the posterior mass function over the possible cluster partitions. See Ward (1963) for additional details.  

\item The \code{"none"} method is typically used in conjunction with the \code{clust} option to specify an initial cluster partition. If \code{"none"} is specified without \code{clust}, a simple algorithm is used to initialize the cluster partition. Otherwise, the cluster partition is initialized using the \code{clust} argument. The posterior statistics are then computed for initialized clusters.
}
}
\item{maxiter}{integer value specifying the maximum number of iterations for the optimization algorithm.}
\item{crit}{numeric scalar constituting a stopping criterion for the \code{"stochastic"} and \code{"gibbs"} optimization methods.}
\item{verbose}{logical value indicating whether the routine should be verbose in printing.}
}

\value{
An instance of the class \code{profLinear} containing the following objects
\item{y}{the numeric vector of observations, where missing observations (\code{NA}) are removed}
\item{x}{the numeric matrix of covariates, where missing covariates (\code{NA}) are removed}
\item{group}{the grouping vector, where missing group values (\code{NA}) are removed}
\item{param}{the list of prior parameters}
\item{clust}{a numeric vector of integers indicating cluster membership for each non-missing observation}
\item{a}{a numeric vector containing the posterior parameter \eqn{a}{a} for each cluster}
\item{b}{a numeric vector containing the posterior parameter \eqn{b}{b} for each cluster}
\item{m}{a list of numeric vectors containing the posterior vector \eqn{m}{m} for each cluster}
\item{s}{a list of numeric matrices containing the posterior matrix \eqn{S}{S} for each cluster}
\item{logp}{the unnormalized log value of the marginal posterior mass function for the cluster partition evaluated at \code{clust}}
}

\details{
This function fits a Dirichlet process mixture of linear models (\acronym{DPMLM}) using the profile method. This method will group the observations into clusters. The clusters are determined by maximizing the marginal posterior distribution over the space of possible clusters. Each cluster has an associated linear model. Notationally, the linear model for cluster \eqn{ k }{ k } has the form \deqn{ y = \gamma^{\prime} x + \epsilon, }{ y = g'x + e, } where \eqn{ y }{ y } and \eqn{ x }{ x } are the observation vector and covariate matrix for a particular cluster, \eqn{ \epsilon }{ e } has a multivariate normal distribution with mean zero and precision matrix \eqn{ \tau I }{ tI }, and \eqn{ \gamma }{ g } is the vector of linear coefficients. In the \acronym{DPMLM}, conditional on the clustering, \eqn{ \gamma }{ g } and \eqn{ \tau }{ t } have a joint normal-gamma posterior distribution of the form \deqn{ p(m, \tau| y, x) = N( \gamma | m, \tau S  ) G( \tau | a, b ), }{ p(m, t| y, x) = N( g | m, tS ) G( t | a, b ), } where \eqn{ N(\cdot) }{ N() } is the multivariate normal density function with mean vector \eqn{ m }{ m } and precision matrix \eqn{ \tau S }{ tS } and \eqn{ G(\cdot) }{ G() } is the gamma density function with shape and scale parameters \eqn{ a }{ a } and \eqn{ b }{ b }. In addition to the cluster indicators, the posterior quantities \eqn{ S }{ S }, \eqn{ m }{ m }, \eqn{ a }{ a }, and \eqn{ b }{ b } are provided for each cluster in the return value.


Missing observations (\code{NA}) are removed automatically and a warning is issued. The return value contains the reduced observation vector and covariate matrix.

}

\seealso{
\code{\link{pci}}
}

\examples{
library(profdpm)

# set up some data
# linear model 0
x0  <- rnorm(50, 0, 3)
y0  <- x0 + rnorm(50, 0, 1) 

# linear model 1
x1  <- rnorm(50, 0, 3)
y1  <- 10 - x1 + rnorm(50, 0, 1)

# add a column of ones to the covariate matrix (intercept)
x   <- cbind(rep(1,100), c(x0, x1))
y   <- c(y0, y1)

# indicate grouping within each linear model
g   <- c( rep(seq(0,4),10), rep(seq(5,9),10) )

# fit the DPMLM
fit <- profLinear(y, x, g)

# plot the resulting fit(s)
plot(x[,2], y)
for(i in 1:length(fit$m)) {
  abline( a=fit$m[[i]][1], b=fit$m[[i]][2] )
}
}
\author{Matt Shotwell}
\references{
Ward, J. H. (1963) Heirarchical Grouping to Optimize an Objective Function. Journal of the American Statistical Association 58:236-244 

MacEachern, S. N. (1994) Estimating Normal Means with Conjugate Style Dirichlet Process Prior. Communications in Statistics B 23:727-741 
}


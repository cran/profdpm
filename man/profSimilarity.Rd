\name{profSimilarity}
\alias{profSimilarity}
\title{Similarity Between Two Instances of the Class profLinear}

\description{
  This function computes the Rand statistic (1971), a metric that quantifies the similarity between two instances of the class \sQuote{profLinear}.
}

\usage{
  profSimilarity(p1, p2)
}

\arguments{
  \item{p1}{an instance of the class \sQuote{profLinear}}
  \item{p2}{an instance of the class \sQuote{profLinear}}
}

\value{
  A real number between 0.0 and 1.0, where 0.0 indicates no similarity and 1.0 indicates identical clustering.
}

\details{
  This function computes the similarity between two instances of the \sQuote{profLinear} class. The two instances should be the result of two independant calls to the function \sQuote{profLinear} with the same dataset. If observation \eqn{ latex }{ k } is clustered together with observation \eqn{ latex }{ j } in both \sQuote{p1} and \sQuote{p2}, then this constitutes an equivalent clustering. Hence, there are \eqn{ latex }{ n(n-1)/2 } possible equivalent clusterings, where \eqn{ latex }{ n } is the number of clustered observations in \sQuote{p1} and \sQuote{p2}. The metric is computed as the number of equivalent clusterings divided by the number possible equivalent clusterings. See Rand (1971) for additional details.
}

\references{
Rand, W. (1971) Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association 66:846-850
}

\author{Matt Shotwell}

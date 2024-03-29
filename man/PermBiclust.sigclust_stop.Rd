% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SC-Biclust-SigClust-stop.R
\name{PermBiclust.sigclust_stop}
\alias{PermBiclust.sigclust_stop}
\title{'SCBiclust' method for identifying means-based biclusters with optional cluster significance testing}
\usage{
PermBiclust.sigclust_stop(
  x,
  nperms = 1000,
  silent = TRUE,
  maxnum.bicluster = 5,
  alpha = 0.05,
  icovest = 1,
  sc = TRUE
)
}
\arguments{
\item{x}{a dataset with n rows and p columns, with observations in rows.}

\item{nperms}{number of \eqn{Beta(\frac{1}{2}, (p-1)/2)} distributed variables generated for each feature (default=1000)}

\item{silent}{should progress be printed? (default=TRUE)}

\item{maxnum.bicluster}{The maximum number of biclusters returned}

\item{alpha}{significance level for \code{\link{sigclust}} test.}

\item{icovest}{Coviariance estimation type for \code{\link{sigclust}} test}

\item{sc}{should the \code{\link{sigclust}} test be used? (default=TRUE)}
}
\value{
The function returns a S3-object with the following attributes:
\itemize{
\item{\code{num.bicluster}}: {The number of biclusters estimated by the procedure.}
\item{\code{x.residual}}: {The data matrix \code{x} after removing the signals}
\item{\code{which.x}}: A list of length \code{num.bicluster} with each list entry containing a 
logical vector denoting if the data observation is in the given bicluster. 
\item{\code{which.y}}: A list of length \code{num.bicluster} with each list entry containing a 
logical vector denoting if the data feature is in the given bicluster. 
}
}
\description{
'SCBiclust' method for identifying means-based biclusters with optional cluster significance testing
}
\details{
Observations in the bicluster are identified such that they maximize the feature-weighted between cluster sum of squares.
Features in the bicluster are identified based on their contribution to the clustering of the observations. 
Feature weights are generated in a similar fashion as  \code{\link{KMeansSparseCluster}} 
except with a modified objective function and no sparsity constraint.
This algoritm uses a numerical approximation to  \eqn{E(\sqrt{B})} where \eqn{B \sim Beta(\frac{1}{2}, (p-1)/2)} as the expected null 
distribution for feature weights. Use of the \code{\link{sigclust}} algorithm  to test the strength of the identified clusters is optional in this
implementation of the algorithm.
}
\examples{
test <- matrix(rnorm(60*180), nrow=60, ncol=180)
test[1:15,1:15] <- test[1:15,1:15]+rnorm(15*15, 2)
test[16:30,51:80] <- test[16:30,51:80]+rnorm(15*30, 3)
PermBiclust.sigclust_stop(test, silent=TRUE)
}
\author{
Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair
}

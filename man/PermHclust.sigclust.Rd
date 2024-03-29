% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SC-HBiclust-SigClust.R
\name{PermHclust.sigclust}
\alias{PermHclust.sigclust}
\title{'SCBiclust' method for identifying hierarchically clustered biclusters}
\usage{
PermHclust.sigclust(
  x = NULL,
  method = c("average", "complete", "single", "centroid"),
  wbound = sqrt(ncol(x)),
  alpha = 0.05,
  dat.perms = 1000,
  dissimilarity = c("squared.distance", "absolute.value"),
  silent = TRUE,
  sigstep = FALSE
)
}
\arguments{
\item{x}{a dataset with n rows and p columns, with observations in rows.}

\item{method}{method for agglomeration. See documentation in \code{\link{hclust}}. (default="average")}

\item{wbound}{the tuning parameter for sparse hierarchical clustering. See documentation in \code{\link{HierarchicalSparseCluster}}. (default=\code{sqrt(ncol(x))})}

\item{alpha}{significance level for \code{\link{sigclust}} test.}

\item{dat.perms}{number of \eqn{Beta(\frac{1}{2}, (p-1)/2)} distributed variables generated for each feature (default=1000)}

\item{dissimilarity}{How should dissimilarity be calculated? (default is "squared.distance").}

\item{silent}{should progress be printed? (default=TRUE)}

\item{sigstep}{Should \code{\link{sigclust}} be used to assess the strength of identified clusters? (default=FALSE)}
}
\value{
The function returns a S3-object with the following attributes:
\itemize{
\item{\code{which.x}}: A list of length \code{num.bicluster} with each list entry containing a 
logical vector denoting if the data observation is in the given bicluster. 
\item{\code{which.y}}: A list of length \code{num.bicluster} with each list entry containing a 
logical vector denoting if the data feature is in the given bicluster. 
}
}
\description{
'SCBiclust' method for identifying hierarchically clustered biclusters
}
\details{
Observations in the bicluster are identified such that they maximize the feature-weighted version of the dissimilarity matrix as implemented in
\code{\link{HierarchicalSparseCluster}}.
Features in the bicluster are identified based on their contribution to the clustering of the observations. 
#' This algoritm uses a numerical approximation to  \eqn{E(\sqrt{B})} where \eqn{B \sim Beta(\frac{1}{2}, (p-1)/2)} as the expected null 
distribution for feature weights.
}
\examples{
test <- matrix(nrow=500, ncol=50)
theta <- rep(NA, 500)
theta[1:300] <- runif(300, 0, pi)
theta[301:500] <- runif(200, pi, 2*pi)
test[1:300,seq(from=2,to=40,by=2)] <- -2+5*sin(theta[1:300])
test[301:500,seq(from=2,to=40,by=2)] <- 5*sin(theta[301:500])
test[1:300,seq(from=1,to=39,by=2)] <- 5+5*cos(theta[1:300])
test[301:500,seq(from=1,to=39,by=2)] <- 5*cos(theta[301:500])
test[,1:40] <- test[,1:40] + rnorm(40*500, 0, 0.2)
test[,41:50] <- rnorm(10*500, 0, 1)
test.PermBiclust <- PermHclust.sigclust(x=test, method='single', dissimilarity='squared.distance')
}
\author{
Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair
}

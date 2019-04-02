#' 'SCBiclust' method for identifying hierarchically clustered biclusters
#'
#' @param x a dataset with n rows and p columns, with observations in rows. 
#' @param method method for agglomeration. See documentation in \code{\link{hclust}}. (default="average")
#' @param wbound the tuning parameter for sparse hierarchical clustering. See documentation in \code{\link{HierarchicalSparseCluster}}. (default=\code{sqrt(ncol(x))})
#' @param alpha significance level for \code{\link{sigclust}} test.
#' @param dat.perms number of \eqn{Beta(\frac{1}{2}, (p-1)/2)} distributed variables generated for each feature (default=1000)
#' @param dissimilarity How should dissimilarity be calculated? (default is "squared.distance").
#' @param silent should progress be printed? (default=TRUE)
#' @param sigstep Should \code{\link{sigclust}} be used to assess the strength of identified clusters? (default=FALSE)
#'
#' @return
#' The function returns a S3-object with the following attributes:
#' \itemize{
#' \item{\code{which.x}}: A list of length \code{num.bicluster} with each list entry containing a 
#' logical vector denoting if the data observation is in the given bicluster. 
#' \item{\code{which.y}}: A list of length \code{num.bicluster} with each list entry containing a 
#' logical vector denoting if the data feature is in the given bicluster. 
#' }
#' 
#' @details
#' Observations in the bicluster are identified such that they maximize the feature-weighted version of the dissimilarity matrix as implemented in
#' \code{\link{HierarchicalSparseCluster}}.
#' Features in the bicluster are identified based on their contribution to the clustering of the observations. 
#' #' This algoritm uses a numerical approximation to  \eqn{E(\sqrt{B})} where \eqn{B \sim Beta(\frac{1}{2}, (p-1)/2)} as the expected null 
#' distribution for feature weights. 
#' @examples
#' test <- matrix(nrow=500, ncol=50)
#' theta <- rep(NA, 500)
#' theta[1:300] <- runif(300, 0, pi)
#' theta[301:500] <- runif(200, pi, 2*pi)
#' test[1:300,seq(from=2,to=40,by=2)] <- -2+5*sin(theta[1:300])
#' test[301:500,seq(from=2,to=40,by=2)] <- 5*sin(theta[301:500])
#' test[1:300,seq(from=1,to=39,by=2)] <- 5+5*cos(theta[1:300])
#' test[301:500,seq(from=1,to=39,by=2)] <- 5*cos(theta[301:500])
#' test[,1:40] <- test[,1:40] + rnorm(40*500, 0, 0.2)
#' test[,41:50] <- rnorm(10*500, 0, 1)
#' test.PermBiclust <- PermHclust.sigclust(x=test, method='single', dissimilarity='squared.distance')
#' @export
#' @name PermHclust.sigclust
#' @author Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair



PermHclust.sigclust <- function(x=NULL, method=c('average', 'complete', 'single', 'centroid'), 
wbound=sqrt(ncol(x)), alpha=0.05, dat.perms=1000, dissimilarity=c('squared.distance', 'absolute.value'),silent=TRUE,sigstep=FALSE)
{

	x <- scale(x)
	sparsehc <- sparcl::HierarchicalSparseCluster(x, silent=TRUE,wbound=wbound, dissimilarity=dissimilarity, method=method)
	
	hcc <- cutree(sparsehc$hc, 2)
	hcc.perm <- hcc

	ws.perm <- matrix(NA, nrow=dat.perms, ncol=ncol(x))
	for(k in 1:dat.perms){
		ws.perm[k,] <- sqrt(rbeta(ncol(x), .5, (ncol(x)-1)/2))
		ws.perm[k,] <- sort(ws.perm[k,])
	}
	sws <- sort(sparsehc$ws)
	 spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])
	data<-x[,which(sparsehc$ws>sws[sig.ndx])]


if(min(table(hcc))<2|sum(sparsehc$ws>sws[sig.ndx])<2){
which.y<-rep(NA,ncol(x))
which.x<-rep(NA,nrow(x))}else if (sigstep==TRUE){
sc.pval <- sigclust::sigclust(data,nsim=1000,labflag=1, label=hcc)@pval 
	if(sc.pval >= alpha){
		message('Warning: no bicluster found \n')
		which.y<-rep(NA,ncol(x))
		which.x<-rep(NA,nrow(x))} else{
which.y <- (sparsehc$ws>sws[sig.ndx])
which.x<-hcc}}else{
which.y <- (sparsehc$ws>sws[sig.ndx])
which.x<-hcc
}
return(list(which.x=which.x, which.y=which.y))

}
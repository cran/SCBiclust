#' 'SCBiclust' method for identifying means-based biclusters with optional cluster significance testing
#'
#' @param x a dataset with n rows and p columns, with observations in rows. 
#' @param nperms number of \eqn{Beta(\frac{1}{2}, (p-1)/2)} distributed variables generated for each feature (default=1000)
#' @param silent should progress be printed? (default=TRUE)
#' @param maxnum.bicluster The maximum number of biclusters returned  
#' @param alpha significance level for \code{\link{sigclust}} test.
#' @param icovest Coviariance estimation type for \code{\link{sigclust}} test
#' @param sc should the \code{\link{sigclust}} test be used? (default=TRUE)
#'

#' @return
#' The function returns a S3-object with the following attributes:
#' \itemize{
#' \item{\code{num.bicluster}}: {The number of biclusters estimated by the procedure.}
#' \item{\code{x.residual}}: {The data matrix \code{x} after removing the signals}
#' \item{\code{which.x}}: A list of length \code{num.bicluster} with each list entry containing a 
#' logical vector denoting if the data observation is in the given bicluster. 
#' \item{\code{which.y}}: A list of length \code{num.bicluster} with each list entry containing a 
#' logical vector denoting if the data feature is in the given bicluster. 
#' }
#' 
#' @details
#' Observations in the bicluster are identified such that they maximize the feature-weighted between cluster sum of squares.
#' Features in the bicluster are identified based on their contribution to the clustering of the observations. 
#' Feature weights are generated in a similar fashion as  \code{\link{KMeansSparseCluster}} 
#' except with a modified objective function and no sparsity constraint.
#' This algoritm uses a numerical approximation to  \eqn{E(\sqrt{B})} where \eqn{B \sim Beta(\frac{1}{2}, (p-1)/2)} as the expected null 
#' distribution for feature weights. Use of the \code{\link{sigclust}} algorithm  to test the strength of the identified clusters is optional in this
#' implementation of the algorithm.
#'  
#' @examples
#' test <- matrix(rnorm(60*180), nrow=60, ncol=180)
#' test[1:15,1:15] <- test[1:15,1:15]+rnorm(15*15, 2)
#' test[16:30,51:80] <- test[16:30,51:80]+rnorm(15*30, 3)
#' PermBiclust.sigclust_stop(test, silent=TRUE)
#' @export
#' @name PermBiclust.sigclust_stop
#' @author Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair





PermBiclust.sigclust_stop <- function(x, nperms=1000, silent=TRUE, maxnum.bicluster=5, alpha=0.05,icovest=1,sc=TRUE)
{
	x <- scale(x)
	spcl <- SqrtKMeansSparseCluster(x, K=2, wbounds=sqrt(ncol(x)), silent=silent)
    if (length(spcl)==1) {spcl <- spcl[[1]]}
	ws.perm <- matrix(NA, nrow=nperms, ncol=ncol(x))
	for(i in 1:nperms){
		ws.perm[i,] <- sqrt(rbeta(ncol(x), .5, (ncol(x)-1)/2))
		ws.perm[i,] <- sort(ws.perm[i,])
	}
	sws <- sort(spcl$ws)
	 spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])

if (min(table(spcl$Cs))<(2)){message('Warning: number of significant observations too small')
sc.pval=1} else 
if (sum(spcl$ws>sws[sig.ndx])<(2)){message('Warning: number of significant features too small')
sc.pval=1
}else{
	data<-x[,which(spcl$ws>sws[sig.ndx])]
sc.pval<-ifelse(sc==FALSE,0,sigclust::sigclust(data,nsim=1000,labflag=1, label=spcl$Cs,icovest=icovest)@pval)}

	
    which.x <- list()
	which.y <- list()
   
	if(sc.pval >= alpha | sum(spcl$ws>sws[sig.ndx])<2){
		which.x[[1]] <- rep(NA, nrow(x))
		which.y[[1]] <- rep(NA, ncol(x))
		return(list(num.bicluster=0, x.residual=x, which.x=which.x, which.y=which.y))
	} else{
		num.bicluster <- 1            ## primary bicluster
	    which.x[[1]] <- spcl$Cs==1
	    if (sum(spcl$Cs==2)<sum(which.x[[1]])) {which.x[[1]] <- spcl$Cs==2}
	    which.y[[1]] <- (spcl$ws>sws[sig.ndx])

		while(num.bicluster < maxnum.bicluster){
			mean.signal <- colMeans(x[which.x[[num.bicluster]], which.y[[num.bicluster]]])
			mu <- colMeans(x[!which.x[[num.bicluster]], which.y[[num.bicluster]]])
			delta <- mean.signal - mu
			diff <- matrix(rep(delta, each=sum(which.x[[num.bicluster]])), nrow=sum(which.x[[num.bicluster]]))
			x[which.x[[num.bicluster]], which.y[[num.bicluster]]] <- x[which.x[[num.bicluster]], which.y[[num.bicluster]]] - diff

			x <- scale(x)
			spcl <- SqrtKMeansSparseCluster(x, K=2, wbounds=sqrt(ncol(x)), silent=silent)
	        if (length(spcl)==1) {spcl <- spcl[[1]]}
			ws.perm <- matrix(NA, nrow=nperms, ncol=ncol(x))
			for(i in 1:nperms){
				ws.perm[i,] <- sqrt(rbeta(ncol(x), .5, (ncol(x)-1)/2))
				ws.perm[i,] <- sort(ws.perm[i,])
			}
		sws <- sort(spcl$ws)
	      spcl.diff <- colMeans(ws.perm)-sws
	      sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])

if (min(table(spcl$Cs))<(2)){message('Warning: number of significant observations too small')
sc.pval=1} else 
if (sum(spcl$ws>sws[sig.ndx])<(2)){message('Warning: number of significant features too small')
sc.pval=1
}else{
	data<-x[,which(spcl$ws>sws[sig.ndx])]
sc.pval<-ifelse(sc==FALSE,0,sigclust::sigclust(data,nsim=1000,labflag=1, label=spcl$Cs,icovest=icovest)@pval)}

			if(sc.pval >= alpha | sum((spcl$ws>sws[sig.ndx]))<2){
				message('Warning: bicluster', num.bicluster+1, 'does not exist \n')
				return(list(num.bicluster=num.bicluster, x.residual=x, which.x=which.x, which.y=which.y))
			} else{
			
				num.bicluster <- num.bicluster+1
				which.x[[num.bicluster]] <- spcl$Cs==1
			    if (sum(spcl$Cs==2)<sum(which.x[[num.bicluster]])) {which.x[[num.bicluster]] <- spcl$Cs==2}
			    which.y[[num.bicluster]] <- (spcl$ws>sws[sig.ndx])
			}
		}
		
		mean.signal <- colMeans(x[which.x[[num.bicluster]], which.y[[num.bicluster]]])
		mu <- colMeans(x[!which.x[[num.bicluster]], which.y[[num.bicluster]]])
		delta <- mean.signal - mu
		diff <- matrix(rep(delta, each=sum(which.x[[num.bicluster]])), nrow=sum(which.x[[num.bicluster]]))
		x[which.x[[num.bicluster]], which.y[[num.bicluster]]] <- x[which.x[[num.bicluster]], which.y[[num.bicluster]]] - diff
	    return(list(num.bicluster=num.bicluster, x.residual=x, which.x=which.x, which.y=which.y))
	}
}



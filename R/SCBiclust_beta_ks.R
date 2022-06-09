#' 'SCBiclust' method for identifying means-based biclusters with Kolmogorov-Smirnov test of feature weights
#'
#' @param x a dataset with n rows and p columns, with observations in rows. 
#' @param nperms number of \eqn{Beta(\frac{1}{2}, (p-1)/2)} distributed variables generated for each feature (default=1000)
#' @param silent should progress be printed? (default=TRUE)
#' @param maxnum.bicluster The maximum number of biclusters returned  
#' @param ks.alpha significance level for Kolmogorov-Smirnov test.
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
#' Observations in the bicluster are identified such that they maximize the feature-weighted square-root of the between cluster sum of squares.
#' Features in the bicluster are identified based on their contribution to the clustering of the observations. 
#' Feature weights are generated in a similar fashion as  \code{\link{KMeansSparseCluster}} 
#' except with a modified objective function and no sparsity constraint.
#' 
#' This algoritm uses a numerical approximation to  \eqn{E(\sqrt{B})} where \eqn{B \sim Beta(\frac{1}{2}, (p-1)/2)} as the expected null 
#' distribution for feature weights. The Kolmogorov-Smirnov test is used to assess if feature weights deviate from the expected null distribution.
#'  
#' @examples
#' test <- matrix(rnorm(100*200), nrow=100, ncol=200)
#' test[1:20,1:20] <- test[1:20,1:20]+rnorm(20*20, 2)
#' test[16:30,51:80] <- test[16:30,51:80]+rnorm(15*30, 3)
#' PermBiclust.beta.ks(test, silent=TRUE)
#' @export
#' @name PermBiclust.beta.ks
#' @author Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair



PermBiclust.beta.ks <- function(x, nperms=1000, silent=TRUE, maxnum.bicluster=5, ks.alpha=0.05)
{
	x <- scale(x)
	spcl <- SqrtKMeansSparseCluster(x, K=2, wbounds=sqrt(ncol(x)), silent=silent)
    if (length(spcl)==1) {spcl <- spcl[[1]]}
	ws.perm <- matrix(NA, nrow=nperms, ncol=ncol(x))
	for(i in 1:nperms){
		ws.perm[i,] <- sqrt(stats::rbeta(ncol(x), .5, (ncol(x)-1)/2))
		ws.perm[i,] <- sort(ws.perm[i,])
	}

	sws <- sort(spcl$ws)
    which.x <- list()
	which.y <- list()
    spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])

if (min(table(spcl$Cs))<2){message('Warning: number of significant observations <2')
ks.pval=1} else 
if (sum(spcl$ws>sws[sig.ndx])<2){message('Warning: number of significant features <2')
ks.pval=1
}else{
  ks.pval <- stats::ks.test(x=sws, y=ws.perm)$p.value}
	
	which.x <- list()
	which.y <- list()

	if(ks.pval >= ks.alpha | sum(spcl$ws>sws[sig.ndx])<2| min(table( spcl$Cs))<2 ){           
		which.x[[1]] <- rep(1, nrow(x))
		which.y[[1]] <- rep(0, ncol(x))
		return(list(num.bicluster=0, x.residual=x, which.x=which.x, which.y=which.y))
	}
	else{
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
				ws.perm[i,] <- sqrt(stats::rbeta(ncol(x), .5, (ncol(x)-1)/2))
				ws.perm[i,] <- sort(ws.perm[i,])
			}

	sws <- sort(spcl$ws)
    spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])

if (min(table(spcl$Cs))<2){message('Warning: number of significant observations <2')
ks.pval=1} else 
if (sum(spcl$ws>sws[sig.ndx])<2){message('Warning: number of significant features <2')
ks.pval=1
}else{
ks.pval <- stats::ks.test(x=sws, y=ws.perm)$p.value}

			if(ks.pval >= ks.alpha | sum((spcl$ws>sws[sig.ndx]))<2| min(table( spcl$Cs))<2){    
				message('Warning: bicluster', num.bicluster+1, 'does not exist \n')       
				return(list(num.bicluster=num.bicluster, x.residual=x, which.x=which.x, which.y=which.y))
			}
			else{
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



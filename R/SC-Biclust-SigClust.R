#' 'SCBiclust' method for identifying means-based biclusters
#'
#' @param x a dataset with n rows and p columns, with observations in rows. 
#' @param nperms number of \eqn{Beta(\frac{1}{2}, (p-1)/2)} distributed variables generated for each feature (default=1000)
#' @param silent should progress be printed? (default=TRUE)
#' @param maxnum.bicluster The maximum number of biclusters returned  
#' @param alpha significance level for \code{\link{sigclust}} test.
#' @param icovest Coviariance estimation type for \code{\link{sigclust}} test
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
#' distribution for feature weights. The \code{\link{sigclust}} algorithm is used to test the strength of the identified clusters.
#'  
#' @examples
#'  test <- matrix(rnorm(60*180), nrow=60, ncol=180)
#' test[1:15,1:15] <- test[1:15,1:15]+rnorm(15*15, 2)
#' test[16:30,51:80] <- test[16:30,51:80]+rnorm(15*30, 3)
#' PermBiclust.sigclust(test, silent=TRUE)
#' @export
#' @name PermBiclust.sigclust
#' @author Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair 



PermBiclust.sigclust <- function(x, nperms=1000, silent=TRUE, maxnum.bicluster=5, alpha=0.05,icovest=1)
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
	 spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])
	data<-x[,which(spcl$ws>sws[sig.ndx])]


if (ncol(as.matrix(data))<2){message('Warning: number of significant features equals',ncol(as.matrix(data)))
sc.pval=1
}else{
sc.pval <- sigclust::sigclust(data,nsim=1000,labflag=1, label=spcl$Cs,icovest=icovest)@pval}

	
    which.x <- list()
	which.y <- list()
   
	if(sc.pval >= alpha | sum(spcl$ws>sws[sig.ndx])<2){
		message('Warning: no bicluster found \n')
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
				ws.perm[i,] <- sqrt(stats::rbeta(ncol(x), .5, (ncol(x)-1)/2))
				ws.perm[i,] <- sort(ws.perm[i,])
			}
		sws <- sort(spcl$ws)
	      spcl.diff <- colMeans(ws.perm)-sws
	      sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])
	      data<-x[,which(spcl$ws>sws[sig.ndx])]

if (ncol(as.matrix(data))<2){message('Warning: number of significant features equals',ncol(as.matrix(data)))
sc.pval=1
}else{
sc.pval <- sigclust::sigclust(data,nsim=1000,labflag=1, label=spcl$Cs,icovest=icovest)@pval}

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

################################################################################
## helper functions

BCSS<-function(x,CS){
  N=by(x,CS,nrow)
  CM=by(x,CS,colMeans)
  M=colMeans(x)
  (N[[1]]*(CM[[1]]-M)^2)+(N[[2]]*(CM[[2]]-M)^2)
}

SqrtUpdateWs <- function (x, Cs, l1bound)
{
  bcss.perfeature <- BCSS(x,Cs)
  return(sqrt(bcss.perfeature)/sqrt(sum(bcss.perfeature)))
}

UpdateCs<-function (x, K, ws, Cs) 
{
  x <- x[, ws != 0]
  z <- sweep(x, 2, sqrt(ws[ws != 0]), "*")
  nrowz <- nrow(z)
  mus <- NULL
  if (!is.null(Cs)) {
    for (k in unique(Cs)) {
      if (sum(Cs == k) > 1) 
        mus <- rbind(mus, apply(z[Cs == k, ], 2, mean))
      if (sum(Cs == k) == 1) 
        mus <- rbind(mus, z[Cs == k, ])
    }
  }
  if (is.null(mus)) {
    km <- stats::kmeans(z, centers = K, nstart = 10)
  }
  else {
    distmat <- as.matrix(stats::dist(rbind(z, mus)))[1:nrowz, (nrowz + 
                                                          1):(nrowz + K)]
    nearest <- apply(distmat, 1, which.min)
    if (length(unique(nearest)) == K) {
      km <- stats::kmeans(z, centers = mus)
    }
    else {
      km <- stats::kmeans(z, centers = K, nstart = 10)
    }
  }
  return(km$cluster)
}

SqrtKMeansSparseCluster <- function (x, K = NULL, wbounds = NULL, nstart = 20, silent = FALSE,
                                     maxiter = 10, centers = NULL)
{
  if (is.null(K) && is.null(centers))
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers)) {
    if (nrow(centers) != K)
      stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if (nrow(centers) == K)
      K <- NULL
  }
  if (!is.null(centers) && ncol(centers) != ncol(x))
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  if (is.null(wbounds))
    wbounds <- seq(1.1, sqrt(ncol(x)), len = 20)
  if (min(wbounds) <= 1)
    stop("wbounds should be greater than 1")
  wbounds <- c(wbounds)
  out <- list()
  if (!is.null(K))
    Cs <- stats::kmeans(x, centers = K, nstart = nstart)$cluster
  if (is.null(K))
    Cs <- stats::kmeans(x, centers = centers)$cluster
  for (i in 1:length(wbounds)) {
    if (length(wbounds) > 1 && !silent)
      message(i, fill = FALSE)
    ws <- rep(1/sqrt(ncol(x)), ncol(x))
    ws.old <- stats::rnorm(ncol(x))
    store.bcss.ws <- NULL
    niter <- 0
    while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 &&
           niter < maxiter) {
      if (!silent)
        message(niter, fill = FALSE)
      niter <- niter + 1
      ws.old <- ws
      if (!is.null(K)) {
        if (niter > 1)
          Cs <- UpdateCs(x, K, ws, Cs)
      }
      else {
        if (niter > 1)
          Cs <- UpdateCs(x, nrow(centers), ws, Cs)
      }
      ws <- SqrtUpdateWs(x, Cs, wbounds[i])
      store.bcss.ws <- c(store.bcss.ws, sum(BCSS(x,Cs)* ws))
    }
    out[[i]] <- list(ws = ws, Cs = Cs, wbound = wbounds[i])
  }
  if (!silent)
    message(fill = TRUE)
  class(out) <- "KMeansSparseCluster"
  return(out)
}



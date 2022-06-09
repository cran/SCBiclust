#' 'SCBiclust' method for identifying variance-based biclusters
#'
#' @param x a dataset with n rows and p columns, with observations in rows. 
#' @param min.size Minimum size of observations included in a valid bicluster (default=\code{max(5,round(nrow(x)/20)))}
#' @param nperms number of \eqn{\chi^2_{n_1}} and \eqn{\chi^2_{n_2}} variables generated for each feature
#'  where \eqn{n_1} and \eqn{n_2} are the number of observations in cluster 1 and cluster 2, respectively. (default=100)
#' @param silent should progress be printed? (default=TRUE)
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
#' Observations in the bicluster are identified such that they maximize the feature-weighted sum of between cluster difference in feature variances.
#' Features in the bicluster are identified based on their contribution to the clustering of the observations. 
#' This algoritm uses a numerical approximation \eqn{log(abs(\chi^2_{n_1}-chi^2_{n_2})+1)} as the expected null 
#' distribution for feature weights. 
#'
#' \code{VarPermBiclust.chisqdiff} will identify at most one variance bicluster. To identify additional biclusters first the feature signal
#' of the identified bicluster should be removed by scaling the variance of elements in the previously identified bicluster, Then 
#' \code{VarPermBiclust.chisqdiff} can be used on the residual data matrix. (see example)
#' @examples
#' test <- matrix(rnorm(100*50, mean=1, sd=2), nrow=100)
#' test[1:30, 1:20] <- matrix(rnorm(30*20, mean=1, sd=15), nrow=30)
#' test.VarPermBiclust <- VarPermBiclust.chisqdiff(test)
#' x=test.VarPermBiclust$which.x
#' y=test.VarPermBiclust$which.y
#' # Code for identifying additional biclusters after removing bicluster signal
#' \donttest{
#' temp <- scale(test)
#' temp[x,y] <-t(t(temp[x,y])*(apply(temp[!x,y],2,sd)/
#'                               apply(temp[x,y],2,sd)))
#' test.VarPermBiclust.2 <- VarPermBiclust.chisqdiff(temp)
#' }
#' @export
#' @name VarPermBiclust.chisqdiff
#' @author Erika S. Helgeson, Qian Liu, Guanhua Chen, Michael R. Kosorok , and Eric Bair





VarPermBiclust.chisqdiff <- function(x, min.size=max(5,round(nrow(x)/20)), nperms=1000, silent=TRUE){
	## assign initial CI based on obs variance
	x <- scale(x)
	CI.old <- rep(1, nrow(x))
	obsvar <- apply(x, 1, stats::var)
	obsvar.ord <- order(obsvar)
	cut <- round(nrow(x)/2)
	ind <- c(1:nrow(x))
	y <- ind[obsvar.ord]
	CI.old[1:cut] <- CI.old[1:cut]+1
	CI.old <- CI.old[order(y)]                                                    ## initial CI
	ws.old <- var.updateWS(x, CI.old)                                             ## initial ws
	BCVarDiff.old <- var.diff(x, CI.old, ws.old)                                  ## sum of between cluster variance difference
	#ws.new <- rep(1/sqrt(ncol(x)), ncol(x))
	CI.new <- var.updateCI.2(x, min.size, ws.old)
	ws.new <- var.updateWS(x, CI.new)

    if((sum(abs(ws.new - ws.old))/sum(abs(ws.old))) <= 1e-04){
		CI.new <- CI.old
	}

	while((sum(abs(ws.new - ws.old))/sum(abs(ws.old))) > 1e-04){
		ws.old <- ws.new
		CI.new <- var.updateCI.2(x, min.size, ws.old)
		ws.new <- var.updateWS(x, CI.new)                                                   ## update ws
	}

	ref.ws <- ws.new
	ws.perm <- matrix(NA, nrow=nperms, ncol=ncol(x))
	for(i in 1:nperms){
		chisqdiff <- rchisqdiff(n=ncol(x), df1=sum(CI.new==1)-1, df2=sum(CI.new==2)-1)
		LogDiff <- log(abs(chisqdiff)+1)
		ws.perm[i,] <- LogDiff/sqrt(sum(LogDiff^2))
	}

	sws <- sort(ref.ws)
    spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-
                             spcl.diff[2:length(spcl.diff)])

    which.x <- CI.new==1
    if (sum(CI.new==2)<sum(which.x)) {which.x <- CI.new==2}
    which.y <- (ref.ws>sws[sig.ndx])

	return(list(which.x=which.x, which.y=which.y))
}


VarPermBiclust.chisqdiff.ks <- function(x, min.size=max(5,round(nrow(x)/20)), nperms=1000, maxnum.bicluster=5, ks.alpha=0.05, silent=TRUE){
	## assign initial CI based on obs variance
	x <- scale(x)
	CI.old <- rep(1, nrow(x))
	obsvar <- apply(x, 1, stats::var)
	obsvar.ord <- order(obsvar)
	cut <- round(nrow(x)/2)
	ind <- c(1:nrow(x))
	y <- ind[obsvar.ord]
	CI.old[1:cut] <- CI.old[1:cut]+1
	CI.old <- CI.old[order(y)]                                                    ## initial CI
	ws.old <- var.updateWS(x, CI.old)                                             ## initial ws
	BCVarDiff.old <- var.diff(x, CI.old, ws.old)                                  ## sum of between cluster variance difference
	CI.new <- var.updateCI.2(x, min.size, ws.old)
	ws.new <- var.updateWS(x, CI.new)
    if((sum(abs(ws.new - ws.old))/sum(abs(ws.old))) <= 1e-04){
		CI.new <- CI.old
	}
	while((sum(abs(ws.new - ws.old))/sum(abs(ws.old))) > 1e-04){
		ws.old <- ws.new
		CI.new <- var.updateCI.2(x, min.size, ws.old)
		ws.new <- var.updateWS(x, CI.new)                                                   ## update ws
	}

	ref.ws <- ws.new
	ws.perm <- matrix(NA, nrow=nperms, ncol=ncol(x))
	for(i in 1:nperms){
		chisqdiff <- rchisqdiff(n=ncol(x), df1=sum(CI.new==1)-1, df2=sum(CI.new==2)-1)
		LogDiff <- log(abs(chisqdiff)+1)
		ws.perm[i,] <- LogDiff/sqrt(sum(LogDiff^2))
		ws.perm[i,] <- sort(ws.perm[i,])
	}

	sws <- sort(ref.ws)
	ks.pval <- stats::ks.test(x=sws, y=ws.perm)$p.value
    which.x <- list()
	which.y <- list()
    spcl.diff <- colMeans(ws.perm)-sws
	sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])
	
	if(ks.pval >= ks.alpha | sum(ref.ws>sws[sig.ndx])<2){
		message('Warning: no variance bicluster found \n')
		which.x[[1]] <- rep(1, nrow(x))
		which.y[[1]] <- rep(0, ncol(x))
		return(list(num.bicluster=0, x.residual=x, which.x=which.x, which.y=which.y))
	}
	
	else{
		num.bicluster <- 1
	    which.x[[1]] <- CI.new==1
	    if (sum(CI.new==2)<sum(which.x[[1]])) {which.x[[1]] <- CI.new==2}
	    which.y[[1]] <- (ref.ws>sws[sig.ndx])
		
		while(num.bicluster < maxnum.bicluster){
				x[which.x[[num.bicluster]],which.y[[num.bicluster]]] <- t(t(x[which.x[[num.bicluster]],which.y[[num.bicluster]]]) 
			                                                        * (apply(x[!which.x[[num.bicluster]],which.y[[num.bicluster]]], 2, stats::sd) 
																	/ apply(x[which.x[[num.bicluster]],which.y[[num.bicluster]]], 2, stats::sd)))
			x <- scale(x)
			CI.old <- rep(1, nrow(x))
			obsvar <- apply(x, 1, stats::var)
			obsvar.ord <- order(obsvar)
			cut <- round(nrow(x)/2)
			ind <- c(1:nrow(x))
			y <- ind[obsvar.ord]
			CI.old[1:cut] <- CI.old[1:cut]+1
			CI.old <- CI.old[order(y)]                                                    ## initial CI
			ws.old <- var.updateWS(x, CI.old)                                             ## initial ws
			BCVarDiff.old <- var.diff(x, CI.old, ws.old)                                  ## sum of between cluster variance difference
			CI.new <- var.updateCI.2(x, min.size, ws.old)
			ws.new <- var.updateWS(x, CI.new)
		    if((sum(abs(ws.new - ws.old))/sum(abs(ws.old))) <= 1e-04){
				CI.new <- CI.old
			}
			while((sum(abs(ws.new - ws.old))/sum(abs(ws.old))) > 1e-04){
				ws.old <- ws.new
				CI.new <- var.updateCI.2(x, min.size, ws.old)
				ws.new <- var.updateWS(x, CI.new)                                                   ## update ws
			}
			
			ref.ws <- ws.new
			ws.perm <- matrix(NA, nrow=nperms, ncol=ncol(x))
			for(i in 1:nperms){
				ws.perm[i,] <- LogDiff/sqrt(sum(LogDiff^2))
				ws.perm[i,] <- sort(ws.perm[i,])
			}
			sws <- sort(ref.ws)
			ks.pval <- stats::ks.test(x=sws, y=ws.perm)$p.value
		    spcl.diff <- colMeans(ws.perm)-sws
			sig.ndx <- which.max(spcl.diff[1:(length(spcl.diff)-1)]-spcl.diff[2:length(spcl.diff)])

			if(ks.pval >= ks.alpha | sum((ref.ws>sws[sig.ndx]))<2){
				message('Warning: variance bicluster', num.bicluster+1, 'does not exist \n')
				return(list(num.bicluster=num.bicluster, x.residual=x, which.x=which.x, which.y=which.y))
			}
			else{
				num.bicluster <- num.bicluster+1
				which.x[[num.bicluster]] <- CI.new==1
			    if (sum(CI.new==2)<sum(which.x[[num.bicluster]])) {which.x[[num.bicluster]] <- CI.new==2}
			    which.y[[num.bicluster]] <- (ref.ws>sws[sig.ndx])
			}
		}
		x[which.x[[num.bicluster]],which.y[[num.bicluster]]] <- t(t(x[which.x[[num.bicluster]],which.y[[num.bicluster]]]) 
		                                                        * (apply(x[!which.x[[num.bicluster]],which.y[[num.bicluster]]], 2, stats::sd) 
																/ apply(x[which.x[[num.bicluster]],which.y[[num.bicluster]]], 2, stats::sd)))
		
		return(list(num.bicluster=num.bicluster, x.residual=x, which.x=which.x, which.y=which.y))
	}
}

################################################################################
## helper functions

var.diff <- function(x, CI, ws) {
  x <- scale(x)
  var1 <- apply(x[CI==1,], 2, stats::var)%*%ws
  var2 <- apply(x[CI==2,], 2, stats::var)%*%ws
  return(log(abs(var1-var2)+1))
}

rchisqdiff <- function(n, df1, df2){
  ## function to generate random variables following the same distribution as (X_df1^2 - X_df2^2)
  chisq.1 <- stats::rchisq(n, df1, ncp=0)/df1
  chisq.2 <- stats::rchisq(n, df2, ncp=0)/df2
  return(chisq.1-chisq.2)
}

var.updateCI <- function(x, min.size, ws) {
  x.Cs <- rep(1, nrow(x))
  n2 <- sample(min.size:(nrow(x)-min.size))[1]
  x.Cs[sample(1:nrow(x))[1:n2]] <- 2
  cur.diff <- var.diff(x, x.Cs, ws)
  last.diff <- 0
  while (cur.diff>last.diff) {
    last.diff <- cur.diff
    for (i in 1:nrow(x)) {
      x.Cs2 <- x.Cs
      if (x.Cs2[i]==1) {
        x.Cs2[i] <- 2
      }
      else {
        x.Cs2[i] <- 1
      }
      if ((sum(x.Cs2==1)>=min.size)&(sum(x.Cs2==2)>=min.size)) {
        cur.diff2 <- var.diff(x, x.Cs2, ws)
        if (cur.diff2>cur.diff) {
          x.Cs <- x.Cs2
          cur.diff <- cur.diff2
        }
      }
    }
  }
  return(x.Cs)
}

var.updateCI.2 <- function(x, min.size, ws) {
  x.Cs <- rep(1, nrow(x))
  obsvar <- apply(x, 1, stats::var)
  obsvar.ord <- order(obsvar)
  cut <- round(nrow(x)/2)
  ind <- c(1:nrow(x))
  y <- ind[obsvar.ord]
  x.Cs[1:cut] <- x.Cs[1:cut]+1
  x.Cs <- x.Cs[order(y)]
  cur.diff <- var.diff(x, x.Cs, ws)
  last.diff <- 0
  while (cur.diff>last.diff) {
    last.diff <- cur.diff
    for (i in 1:nrow(x)) {
      x.Cs2 <- x.Cs
      if (x.Cs2[i]==1) {
        x.Cs2[i] <- 2
      }
      else {
        x.Cs2[i] <- 1
      }
      if ((sum(x.Cs2==1)>=min.size)&(sum(x.Cs2==2)>=min.size)) {
        cur.diff2 <- var.diff(x, x.Cs2, ws)
        if (cur.diff2>cur.diff) {
          x.Cs <- x.Cs2
          cur.diff <- cur.diff2
        }
      }
    }
  }
  return(x.Cs)
}


var.updateWS <- function(x, CI){         ## update ws based on log of between cluster variance difference
  ColVar1 <- apply(x[CI==1,], 2, stats::var)
  ColVar2 <- apply(x[CI==2,], 2, stats::var)
  LogVarDiff <- log(abs(ColVar1-ColVar2)+1)
  ws.LogVarDiff <- LogVarDiff/sqrt(sum(LogVarDiff^2))
  return(ws.LogVarDiff)
}
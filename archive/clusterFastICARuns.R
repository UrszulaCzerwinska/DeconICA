#' @importFrom foreach "%dopar%"
#' This function runs the fastICA algorithm several times with random initializations.
#' The obtained components are clustered and
#' the medoids of these clusters are used as the final estimates. The returned estimates are ordered by decreasing Iq values which measure the compactness of the clusters (see details).
#'
#'
#' This function implements in R fastICA iterations followed by a clustering step, as defined in the matlab package 'icasso'.
#' Among the indices computed by icasso,
#' only the Iq index is currently computed. As defined in 'icasso', the Iq index measures the difference between the intra-cluster similarity and the extra-cluster similiarity.
#' No visualization of the clusters is yet available.
#'
#' If \code{bootstrap=TRUE} a bootstrap (applied to the observations) is used to perturb the data before each iteration, then function \code{fastICA} is applied with random initializations.
#'
#' By default, in 'icasso', agglomerative hierarchical clustering with average linkage is performed. To use the same clustering, please use \code{funClus="hclust"} and \code{method="average"}. But this function also allows you to apply the clustering of your choice among \code{kmeans, pam, hclust, agnes} by specifying \code{funClus} and adding the adequat additional parameters.
#'
#' See details of the functions \code{\link[fastICA]{fastICA}}.
#' @title Run of fastICA and JADE algorithms
#' @param X A data matrix with n rows representing observations (e.g genes) and p columns representing variables (e.g samples).
#' @param nbComp The number of components to be extracted.
#' @param nbIt The number of iterations of FastICA
#' @param alg.type If \code{alg.type="parallel"} the components are extracted simultaneously (the default), if \code{alg.type="deflation"} the components are extracted one at a time, see \code{\link[fastICA]{fastICA}}.
#' @param fun The functional form of the G function used in the approximation to neg-entropy (see 'details' of the help of function \code{\link[fastICA]{fastICA}}).
#' @param row.norm a logical value indicating whether rows of the data matrix \code{X} should be standardized beforehand (see help of function \code{fastICA})
#' @param maxit The maximum number of iterations to perform.
#' @param tol A positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged.
#' @param funClus The clustering function to be used to cluster the estimates
#' @param bootstrap if TRUE the data is bootstraped before each fastICA iteration, else (default) only random initializations are done
#' @param ... Additional parameters for code{funClus}
#' @return A list consisting of: \describe{\item{A}{the estimated mixing matrix} \item{S}{the estimated source matrix}, item{W}{the estimated unmixing matrix}, \item{Iq}{Iq indices.}}
#' @author Anne Biton
#' @export
#' @examples
#' ## generate a data
#' set.seed(2004);
#' M <- matrix(rnorm(5000*6,sd=0.3),ncol=10)
#' M[1:100,1:3] <- M[1:100,1:3] + 2
#' M[1:200,1:3] <- M[1:200,4:6] +1
#'
#' ## Random initializations are used for each iteration of FastICA
#' ## Estimates are clustered using hierarchical clustering with average linkage
#' res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
#'                           nbIt=3, funClus="hclust", method="average")
#'
#' ## Data are boostraped before each iteration and random initializations
#' ## are used for each iteration of FastICA
#' ## Estimates are clustered using hierarchical clustering with ward
#' res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
#'                           nbIt=3, funClus="hclust", method="ward")
#'
#'
clusterFastICARuns <- function(X, nbComp, nbIt=100, alg.type = c("deflation", "parallel"), fun = c("logcosh","exp"), maxit = 500, tol = 10^-6, funClus= c("hclust","agnes","pam","kmeans"), row.norm = FALSE, bootstrap=FALSE, ...) {

  message(paste("FastICA iteration 1"))
  resICA <- fastICA2(X = X, n.comp = nbComp, alg.typ = alg.type, fun = fun, maxit = maxit, tol = tol, w.distribution = 'gaussian')

  ## whitening matrix
  whit <- resICA$K
  ## dewhitening matrix
  dewhit <- solve(t(resICA$K)%*%resICA$K) %*% t(resICA$K)

  it <-  clus <- NULL
  ## run fastICA x times after data are bootstrapped (at the gene level)
  allW <-
    foreach::foreach(it=2:nbIt, .combine=cbind) %dopar% {
      message(paste("FastICA iteration "),it)
      if (bootstrap)
        Xbis <- X[sample(1:ncol(X), replace=TRUE),]
      else
        Xbis <- X
      res <- fastICA2(X = Xbis, n.comp = nbComp, alg.typ = alg.type, fun = fun, maxit = maxit, tol = tol, w.distribution = 'gaussian')
      res$W
    }


  allW <- cbind(resICA$W,allW)
  ## project W in original space (dewhitening)
  allWdewith <- t(allW)%*%dewhit
  allWdewith <- (apply(allWdewith,2,scale))


  ## compute similarity between W = absolute correlation values
  sim <- abs(cor(t(allWdewith)))
  dsim <- 1-sim
  centrotypes <- c()

  switch(funClus,
         hclust={
           resClus <- stats::hclust(d=stats::as.dist(dsim), ...)
           partition <- stats::cutree(resClus, k=nbComp)
         },
         agnes={
           resClus <- cluster::agnes(x=dsim, diss=TRUE, ...)
           partition <-  stats::cutree(stats::as.hclust(resClus), k=nbComp)
         },
         pam={
           resClus <- cluster::pam(x=dsim, diss=TRUE, k=nbComp, keep.diss=FALSE, ...)
           partition <- resClus$clustering
           centrotypes <- resClus$medoids
         },
         kmeans={
           resClus <- stats::kmeans(x=allW, centers=nbComp, ...)
           partition <- resClus$cluster
         }
  )

  ## compute Iq indices and extract centrotypes
  # Iq=avg(intra-cluster similarity) - avg(extra-cluster similarity)

  Iq <-
    foreach::foreach(clus=unique(partition), .combine=c) %dopar% {
      indC <- which(partition==clus)
      if (length(indC)>1) {
        if (funClus != "pam")
          centrotypes <- c(centrotypes,indC[which.max(apply(sim[indC,indC],1,sum))])
        internalSim <- mean(sim[indC,indC])
      }
      else {
        if (funClus != "pam")
          centrotypes <- c(centrotypes,indC)
        internalSim <- sim[indC,indC]
      }
      externalSim <- mean(sim[indC,setdiff(1:ncol(sim),indC)])
      iq <- internalSim-externalSim

      return(iq)
    }


  ## extract W including the centrotypes of each cluster
  W <- whit %*% allW[,centrotypes]
  A <- solve(t(W)%*%W) %*% t(W)
  S <- as.matrix(X)%*%as.matrix(W)
  rownames(S) <- rownames(X)
  colnames(A) <- colnames(X)

  orderIq <- order(Iq, decreasing=TRUE)
  return(list(A=A[orderIq,],S=S[,orderIq],W=W[,orderIq],Iq=Iq[orderIq]))

}

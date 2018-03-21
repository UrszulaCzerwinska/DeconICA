
# The fastICA2 is a slightly modified version of fastICA http://cran.r-project.org/web/packages/fastICA/index.html
fastICA2<-function (X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh","exp"),
                    w.distribution = c("uniform","gaussian","beta"), alpha = 1, maxit = 200, tol = 1e-04, verbose = FALSE)
{
  # The additional parameter is tot_var, the variance of dataset X preserved for ICA
  # and w.distribution, from which random distribut
  # In order to improve the readability, the variables n.comp,row.norm & verbose,w.init were ommited

  ### Input X, same as original fastICA
  dd <- dim(X)
  d <- dd[dd != 1L]
  if (length(d) != 2L)
    stop("data must be matrix-conformal")
  X <- if (length(d) != length(dd))
    matrix(X, d[1L], d[2L]) else
      as.matrix(X)


  # Since we want to keep at least 3 PCs, so the pourcentage should not be smaller than cum_var_3

  ### Input alg.typ & fun & w.distribution, same as original fastICA
  #alg.typ <- match.arg(alg.typ)
 # fun <- match.arg(fun)
  #w.distribution<- match.arg(w.distribution)

  ### Input alpha, same as original fastICA
  if (alpha < 1 || alpha > 2)
    stop("alpha must be in range [1,2]")

  ### Calculate the number of PCs preserved from tot_var
  #n.comp <- match.arg(n.comp)

  ### Data pretreatment: denoising
  # Unlike original script, PCA was used for denoising instead of SVD on covariance matrix
  prx <- prcomp(X,scale=F)
  D <- diag(c(1/prx$sdev))
  K <- D %*% t(prx$rotation)
  K <- matrix(K[1:n.comp,],n.comp,ncol(X))
  X=t(X)
  Xd <- K %*% X

  ### Generation of w.init by the defined distribution
  if (w.distribution=='uniform'){w.init=matrix(runif(n.comp*n.comp,0,n.comp),n.comp,n.comp)}
  if (w.distribution=='gaussian'){w.init=matrix(rnorm(n.comp*n.comp),n.comp,n.comp)}
  if (w.distribution=='beta'){w.init=matrix(rbeta(n.comp*n.comp,1,n.comp),n.comp,n.comp)}

  ### FastICA algorithm applied on the denoised matrix Xd, same as original fastICA
  a <- if (alg.typ == "deflation")
    ica.R.def(
      Xd,
      n.comp,
      tol = tol,
      fun = fun,
      alpha = alpha,
      maxit = maxit,
      verbose = verbose,
      w.init = w.init
    ) else if (alg.typ == "parallel")
    ica.R.par(
      Xd,
      n.comp,
      tol = tol,
      fun = fun,
      alpha = alpha,
      maxit = maxit,
      verbose = verbose,
      w.init = w.init
    )

  ### Calculation of source and loading matrix & function output, same as original fastICA
  w <- a %*% K
  S <- w %*% X # Source matrix S = a*K*X
  A <- t(w) %*% solve(w %*% t(w)) # Loading matrix is the pseudo inverse of matrix w
  ### Output
  return(list(X = t(X), K = t(K), W = t(a), A = t(A), Xd = t(Xd), S = t(S), W0 = w.init, IC = n.comp))
}

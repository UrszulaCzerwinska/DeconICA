### removes mean row-wise
# x - row
.center_rowmeans <- function(x) {
  xcenter <- rowMeans(x)
  x - rep(xcenter, ncol(x))
}

### it orient ICs so that the long tail is positive
# S - ica S matrix
.orient_funct <- function(S) {
  orient <-
    apply(S, 2, function(x) {
      if (min(x) < -3 & max(x) > 3) {
        ifelse (sum(x > 3)  < sum(x < -3), -1, 1)
      } else {
        ifelse (sum(x > 2)  < sum(x < -2), -1, 1)
      }
    })
  S <- as.matrix(S)  %*% diag(orient)
  return(S)
}

### if there is no long tails orient by max correlation (so that the max corr are +)
# S - ica S matrix
# r - correlation matrix
.orient_max <- function(S, r) {
  col.as <-
    apply(data.frame(r), 1, function(ic)
      ifelse(abs(min(ic)) > max(ic), -1, 1))
  S <- as.matrix(S)  %*% diag(col.as)
  return(S)
}

### computes variance of each row
# x - row
.rowVars <- function (x, na.rm = TRUE) {
  x <- as.matrix(x)
  sqr <- function(x)
    x * x
  n <- rowSums(!is.na(x))
  n[n <= 1] <- NA
  rowSums(sqr(x - rowMeans(x, na.rm = na.rm)), na.rm = na.rm) / (n - 1)
}


### computes variance of each column
# x  - column
.colVars <- function (x, na.rm = TRUE) {
  .rowVars(t(x))
}

### computes cumulative variance of ICA
# pca - pca object
.cumVar <- function(pca) {
  cumsum(pca$sdev ^ 2 / sum(pca$sdev ^ 2))
}

### computes norm of a vector
# x - vector
.norm_vec <- function(x)
  sqrt(sum(x ^ 2))

### remove multiplicated genes based on expression values
# X - gene matrix, first column is gene names
.remove_duplicates <- function(X) {
  message("for duplicated genes, entry with higher variance will be kept")
  X[["var"]] <- .rowVars(X[, 2:ncol(X)])
  X.o <-  X[order(X[, 1], -X[["var"]]), ]
  X <- X.o[!duplicated(X.o[, 1]), ]
  X[["var"]] <- NULL
  return (X)
}

### wrapper of intersect function
.intersect.genes <-
  function(res, component)
    intersect(res[, 1], component[, 1])

### selects genes from a list of metagenes to build correlation matrix
.corr_matrix <- function (y, res, name) {
  res[[name]] <- NULL
  res[match(.intersect.genes(res, y), res[, 1]), name] <-
    y[match(.intersect.genes(res, y), y[, 1]), 2]
  return(res)
}

### double checks if there is enough genes to corelate give input condition
# n - matrix
# n.genes.interest - min number of genes to correlate
.verify.n <-
  function(n, n.genes.intersect)
    apply(n, 1, function(ic)
      any(ic > n.genes.intersect))

### computes fraction of variance explained for an IC
# S - ica S matrix in norm vectors,
# X input matrix
# i - component
.fev <- function(S.norm, X, i) {
  stats::var(S.norm[, i]) / sum(.colVars(X))
}


### apply distribution test over rows of a matrix
# df - matrix of gene expression
# array1 - names of columns in group/condition 1 or indexes
# array2 - names of columns in group/condition 2 or indexes
# test - test function

## example
# test <- METABRIC[1:10,1:20]
# array1 <- sample(colnames(test), size =10)
# array2 <- setdiff(colnames(test), array1)
# .compute_x_test_for_many(test, array1, array2)
.compute_x_test_for_many <-
  function (df, array1, array2, test = "t.test") {
    res <-
      apply(df, 1, function(x)
        do.call(test, list(x[array1], x[array2])))
    p.value <- sapply(res, function(x)
      x[[3]])
    stat <- sapply(res, function(x)
      x[[1]])
    res.df <- data.frame(cbind(p.value, stat))
    row.names(res.df) <- row.names(df)
    return(res.df)
  }

### geometric mean
# x- vector
gm_mean <- function(x, na.rm = TRUE, zero.propagate = FALSE){
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero.propagate){
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}

### harmonic mean
# a - vector
harmonic_mean <- function(a, ...){
  1/mean(1/a, ...)
}


### positive gaussian
# p - vector of probabilities.
# mean - vector of means
# sd - vector of standard deviations.
# lwr - lower bound
# upr - upper bound
# rounding - round() function parameter

.pos.gaussian <- function(p, mean, sd, lwr, upr, rounding= Inf) {
  samp <- round(stats::rnorm(p, mean, sd), rounding)
  samp[samp < lwr] <- lwr
  samp[samp > upr] <- upr
  samp
}

#' Verify if data is in log scale
#'
#' @param x \code{data.frame} or \code{matrix}
#'
#' @return
#' \code{TRUE} or \code{FALSE}
#' @export
#'
#' @examples
#' M <- matrix(sample(-1:14, 100, replace = TRUE),10,10, byrow = TRUE)
#' is_logscale(M)
#' M2 <- 2^M
#' is_logscale(M2)
is_logscale <- function(x){
  qx <- as.numeric(stats::quantile(x, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  !LogC
}


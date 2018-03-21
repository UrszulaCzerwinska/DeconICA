#'Decompose dataset with ICA.
#'
#'This is a wrapper of \code{\link[fastICA]{fastICA}}. It allows compute number
#'of ICs optimal for over decomposition for immune deconvolution.
#'
#'@param X a data matrix with \code{n} rows representing observations and
#'  \code{p} columns representing variables, place gene names in the first
#'  column and select \code{with.names = TRUE} to keep gene names for further
#'  analysis
#'@param row.center if \code{TRUE} substract row mean from data
#'@param with.names if first column of X is row.names please indicate
#'  \code{TRUE}, in case of duplicated names, the transcript with highest
#'  variance will be kept, names need to be HUGO names, if names are not
#'  provided at this step, you can provide them later
#'@param gene.names character vector of row names - gene names
#'@param R if TRUE (default) the R version of fastICA is running, else the
#'  matlab version (you need to provide parametrs of your matlab engine)
#'@param optimal check \code{TRUE} to let select best number of components for
#'  deconvolution, for datasets >120 columns, n.comp will be set to 100, if <120
#'  then number of components will be selected according to Kaiser Rule (90
#'  percent of variance explained)
#'@inheritParams fastICA::fastICA
#'@param isLog if data is in log \code{TRUE} if data is in counts \code{FALSE}
#'@param path_global only if \code{R = FALSE}, the global path where files will
#'  be written, current directory by default
#'@param matlbpth only if \code{R = FALSE},  the path to matlab engine, for Mac
#'  : MATLAB_XXX.app/Contents/MacOS/MATLAB_maci64
#'@param fasticapth path to repository of source matlab
#'@param ... other possible parameters from \code{\link[fastICA]{fastICA}}
#'@return A list containing the following components as in
#'  \code{\link[fastICA]{fastICA}} \describe{
#'  \item{X}{pre-processed data matrix
#'  (after PCA)}
#'  \item{K}{pre-whitening matrix that projects data onto the first
#'  n.comp principal components.}
#'  \item{W}{estimated un-mixing matrix (see definition in details)}
#'  \item{A}{estimated mixing matrix}
#'  \item{S}{estimated source matrix}
#'   \item{names}{if  \code{with.names = TRUE} will contain row
#'  names list}
#'  \item{counts}{if  \code{isLog = FALSE} will contain initial
#'  matrix without duplicated genes}
#'  \item{log.counts}{initial matrix without
#'  duplicated genes in log2(x+1) before centering}
#'  \item{samples}{sample names as provided}
#'  }
#' @export
#' @examples
#' \dontrun{
#' # numerical matrix
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- data.frame(S %*% A)
#' run_fastica(X, row.center = TRUE, n.comp = 2, optimal = FALSE)
#' # matrix with gene names
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- data.frame(S %*% A)
#' names <- paste("A",1:nrow(X), sep="")
#' X <- cbind(names,X)
#' run_fastica(X, row.center = TRUE, n.comp = 2, optimal = FALSE, with.names = TRUE)
#' }
#'@seealso \code{\link[fastICA]{fastICA}}
#'  \url{https://cran.r-project.org/web/packages/fastICA/index.html}
run_fastica <-
  function(X,
           optimal = TRUE,
           row.center = TRUE,
           with.names = FALSE,
           gene.names = NULL,
           alg.typ = "parallel",
           method = "C",
           n.comp = 100,
           isLog = TRUE,
           R = TRUE,
           path_global = getwd(),
           matlbpth = "/Applications/MATLAB_R2016a.app/Contents/MacOS/MATLAB_maci64",
           fasticapth = "/Users/ulala/Documents/CURIE/BIODICA-master/BIODICA/bin/fastica++",
           ...) {
    #get.dataset.name
    df.name <- deparse(substitute(X))
    if (!is.null(gene.names)) {
      if (length(gene.names) != nrow(X)) {
        stop("wrong gene.names")
      } else {
        X <- data.frame(gene.names = gene.names, X)
        with.names <- TRUE
      }
    }
    # If gene names included, check if there are duplicates
    if (with.names) {
      if (any(duplicated(X[, 1]))) {
        X <- .remove_duplicates(X)
      }
      names <- X[, 1]
      X <- X[, 2:ncol(X)]
    }
    #define colnames
    if (!is.null(colnames(X))) {
      samples <- colnames(X)
    } else {
      samples <- make.names(seq(1, ncol(X), 1))
    }
    # if matrix is wide it might mean user didn't transpose the matrix
    if (nrow(X) < ncol(X))
      warning(
        "row length is less than column length :
        place genes in rows and samples in columns for deconvolution"
      )
    if (!isLog) {
      counts  <- X
      X <- log2(X + 1)
      message("transforming data to log2")
    }
    else {
      log.counts <- X
    }
    # center rows if needed
    if (row.center)
      X <- .center_rowmeans(X)
    if (all.equal(round(mean(as.matrix(X[1, ])), 2), 0) != TRUE)
      message("your data is not mean centerd by rows (genes)")

    if (optimal) {
      if (ncol(X) < 120) {
        message("applying Kaiser Rule")
        pca <- stats::prcomp(X, center = FALSE, scale. = FALSE)
        .cumVar(pca)
        n.pca <- which(.cumVar(pca) > 0.9)[1]
        n.comp <- n.pca
        message(paste("optimal number of components should be : ", n.pca, sep =
                        ""))
      } else {
        n.comp <- 100
      }
    }
    # you cannot get more signals than columns
    if (n.comp > ncol(X)) {
      n.comp <- ncol(X)
      message(paste(
        "number of independent components was lowered to max (",
        n.comp,
        ")",
        sep =
          ""
      ))
    }
    #running ICA with R::fastica
    if (R) {
      # denoise with PCA
      message("running PCA")
      pca <- stats::prcomp(X, center = FALSE, scale. = FALSE)
      X <- pca$x[, 1:ncol(X)]
      # run ICA
      message(paste("running ICA for ", n.comp, " components", sep = ""))
      X.ica <-
        fastICA::fastICA(
          X = X,
          n.comp = n.comp,
          alg.typ = alg.typ,
          method = method,
          maxit = 1000,
          tol = 1e-09,
          ...
        )
    } else {
      #running ICA with Matlab with icasso stabilisation
      #detect operating system
      sys <- switch(Sys.info()[["sysname"]],
             Windows = "Windows",
             Linux  = "Linux",
             Darwin = "Mac")
      message("Running MATLAB ICA")
      X.ica <-
        .doICA(
          X,
          names = names,
          samples = samples,
          path_global = path_global,
          n = n.comp,
          name = df.name,
          corr_folder = "CORRELATION",
          matlbpth,
          fasticapth
        )
    }
    # add names if with.names
    message("adding names to the object")
    if (with.names) {
      X.ica[["names"]] <- as.character(names)
    }
    #add sample names
    message("adding sample names to the object")
    X.ica[["samples"]] <- samples
    #add counts if counts provided
    if (!isLog) {
      message("adding raw.counts to the object")
      X.ica[["counts"]] <- counts
    }
    #add log.counts
    message("adding counts in log to the object")
    X.ica[["log.counts"]] <- log.counts
    return(X.ica)
    message("ready")
  }

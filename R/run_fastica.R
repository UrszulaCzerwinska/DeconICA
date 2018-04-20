#'Decompose dataset with ICA.
#'
#'This is a wrapper of \code{\link[fastICA]{fastICA}}. It allows compute number
#'of ICs overdecompose for over decomposition for immune deconvolution.
#'
#'@param X a data matrix with \code{n} rows representing observations and
#'  \code{p} columns representing variables, place gene names in the first
#'  column and select \code{with.names = TRUE}
#'@param row.center if \code{TRUE} subtract row mean from data
#'@param with.names if first column of X is row.names please indicate
#'  \code{TRUE}, in case of duplicated names, the transcript with highest
#'  variance will be kept, names need to be HUGO names, if names are not
#'  provided at this step, you can provide them later
#'@param gene.names character vector of row names - gene names
#'@param samples if samples names different from column names
#'@param R if TRUE (default) the R version of fastICA is running, else the
#'  matlab version (you need to provide parameters of your matlab engine)
#'@param overdecompose check \code{TRUE} to let select best number of components for
#'  deconvolution, for datasets >120 columns, n.comp will be set to 100, if <120
#'  then number of components will be selected according to Kaiser Rule (90
#'  percent of variance explained)
#'@inheritParams fastICA::fastICA
#'@param isLog if data is in log \code{TRUE} if data is in counts \code{FALSE}
#'@param path_global only if \code{R = FALSE}, the global path where files will
#'  be written, current directory by default
#'@param matlbpth only if \code{R = FALSE},  the path to matlab engine, it uses
#'\code{\link[matlabr:get_matlab]{get_matlab}} to find path to your matlab automatically
#'@param export.corr \code{TRUE} if you need to export \code{S} matrix in a specific
#'format for correlation in external java app
#'@param fasticapth path to repository of source matlab code, it is set by default
#'as coming with the package
#'@param name important for Matlab version, defines the name of your files
#'@param ... other possible parameters for \code{\link[fastICA]{fastICA}}
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
#'
#' # numerical matrix
#' S <- matrix(runif(10000), 10, 2)
#' A <- matrix(sample(-3:3, 16, replace = TRUE),2,8, byrow = TRUE)
#' X <- data.frame(S %*% A)
#' run_fastica(X, row.center = TRUE, n.comp = 2, overdecompose = FALSE)
#' #matlab
#' \dontrun{
#'run_fastica(X, row.center = TRUE, n.comp = 3, overdecompose = FALSE, R = FALSE)
#'}
#' # matrix with gene names
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- data.frame(S %*% A)
#' names <- paste("A",1:nrow(X), sep="")
#' X <- cbind(names,X)
#' run_fastica(X, row.center = TRUE, n.comp = 2, overdecompose = FALSE, with.names = TRUE)
#'
#'@seealso \code{\link[fastICA]{fastICA}}
#'  \url{https://cran.r-project.org/web/packages/fastICA/index.html}
run_fastica <-
  function(X,
           overdecompose = TRUE,
           row.center = TRUE,
           with.names = FALSE,
           gene.names = NULL,
           samples = NULL,
           alg.typ = "parallel",
           method = "C",
           n.comp = 100,
           isLog = TRUE,
           R = TRUE,
           path_global = getwd(),
           matlbpth = NULL,
           fasticapth = paste0(path.package("deconica", quiet = TRUE), "/fastica++"),
           export.corr = FALSE,
           name = NULL,
           ...) {
    #get.dataset.name
    if(is.null(name)){
      df.name <- deparse(substitute(X))
      } else {
        df.name <- name }

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
    } else {
      names <- make.names(seq(1, nrow(X), 1))
    }

    #define colnames
    if(is.null(samples)){
      if (!is.null(colnames(X))) {
        samples <- colnames(X)
      } else {
        samples <- make.names(seq(1, ncol(X), 1))
      }
    } else {
      samples <- samples
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
      log.counts <- X
    }
    else {
      log.counts <- X
    }
    # center rows if needed
    if (row.center)
      X <- .center_rowmeans(X)
    if (all.equal(round(mean(as.matrix(X[1,])), 2), 0) != TRUE)
      message("your data is not mean centerd by rows (genes)")

    if (overdecompose) {
      if (ncol(X) < 120) {
        message("applying Kaiser Rule")
        pca <- stats::prcomp(X, center = FALSE, scale. = FALSE)
        .cumVar(pca)
        n.pca <- which(.cumVar(pca) > 0.9)[1]
        n.comp <- n.pca
        message(paste("overdecompose number of components should be : ", n.pca, sep =
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
      message("Running MATLAB ICA")
      X.ica <-
        doICA(
          X,
          names = names,
          samples = samples,
          path_global = path_global,
          n = n.comp,
          name = df.name,
          export.corr = export.corr,
          corr_folder = "CORRELATION",
          matlbpth = matlbpth,
          fasticapth = fasticapth
        )
    }
    colnames(X.ica$A) <- samples
    row.names(X.ica$A) <- paste0("IC", 1:nrow(X.ica$A))
    colnames(X.ica$S) <- paste0("IC", 1:ncol(X.ica$S))
    row.names(X.ica$S) <- names
    # add names if with.names
    message("adding names to the object")
    X.ica[["names"]] <- as.character(names)
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



#' Repriduce process of run_fastica
#'
#' Applies preprocessig of \code{\link{run_fastica}} but instead of running ICA
#' it imports matlab output files. It is handy if you run the matlab idenpendly
#' or if you lost R session data
#'
#' @inheritParams run_fastica
#' @param import imports data only if TRUE
#'
#' @return an object as \code{\link{run_fastica}}
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # numerical matrix
#' S <- matrix(runif(10000), 10, 2)
#' A <- matrix(sample(-3:3, 16, replace = TRUE),2,8, byrow = TRUE)
#' X <- data.frame(S %*% A)
#' #matlab
#'run_fastica(X, row.center = TRUE, n.comp = 2, overdecompose = FALSE, R = FALSE)
#' run_fastica_import(X, row.center = TRUE, n.comp = 2, overdecompose = FALSE, import=TRUE)
#' }
#'
#'
run_fastica_import <-
  function(X,
           overdecompose = TRUE,
           row.center = TRUE,
           with.names = FALSE,
           gene.names = NULL,
           n.comp = 100,
           isLog = TRUE,
           import = TRUE,
           path_global = getwd(),
           name = NULL,
           ...) {
    #get.dataset.name
    if(is.null(name)){
      df.name <- deparse(substitute(X))
    } else {
      df.name <- name }

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
    } else {
      names <- make.names(seq(1, nrow(X), 1))
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
      log.counts <- X
    }
    else {
      log.counts <- X
    }
    # center rows if needed
    if (row.center)
      X <- .center_rowmeans(X)
    if (all.equal(round(mean(as.matrix(X[1,])), 2), 0) != TRUE)
      message("your data is not mean centerd by rows (genes)")

    if (overdecompose) {
      if (ncol(X) < 120) {
        message("applying Kaiser Rule")
        pca <- stats::prcomp(X, center = FALSE, scale. = FALSE)
        .cumVar(pca)
        n.pca <- which(.cumVar(pca) > 0.9)[1]
        n.comp <- n.pca
        message(paste("overdecompose number of components should be : ", n.pca, sep =
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
    if(import){
      name <- paste0(df.name, "_", n.comp)
      X.ica <- import_ICA_res (name = name, ncomp = n.comp, path_global_1 = paste0(path_global,"/", df.name, "_", n.comp,"/"))
      #running ICA with Matlab with icasso stabilisation
      #detect operating system
    } else {
      X.ica <- NULL
    }
    colnames(X.ica$A) <- samples
    row.names(X.ica$A) <- paste0("IC", 1:nrow(X.ica$A))
    colnames(X.ica$S) <- paste0("IC", 1:ncol(X.ica$S))
    row.names(X.ica$S) <- names
    # add names if with.names
    message("adding names to the object")
    X.ica[["names"]] <- as.character(names)
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

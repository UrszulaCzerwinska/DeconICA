#' Call \code{doICA} matlab function
#'
#' function used inside \code{\link{run_fastica}} to run \code{fastICA} with
#' \code{icasso} stabilization. Matlab engine is necessary
#'
#' @param df.scaled.t scaled numerical data matrix
#' @param names gene names, no duplicates
#' @param samples sample names
#' @param path_global path where files will be saved
#' @param n number of components
#' @param name \code{FALSE} by default, name of dataset is used, you can put
#'   your name
#' @param export.corr \code{FALSE} by default, if you want to use a java
#'   correlation function later or select \code{TRUE}
#' @param corr_folder \code{"CORRELATION"} by default, only if you selected
#'   \code{export.corr = TRUE}
#' @param matlbpth is found automatically with \code{\link{get_matlab}} function,
#'   replace if not functional
#' @param fasticapth path to \code{fastica++} repository with MATLAB scripts
#'
#' @return it returns \code{A}, \code{S} matrices of ICA and \code{names} and
#' \code{samples} for coherence
#' @export
#'
#' @seealso \code{\link[matlabr]{get_matlab}}, \code{\link{run_fastica}},
#' \code{\link{export_for_ICA}},\code{\link[matlabr]{run_matlab_code}},
#' \code{\link{import_ICA_res}}, code{\link{export_for_correlation_java}}
#' @examples
#'\dontrun{
#'data(Example_ds)
#'res.pre <-
#'  prepare_data_for_ica(Example_ds[, -1], names = Example_ds[, 1])
#'res.do <- doICA(
#'   df.scaled.t =  res.pre$df.scaled,
#'   names = res.pre$names,
#'   samples = res.pre$samples,
#'   path_global = getwd(),
#'   n = 5,
#'   name = "test",
#'   export.corr = FALSE
#')
#'}
#'
doICA <-
  function(df.scaled.t,
           names,
           samples,
           path_global = getwd(),
           n,
           name = FALSE,
           export.corr = FALSE,
           corr_folder = "CORRELATION",
           matlbpth = NULL,
           fasticapth = paste0(path.package("deconica", quiet = TRUE), "/fastica++")) {
    if (is.null(matlbpth) &
        (!(matlabr::have_matlab())))
      stop("Matlab could not be found on your disk \n provide path to your matlab in \'matlbpth\'")
    path.init <- getwd()
    res.exp <- export_for_ICA(
      df.scaled.t = df.scaled.t,
      names = names,
      samples = samples,
      path_global = path_global,
      name = name,
      n = n
    )
    path_global_1 <- res.exp$path_global_1
    name <- res.exp$name
    fun = "doICA"

    ncomp = n

    path = paste("'", path_global_1, "'", sep = "")

    file = paste("'", paste(name, "_numerical.txt", sep = ""), "'", sep =
                   "")

    fasticapth = fasticapth

    cd <- paste0("cd('", fasticapth, "');")

    cmd <- paste0(fun, "(", path, ",", file, ",", ncomp, ")")

    run_matlab_code_2(c(cd, cmd), matcmd = matlbpth)

    res.imp <-
      import_ICA_res(name = name,
                   ncomp = ncomp,
                   path_global_1 = path_global_1)
    if (export.corr)
      export_for_correlation_java(
        path_global_1 = path_global_1,
        corr_folder = corr_folder,
        names =  names,
        S = res.imp$S,
        samples = samples,
        A = res.imp$A,
        ncomp = ncomp,
        name = name
      )
    setwd(path.init)
    return(list(
      A = t(res.imp$A),
      S = res.imp$S,
      names = names,
      samples = samples
    ))
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Import results of ICA
#'
#' imports files run in Matlab or precomputed
#'
#' @param name name of the dataset
#' @param ncomp number of components
#' @param path_global_1 absolute path of the files
#'
#' @return
#' imports \code{A} and \code{S} ICA matrix
#' @seealso \code{\link{run_fastica}}, \code{\link{export_for_ICA}},
#' \code{\link{doICA}}, \code{\link{export_for_correlation_java}}
#' @examples
#'\dontrun{
#' data(Example_ds)
#'res.pre <-
#'  prepare_data_for_ica(Example_ds[, -1], names = Example_ds[, 1])
#'res.do <- doICA(
#'   df.scaled.t =  res.pre$df.scaled,
#'   names = res.pre$names,
#'   samples = res.pre$samples,
#'   path_global = getwd(),
#'   n = 5,
#'   name = "test",
#'   export.corr = FALSE
#')
#' import_ICA_res("test_5", 5, paste0(getwd(),"/test_5/"))
#' }
import_ICA_res <- function(name, ncomp, path_global_1) {
  a.imp = paste(paste("A", paste0(name, "_numerical.txt", sep = ""), ncomp, sep =
                        "_"), "num", sep = ".")
  a.imp.path = paste(path_global_1, a.imp, sep = "")
  A = utils::read.delim(a.imp.path, header = FALSE, sep = "\t")
  s.imp = paste(paste("S", paste(name, "_numerical.txt", sep = ""), ncomp, sep =
                        "_"), "num", sep = ".")
  s.imp.path = paste(path_global_1, s.imp, sep = "")
  S = utils::read.delim(s.imp.path, header = FALSE, sep = "\t")

  A = A[, 1:(ncol(A) - 1)]
  S = S[, 1:(ncol(S) - 1)]
  return(list(A = A, S = S))
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Exports \code{S} ICA matrix in a specific format
#'
#' needed for an external function in java
#'
#' @param corr_folder export folder name \code{"CORRELATION"} by default
#' @param names gene names
#' @param S \code{S} ICA matrix
#' @param samples sample names
#' @param A \code{A} ICA matrix
#' @param ncomp number of computed components
#' @param name name of the dataset
#' @param path_global_1 absolute path
#'
#' @return
#' saves on the drive in \code{corr_folder} exported files
#' @export
#'
#' @seealso \code{\link{run_fastica}}, \code{\link{import_ICA_res}},
#' \code{\link{doICA}}, \code{\link{export_for_ICA}}
#' @examples
#' \dontrun{
#'data(Example_ds)
#'res.pre <-
#'  prepare_data_for_ica(Example_ds[, -1], names = Example_ds[, 1])
#'res.do <- doICA(
#'   df.scaled.t =  res.pre$df.scaled,
#'   names = res.pre$names,
#'   samples = res.pre$samples,
#'   path_global = getwd(),
#'   n = 5,
#'   name = "test",
#'   export.corr = FALSE
#')
#' export_for_correlation_java(
#'S = res.do$S,
#'A = t(res.do$A),
#'names =  res.do$names,
#'samples = res.do$samples,
#'name = "test",
#'ncomp = 5
#')
#'}
export_for_correlation_java <-
  function(corr_folder = "CORRELATION",
           names,
           S,
           samples,
           A,
           ncomp,
           name,
           path_global_1 = getwd()) {
    path1 <- paste0(path_global_1,"/", corr_folder,"/" )
    dir.create(path1)
    utils::write.table(
      cbind(names, S),
      file = paste(path1, paste(
        paste(name, ncomp, sep = "."), "_S", ".xls", sep = ""
      ), sep = "") ,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
    utils::write.table(
      data.frame(SAMPLES = samples, A),
      file = paste0(
        path_global_1,
        paste(name, ncomp, "full_samples", "txt", sep = ".")
      ) ,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
    utils::write.table(
      data.frame(GENES = names, cbind(S)),
      file = paste0(
        path_global_1,
        paste(name, ncomp, "full_genes", "txt", sep = ".")
      ) ,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )

  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Export files
#'
#' export files in right format to run fastICA in MATLAB or BiodICA
#'
#' @param df.scaled.t scaled numerical matrix
#' @param names gene names, vector of character string
#' @param samples sample names, vector of character string
#' @param path_global path to export files, current directory by default
#' @param name name of the dataset
#' @param n number of components
#'
#' @return
#' writes files on the drive in indicated location
#' @seealso \code{\link{run_fastica}}, \code{\link{import_ICA_res}},
#' \code{\link{doICA}}, \code{\link{export_for_correlation_java}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(Example_ds)
#' res.pre <-
#'  prepare_data_for_ica(Example_ds[, -1], names = Example_ds[, 1])
#' export_for_ICA(res.pre$df.scaled,
#'res.pre$names,
#'res.pre$samples,
#'path_global = getwd(),
#'name ="test",
#'n = 5)
#' }
export_for_ICA <-
  function(df.scaled.t,
           names,
           samples,
           path_global = getwd(),
           name = FALSE,
           n = "") {
    setwd(path_global)
    if (name != FALSE) {
      name <- paste0(name, "_", n)
    } else {
      name <- paste0(deparse(substitute(df.scaled.t)), "_", n)
    }

    dir.create(name)

    path_global_1 = paste0(path_global, "/", name, "/")

    utils::write.table(
      df.scaled.t,
      file = paste0(path_global_1, name, "_numerical.txt"),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    utils::write.table(
      names,
      file = paste0(path_global_1, name, "_ids.txt"),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )

    utils::write.table(
      samples,
      file = paste0(path_global_1, name, "_samples.txt"),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    return(list(path_global_1 = path_global_1, name = name))
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Formats data for ICA in MATLAB
#'
#' @param df numerical data matrix
#' @param names gene names character vector
#' @param samples if not provided column names will be used
#'
#' @return
#' \item{\code{df.scaled}}{scaled data without duplicates}
#' \item{\code{names}}{gene names without duplicates}
#' \item{\code{non.scaled}}{non scaled data without duplicates}
#' \item{\code{samples}}{sample names}
#'
#' @seealso \code{\link{run_fastica}}, \code{\link{import_ICA_res}},
#' \code{\link{doICA}}, \code{\link{doICABatch}}
#' @export
#'
#' @examples
#'
#' data(Example_ds)
#' prepare_data_for_ica(Example_ds[, -1], names = Example_ds[, 1])
prepare_data_for_ica <- function(df, names, samples = NULL) {
  X <- data.frame(Genes = names, df)
  if (any(duplicated(X[, 1]))) {
    X <- .remove_duplicates(X)
    names <- X[, 1]
    X <- as.matrix(X[, 2:ncol(X)])
    X.ndup <- X
  } else {
    names <- X[, 1]
    X <- as.matrix(X[, 2:ncol(X)])
    X.ndup <- X
  }
  if (is.null(samples)) {
    samples <- as.character(colnames(df))
  } else {
    samples <- samples
  }
  # center rows if needed
  X <- .center_rowmeans(X)
  if (all.equal(round(mean(as.matrix(X[1, ])), 2), 0) != TRUE)
    message("your data is not mean centerd by rows (genes)")
  return(
    list(
      df.scaled = X,
      non.scaled = X.ndup ,
      names = as.character(names),
      samples = as.character(samples)
    )
  )
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' \code{doBatchICA}
#'
#' prepares the data (scales and removes duplicates), runs \code{doBatchICA.m} MATLAB script
#'
#' @param df numerical data matrix
#' @param vec vector of values for which ICA should be computed
#' @param path_global path were files will be saved, current directory by default
#' @param names gene names
#' @param samples sample names
#' @param name name of the dataset, if not provided, name of R variable
#' @param matlbpth path to matlab, found automatically with \code{\link{get_matlab}}
#' @param fasticapth path to \code{fastica++}
#'
#' @return
#' plots of stability and MSTD if possible
#' @export
#'
#' @examples
#' \dontrun{
#' data(Example_ds)
#' doICABatch(
#'   Example_ds[, -1],
#'   seq(2, 4, 1),
#'   names = Example_ds[, 1],
#'   samples = colnames(Example_ds[, -1]),
#'   name = "test",
#'   fasticapth = paste0(path.package("deconica", quiet = FALSE), "/fastica++")
#' )
#' }
doICABatch <-
  function(df,
           vec,
           path_global = getwd(),
           names,
           samples,
           name = NULL,
           matlbpth = NULL,
           fasticapth = paste0(path.package("deconica", quiet = TRUE), "/fastica++")) {
    if (is.null(matlbpth) &
        (!(matlabr::have_matlab())))
      stop("Matlab could not be found on your disk \n provide path to your matlab in \'matlbpth\'")
    res.pre <-
      prepare_data_for_ica(df = df,
                           names = names,
                           samples  = samples)

    n = paste(vec[1], vec[length(vec)], sep = "_")
    path.init <- getwd()
    res.exp <- export_for_ICA(
      df.scaled.t = res.pre$df.scaled,
      names = res.pre$names,
      samples = res.pre$samples,
      path_global = path_global,
      name = name,
      n = n
    )
    path_global_1 <- res.exp$path_global_1
    name <- res.exp$name
    ncomp = paste("[", paste(as.character(vec), collapse = " "), "]", sep =
                    " ")
    fun = "doICABatch"
    path = paste("'", path_global_1, "'", sep = "")

    file = paste("'", paste(name, "_numerical.txt", sep = ""), "'", sep =
                   "")
    matlbpth = matlbpth
    fasticapth = fasticapth
    cd <- paste0("cd('", fasticapth, "');")
    cmd <- paste0(fun, "(", path, ",", file, ",", ncomp, ")")
    run_matlab_code_2(c(cd, cmd), matcmd = matlbpth)
    t.imp.path = paste(path_global_1, "avg_stability.plot.txt", sep = "")
    T = utils::read.delim(t.imp.path, header = FALSE, sep = "\t")
    row.names(T) = as.character(T[, 1])
    T = T[order(T[, 1]),]
    bp <-
      graphics::barplot(
        t(T[, 2, drop = FALSE]),
        col = 'blue',
        ylim = c(0, 1),
        main = paste0(name, " ", "stability")
      )
    im1 <- paste0(path_global_1, "_MSTD_estimate.png")
    im2 <- paste0(path_global_1, "_AverageStability.png")
    p2 <- p3 <- NULL
    if (file.exists(im1)) {
      img1 <- png::readPNG(im1)
      p2 <- grid::grid.raster(img1)
    }
    if (file.exists(im2)) {
      img2 <- png::readPNG(im2)
      p3 <- grid::grid.raster(img1)
    }
    return(list(
      stability = T,
      stability.p = bp,
      MSTD.p = p2,
      Av.stability.p = p3
    ))
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' @title Run matlab script
#'
#' @description This function runs a matlab script, and
#' returns exit statuses, slighly modified from matlabr
#' @param fname Filename of matlab script (.m file)
#' @param verbose print diagnostic messages
#' @param ... Options passed to \code{\link{system}}
#' @param matcmd path to matlab engine
#' @inheritParams matlabr::get_matlab
#' @export
#' @return Exit status of matlab code
run_matlab_script_2 = function(fname,
                               matcmd = NULL,
                               verbose = TRUE,
                               desktop = FALSE,
                               splash = FALSE,
                               display = FALSE,
                               wait = TRUE,
                               ...) {
  stopifnot(file.exists(fname))
  if (is.null(matcmd))
    matcmd <- matlabr::get_matlab(
      desktop = desktop,
      splash = splash,
      display = display,
      wait = wait
    )
  cmd <- paste0(
    ' "',
    "try, run('",
    fname,
    "'); ",
    "catch err, disp(err.message); ",
    "exit(1); end; exit(0);",
    '"'
  )
  cmd <- paste0(matcmd, cmd)
  if (verbose) {
    message("Command run is:")
    message(cmd)
  }
  x <- system(cmd, wait = wait, ...)
  return(x)
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' @title Runs matlab code
#'
#' @description This function takes in matlab code, where
#' the last line must end with a ;, and returns the exit
#' status
#' @param code Character vector of code.
#' @param endlines Logical of whether the semicolon (;) should be
#' pasted to each element of the vector.
#' @param verbose Print out filename to run
#' @param add_clear_all Add \code{clear all;} to the beginning of code
#' @param paths_to_add Character vector of PATHs to add to the
#' script using \code{\link{add_path}}
#' @param ... Options passed to \code{\link{run_matlab_script}}
#' @param matcmd path to matlab engine
#' @export
#' @return Exit status of matlab code
#' @examples
#' if (matlabr::have_matlab()){
#'    run_matlab_code_2("disp(version)")
#'    run_matlab_code_2("disp(version)", paths_to_add = "~/")
#'    run_matlab_code_2(c("disp('The version of the matlab is:')", "disp(version)"))
#'    run_matlab_code_2(c("x = 5", "disp(['The value of x is ', num2str(x)])"))
#' }
run_matlab_code_2 = function(code,
                             matcmd = NULL,
                             endlines = TRUE,
                             verbose = TRUE,
                             add_clear_all = FALSE,
                             paths_to_add = NULL,
                             ...) {
  code <- c(ifelse(add_clear_all, "clear all;", ""),
            paste0("cd('", getwd(), "');"),
            code)
  if (!is.null(paths_to_add)) {
    paths_to_add = add_path(paths_to_add)
    code = c(code, paths_to_add)
  }
  sep <- ifelse(endlines, ";", " ")
  code <- paste0(code, sep = sep, collapse = "\n")
  code <- gsub(";;", ";", code)
  #   cmd <- paste(' "try \n')
  #   cmd <- paste(cmd, code)
  #   cmd <- paste(cmd, "\n catch err \n disp(err.message); \n exit(1); \n")
  #   cmd <- paste0(cmd, 'end; \n exit(0);"')
  #   cmd = gsub("\n", ";", cmd)
  #   cmd = paste0(matcmd, cmd)
  cmd <- code
  fname <- tempfile(fileext = ".m")
  cat(cmd, file = fname)
  if (verbose) {
    message(paste0("Script created: ", fname))
  }
  x <-
    run_matlab_script_2(fname, matcmd = matcmd, verbose = verbose, ...)
  return(x)
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Create PATHs to add to MATLAB PATHs
#'
#' @param path path to add
#'
#' @return A character vector
#' @examples
#' add_path("~/")
#' @export
add_path = function(path) {
  path = sapply(path, function(x) {
    paste0("addpath('", path, "');")
  })
  path = unname(unlist(path))
  return(path)
}

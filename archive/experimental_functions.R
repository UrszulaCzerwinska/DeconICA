#
#
# head(Example_ds)
# dim(Example_ds[, -1])
# dim(Example_ds)
# res.pre <-
#   .prepare_data_for_ica(Example_ds[, -1], names = Example_ds[, 1])
#
# res.do <- .doICA(
#   df.scaled.t =  res.pre$df.scaled,
#   names = res.pre$names,
#   samples = res.pre$samples,
#   path_global = getwd(),
#   n = 5,
#   name = "test",
#   export.corr = TRUE
# )
# str(res.do)

doICA <-
  function(df.scaled.t,
           names,
           samples,
           path_global = getwd(),
           n,
           name = FALSE,
           export.corr = FALSE,
           corr_folder = "CORRELATION",
           matlbpth = get_matlab(),
           fasticapth = paste0(path.package("deconica", quiet = TRUE), "/inst/fastica++")) {
    path.init <- getwd()
    res.exp <- .export_for_ICA(
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
    matlbpth = matlbpth
    fasticapth = fasticapth

    cd <- paste0("cd('", fasticapth, "');")

    cmd <- paste0(fun, "(", path, ",", file, ",", ncomp, ")")

    run_matlab_code(c(cd, cmd))

    res.imp <-
      importICAres(name = name,
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

importICAres <- function(name, ncomp, path_global_1) {
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

export_for_correlation_java <-
  function(corr_folder,
           names,
           S,
           samples,
           A,
           ncomp,
           name,
           path_global_1) {
    path1 <- paste(path_global_1, corr_folder, "/")
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


# data(Example_ds)
# doICABatch(
#   Example_ds[, -1],
#   seq(2, 4, 1),
#   names = Example_ds[, 1],
#   samples = colnames(Example_ds[, -1]),
#   name = "test",
#   fasticapth = paste0(path.package("deconica", quiet = FALSE), "/fastica++")
# )

doICABatch <-
  function(df,
           vec,
           path_global = getwd(),
           names,
           samples,
           name,
           matlbpth = get_matlab(),
           fasticapth = paste0(path.package("deconica", quiet = TRUE), "/inst/fastica++")) {
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

    run_matlab_code(c(cd, cmd))

    t.imp.path = paste(path_global_1, "avg_stability.plot.txt", sep = "")
    T = read.delim(t.imp.path, header = FALSE, sep = "\t")
    row.names(T) = as.character(T[, 1])
    T = T[order(T[, 1]), ]
    bp <-
      barplot(
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


#' @title Find matlab path
#'
#' @description This tries to find matlab's path using a system which
#' command, and then, if not found, looks at \code{getOption("matlab.path")}.  If not path is found, it fails.
#' @param try_defaults (logical) If \code{matlab} is not found from
#' \code{Sys.which}, and \code{matlab.path} not found, then try some
#' default PATHs for Linux and OS X.
#' @param desktop Should desktop be active for MATLAB?
#' @param splash Should splash be active for MATLAB?
#' @param display Should display be active for MATLAB?
#' @param wait Should R wait for the command to finish.  Both
#' passed to \code{\link{system}} and adds the \code{-wait} flag.
#' @author person(given = "John",
#'family = "Muschelli",
#'role = c("aut", "cre"),
#' email = "muschellij2@gmail.com")
#' @source https://github.com/muschellij2/matlabr/
#' @export
#' @return Character of command for matlab
#' @examples
#' if (have_matlab()) {
#' get_matlab()
#' }
get_matlab = function(
  try_defaults = TRUE,
  desktop = FALSE,
  splash = FALSE,
  display = FALSE,
  wait = TRUE){
  # find.matlab <- system("which matlab", ignore.stdout=TRUE)
  mat = paste0(
    "matlab",
    ifelse(
      .Platform$OS.type %in% "windows",
      ".exe",
      "")
  )
  find.matlab = as.numeric(Sys.which(mat) == "")
  myfunc = function(x, name) {
    x = as.logical(x)
    ifelse(x, "", paste0("-no", name))
  }
  desktop = myfunc(desktop, "desktop")
  splash = myfunc(splash, "splash")
  display = myfunc(display, "display")
  wait = ifelse(
    .Platform$OS.type %in% "windows",
    ifelse(wait, "-wait", ""),
    "")

  matcmd <- paste0(mat, " ",
                   wait, " ",
                   desktop, " ",
                   splash, " ",
                   display, " -r ")

  if (find.matlab != 0) {
    mpath = getOption("matlab.path")

    ####################################
    # Trying defaults
    ####################################
    if (is.null(mpath)) {
      if (try_defaults) {
        this.year = as.numeric(format(Sys.Date(), "%Y"))
        years = seq(this.year + 1, this.year - 5, by = -1)
        mac_ends = c(outer(years, c("a", "b"), paste0))
        def_paths = c(
          "/usr/local/bin",
          "/usr/bin",
          paste0("/Applications/MATLAB_R", mac_ends, ".app/bin"),
          paste0("C:/Program Files/MATLAB/R", mac_ends, "/bin"),
          paste0("D:/Program Files/MATLAB/R", mac_ends, "/bin")
        )
        for (ipath in def_paths) {
          def_path = file.path(ipath, mat)
          if (file.exists(def_path)) {
            warning(paste0("Setting matlab.path to ", ipath))
            options(matlab.path = ipath)
            mpath = ipath
            break;
          } # end def_path
        } # end loop
      } # end try_defaults
    } # end null mpath

    stopifnot(!is.null(mpath))
    stopifnot(file.exists(mpath))
    mpath = shQuote(mpath)
    matcmd <- file.path(mpath, matcmd)
  }
  return(matcmd)
}

#' @title Logical check if MATLAB is accessible
#'
#' @description Uses \code{\link{get_matlab}} to check if
#' MATLAB's path accessible
#' @export
#' @return Logical \code{TRUE} is MATLAB is accessible, \code{FALSE} if not
#' @examples
#' have_matlab()
have_matlab = function(){
  x = suppressWarnings(try(get_matlab(), silent = TRUE))
  return(!inherits(x, "try-error"))
}








#' @title Run matlab script
#'
#' @description This function runs a matlab script, and
#' returns exit statuses
#' @param fname Filename of matlab script (.m file)
#' @param verbose print diagnostic messages
#' @param ... Options passed to \code{\link{system}}
#' @inheritParams get_matlab
#' @export
#' @return Exit status of matlab code
run_matlab_script = function(
  fname,
  verbose = TRUE,
  desktop = FALSE,
  splash = FALSE,
  display = FALSE,
  wait = TRUE,
  ...){
  stopifnot(file.exists(fname))
  matcmd = get_matlab(
    desktop = desktop,
    splash = splash,
    display = display,
    wait = wait)
  cmd = paste0(' "', "try, run('", fname, "'); ",
               "catch err, disp(err.message); ",
               "exit(1); end; exit(0);", '"')
  cmd = paste0(matcmd, cmd)
  if (verbose) {
    message("Command run is:")
    message(cmd)
  }
  x <- system(cmd, wait = wait, ...)
  return(x)
}

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
#' @export
#' @return Exit status of matlab code
#' @examples
#' if (have_matlab()){
#'    run_matlab_code("disp(version)")
#'    run_matlab_code("disp(version)", paths_to_add = "~/")
#'    run_matlab_code(c("disp('The version of the matlab is:')", "disp(version)"))
#'    run_matlab_code(c("x = 5", "disp(['The value of x is ', num2str(x)])"))
#' }
run_matlab_code = function(
  code, endlines = TRUE, verbose = TRUE,
  add_clear_all = FALSE,
  paths_to_add = NULL,
  ...){
  # matcmd = get_matlab()
  code = c(ifelse(add_clear_all, "clear all;", ""), code)
  if (!is.null(paths_to_add)) {
    paths_to_add = add_path(paths_to_add)
    code = c(code, paths_to_add)
  }
  sep = ifelse(endlines, ";", " ")
  code = paste0(code, sep = sep, collapse = "\n")
  code = gsub(";;", ";", code)
  #   cmd <- paste(' "try \n')
  #   cmd <- paste(cmd, code)
  #   cmd <- paste(cmd, "\n catch err \n disp(err.message); \n exit(1); \n")
  #   cmd <- paste0(cmd, 'end; \n exit(0);"')
  #   cmd = gsub("\n", ";", cmd)
  #   cmd = paste0(matcmd, cmd)
  cmd = code
  fname = tempfile(fileext = ".m")
  cat(cmd, file = fname)
  if (verbose) {
    message(paste0("Script created: ", fname))
  }
  x = run_matlab_script(fname, verbose = verbose, ...)
  return(x)
}




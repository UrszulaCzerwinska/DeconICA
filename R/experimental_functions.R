####TO DO - add system dependent behaviour for path
#windows cd 'C:\Program Files\MATLAB\R2017a\examples\matlab_featured'
#to call matlab : path to matlab exe -nosplash -r mfile -logfile c:\temp\logfile


.doICA <-
  function(df.scaled.t,
           names,
           samples,
           path_global,
           n,
           name = FALSE,
           corr_folder = "CORRELATION",
           matlbpth = "/Applications/MATLAB_R2016a.app/Contents/MacOS/MATLAB_maci64",
           fasticapth = "/Users/ulala/Documents/CURIE/BIODICA-master/BIODICA/bin/fastica++") {
    setwd(path_global)
    if (name != FALSE) {
      name <- name
    } else {
      name <- deparse(substitute(df.scaled.t))
    }

    dir.create(paste(name, n, sep = "_"))

    path_global_1 = paste(path_global, "/", paste(name, n, sep = "_"), "/", sep =
                            "")
    setwd(path_global_1)

    utils::write.table(
      df.scaled.t,
      file = paste(name, "_numerical.txt", sep = ""),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    utils::write.table(
      names,
      file = paste(name, "_ids.txt", sep = ""),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )

    utils::write.table(
      samples,
      file = paste(name, "_samples.txt", sep = ""),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )

    fun = "doICA"

    ncomp = n

    path = paste("'", path_global_1, "'", sep = "")

    file = paste("'", paste(name, "_numerical.txt", sep = ""), "'", sep =
                   "")
    matlbpth = "/Applications/MATLAB_R2016a.app/Contents/MacOS/MATLAB_maci64"
    fasticapth = "/Users/ulala/Documents/CURIE/BIODICA-master/BIODICA/bin/fastica++"
    start = paste(matlbpth, " -nodisplay -r \"cd('", fasticapth, "');", sep =
                    "")
    cmd = paste(start, " ", fun, "(", path, ",", file, ",", ncomp, "); quit\"", sep =
                  "")

    system(cmd)
    a.imp = paste(paste("A", paste(name, "_numerical.txt", sep = ""), ncomp, sep =
                          "_"), "num", sep = ".")
    a.imp.path = paste(path_global_1, a.imp, sep = "")
    A = utils::read.delim(a.imp.path, header = FALSE, sep = "\t")
    s.imp = paste(paste("S", paste(name, "_numerical.txt", sep = ""), ncomp, sep =
                          "_"), "num", sep = ".")
    s.imp.path = paste(path_global_1, s.imp, sep = "")
    S = utils::read.delim(s.imp.path, header = FALSE, sep = "\t")

    ics = paste("IC", c(1:ncomp), sep = "")
    colnames(S) = ics
    colnames(A) = ics

    A = A[, 1:(ncol(A) - 1)]
    S = S[, 1:(ncol(S) - 1)]

    path1 <- paste("../", corr_folder, "/", sep = "")
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
    utils:: write.table(
      data.frame(SAMPLES = samples, A),
      file = paste(name, ncomp, "full_samples", "txt", sep = ".") ,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
    utils::write.table(
      data.frame(GENES = names, cbind(S)),
      file = paste(name, ncomp, "full_genes", "txt", sep = ".") ,
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
    return(list(A = t(A), S = S))
  }


.export_for_ICA <-
  function(df.scaled.t,
           names,
           samples,
           path_global = getwd(),
           name = FALSE) {
    setwd(path_global)
    if (name != FALSE) {
      name <- name
    } else {
      name <- deparse(substitute(df.scaled.t))
    }

    dir.create(name)

    path_global_1 = paste(path_global, "/", name, "/", sep =
                            "")
    setwd(path_global_1)

    utils::write.table(
      df.scaled.t,
      file = paste(name, "_numerical.txt", sep = ""),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    utils::write.table(
      names,
      file = paste(name, "_ids.txt", sep = ""),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )

    utils::write.table(
      samples,
      file = paste(name, "_samples.txt", sep = ""),
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )

  }


.prepare_data_for_ica <-function(df, names) {
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

  # center rows if needed
  X <- .center_rowmeans(X)
  if (all.equal(round(mean(as.matrix(X[1, ])), 2), 0) != TRUE)
    message("your data is not mean centerd by rows (genes)")
  return(list(df.scaled = X, non.scaled = X.ndup ,names = as.character(names)))
}

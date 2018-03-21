#' Test impact of each Independ Component
#'
#' This function is applying distribution statistical test (i.e. \code{t.test},
#' \code{wilcox.test}) to evaluate which ICs have highest impact on differences
#' between samples
#'
#' @param A result of \code{\link{run_fastica}} the \code{A} matrix
#' @param sample.names names of samples, should correspond to number of columns
#'   of \code{A}
#' @param quant quantiles to use, in form of \code{c(x, y)}
#' @param test.type test of distributions to perform
#' @param isLog by default \code{NULL}, if \code{X} is not conts but log,
#'   provide the base of log, for natural logarithm use \code{exp(1)}
#' @param X.counts expression data
#' @param return if you want to return p.values select
#' @param wide should the output matrix be in wide format (FALSE preferable for
#'   plotting)
#' @param thr threshold of maximal p.value consiered 0.1 by default
#'
#' @return returns a matrix (in long or wide) format
#' @export
#'
#' @examples
#'
#'
dist_test_samples <-
  function(A,
           sample.names,
           quant = c(0.1, 0.9),
           X.counts,
           test.type,
           thr = 0.1,
           isLog = NULL,
           return = "p.value",
           wide = TRUE) {
    if (nrow(A) > 30)
      message(
        "this procedure is time consuming\n...............\nin average it
        takes 5min - 2h30\n...............\nbe patient"
      )
    df <- as.matrix(t(A))
    row.names(df) <- sample.names
    colnames(X.counts) <- sample.names
    if (!is.null(isLog)) {
      X.counts <- isLog ^ X.counts - 1
    }
    fun <-
      function(sel.col,
               quant,
               X.counts,
               test.type,
               return,
               ...) {
        #
        #
        qt <- stats::quantile(sel.col, quant)
        group1 <- which(sel.col < qt[1])
        group2 <- which(sel.col > qt[2])

        exprs.group1 <- X.counts[, names(group1)]
        exprs.group2 <- X.counts[, names(group2)]
        mat <- cbind(exprs.group1, exprs.group2)
        if (test.type == "exactTest") {
          mat <- data.frame(mat)
          DGE <-
            edgeR::DGEList(mat, group = rep(c("Cond1type", "Cond2type"),
                                            c(length(group1), length(group2))))
          # Analysis using common dispersion
          disp <-
            edgeR::estimateCommonDisp(DGE) # Estimating the common dispersion
          #tested <- exactTest(disp,pair=c("Normal","Tumor")) # Testing
          tested <-
            edgeR::exactTest(disp, pair = c("Cond1type", "Cond2type")) # Testing
          return(tested$table$PValue[order(tested$table$PValue)])
        } else {
          if (test.type == "t.test") {
            res <- apply(mat, 1, function(x)
              stats::t.test(x[1:length(group1)], x[length(group1):length(x)]))
          }
          else {
            FUN <- match.fun(test.type)
            res <- apply(mat, 1, function(x)
              FUN(x[1:length(group1)], x[length(group1):length(x)]), ...)
          }
          p.value <- sapply(res, function(x)
            x[[3]])
          stat <- sapply(res, function(x)
            x[[1]])
          res.df <- data.frame(p.value, stat)
          row.names(res.df) <- row.names(mat)
          switch(return,
                 p.value = return(sort(
                   p.value, decreasing = FALSE, na.last = TRUE
                 )),
                 stat =  return(sort(
                   stat, decreasing = FALSE, na.last = TRUE
                 )),
                 both =  return(res.df))
        }

      }

    res.df <-
      apply(
        X = df,
        MARGIN = 2,
        FUN = fun,
        quant = quant,
        X.counts = X.counts,
        test.type = test.type,
        return = return
      )


    if (return == "both") {
      return(res.df)
    } else {
      if (!is.null(thr)) {
        keep <- apply(res.df, 1, function(x)
          all(x < thr))
        res.df <- stats::na.omit(res.df[keep, ])
      }
      if (wide) {
        return(res.df)
      } else {
        return(reshape::melt(data.frame(id = as.numeric(
          1:nrow(res.df)
        ), res.df), id = "id"))
      }
    }


  }

# res.ttest <- dist_test_samples(A = ica.matlab.BRCATCGA$A[1:2, 1:30],
#                                sample.names = ica.matlab.BRCATCGA$samples[1:30]  ,
#                                quant = c(0.1,0.9),
#                                X.counts = ica.matlab.BRCATCGA$log.counts[1:30],
#                                test.type = "t.test",
#                                isLog = 2,
#                                return = "p.value",
#                                wide = FALSE)
# ###### PLOT
#  ggplot2::ggplot(res.ttest , ggplot2::aes(x = value, color = variable)  ) +
#    ggplot2::stat_density(position = "identity", geom = "line") + ggplot2::theme_bw()
#
#  ggplot2::ggplot(data=res.ttest,  ggplot2::aes(y=value, x=id, color = variable)) +
#    ggplot2::geom_line()
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Compute variance explained by each Independent Component
#'
#' @param S result of \code{\link{run_fastica}} the \code{S} matrix
#' @param A result of \code{\link{run_fastica}} the \code{A} matrix
#' @param X data, either post-PCA data of \code{\link{run_fastica}} \code{X} matrix
#' @param n number of top ICs if \code{n} = "all" then fev for all ICs is returned
#'
#' @return
#' returns a data frame with n top ICs numbers ranked by their fraction of variance explained
#' @export
#'
#' @examples
#'
#'
most_variant_IC <- function(S, A, X, n = 5) {
  message("computing fraction of variance explained by all ICs")
  S <- as.matrix(S)
  A <- as.matrix(A)
  X <- as.matrix(X)
  A.vec <- t(apply(A, 1, .norm_vec))
  S.norm <- S %*% diag(array(A.vec))
  fev <- as.data.frame(.colVars(S.norm) / sum(.colVars(X)))
  row.names(fev) <- paste("IC", 1:ncol(S), sep = "")
  if (n == "all") {
    res <- fev[order(-fev), , drop = FALSE]
  } else {
    res <- fev[order(-fev)[1:n], , drop = FALSE]
  }
  colnames(res) <- "fev"
  return(res)
}

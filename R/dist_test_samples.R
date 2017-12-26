#' Title
#'
#' @param assign result of \code{\link{assign_metagenes}}
#' @param A result of \code{\link{run_fastica}} the \code{A} matrix
#' @param sample.names names of samples, should correspond to number of columns
#'   of \code{A}
#' @param quant quantiles to use, in form of \code{c(x, y)}
#' @param normalized.count.df matrix of counts used for \code{\link{run_fastica}} if
#'   data in log provided, it will be un-log to counts, provide \code{isLog},
#'   best to use \code{raw.counts} or \code{log.counts} from  \code{\link{run_fastica}} to avoid problem
#'   with duplicated gene names
#' @param test.type test of distributions to perform
#' @param wide if \code{TRUE} returns n by n.ics matrix, if \code{FALSE} returns
#'   long format convininient for \code{ggplot}
#' @param isLog by default \code{NULL}, if \code{normalized.count.df} is not
#'   conts but log, provide the base of log, for natural logarithm use
#'   \code{exp(1)}
#'
#' @return returns a matrix (in long or wide) format
#' @export
#'
#' @examples
#'
#'
dist_test_samples <- function(assign, A, sample.names, quant, normalized.count.df,test.type, wide = TRUE, isLog = NULL) {
  important.ics <- as.numeric(gsub("IC", "",assign[,2]))
  df <-t(A[important.ics, ])
  colnames(df) <- assign [,1]
  row.names(df) <- sample.names
  if (is.null(!isLog)) {
    normalized.count.df <- isLog^normalized.count.df - 1
  }
  fun <- function(sel.col,quant, normalized.count.df,test.type) {


    qt<- stats::quantile(sel.col, quant)
    group1 <- which(sel.col < qt[1])
    group2 <- which(sel.col > qt[2])

    exprs.group1 <- normalized.count.df[,names(group1)]
    exprs.group2 <- normalized.count.df[,names(group2)]

    mat <- data.frame(exprs.group1, exprs.group2 )


    if (test.type == "exactTest") {

      DGE <- edgeR::DGEList(mat,group = rep(c("Cond1type","Cond2type"),
                                            c(length(group1),length(group2))))
      # Analysis using common dispersion
      disp <- edgeR::estimateCommonDisp(DGE) # Estimating the common dispersion
      #tested <- exactTest(disp,pair=c("Normal","Tumor")) # Testing
      tested <- edgeR::exactTest(disp,pair = c("Cond1type","Cond2type")) # Testing
      return(tested$table$PValue[order(tested$table$PValue)])
    } else {
      test <- .compute_x_test_for_many(mat, names(group1), names(group2), test = stats::test.type)
      return(test$pvalue[order(test$pvalue)])
    }

  }
  #res.df <- data.frame(matrix(rep(NA, ncol(df)*nrow(normalized.count.df)), ncol = ncol(df), nrow= nrow(normalized.count.df) ))
  res.df <- apply(df[1:100,], 2, fun, quant = quant, normalized.count.df = normalized.count.df, test.type = test.type )
  if (wide) {
    return(res.df)
  } else {
    return(reshape::melt(data.frame(id = as.numeric(1:nrow(res.df)), res.df), id = "id" ))
      }
}


# res.ttest <- dist_test_samples(assign =  BRCAMatrixTP.assign,
#                                ica = BRCAMatrixTP.ica,
#                                sample.names =   colnames(BRCAMatrixTP.names) [2:ncol(BRCAMatrixTP.names)],
#                                quant = c(0.1,0.9),
#                                normalized.count.df = BRCAMatrixTP.red,
#                                test.type = "exactTest")
#

##### PLOT
# ggplot2::ggplot(res.df.long, ggplot2::aes(x = value, color = variable)  ) +
#   ggplot2::stat_density(position = "identity", geom = "line") + ggplot2::theme_bw()

#' @importFrom magrittr "%>%"
#'
#'
#'

.center_rowmeans <- function(x) {
  xcenter <- rowMeans(x)
  x - rep(xcenter, ncol(x))
}

.orient_funct <- function(S) {
  orient <-
    apply(S, 2, function(x)
      ifelse (sum(x > 3)  < sum(x < -3), -1, 1))
  S <- as.matrix(S)  %*% diag(orient)
  return(S)
}

.rowVars <- function (x, na.rm = TRUE) {
  x <- as.matrix(x)
  sqr <- function(x)
    x * x
  n <- rowSums(!is.na(x))
  n[n <= 1] <- NA
  rowSums(sqr(x - rowMeans(x, na.rm = na.rm)), na.rm = na.rm) / (n - 1)
}

.cumVar <- function(pca) {
  cumsum(pca$sdev ^ 2 / sum(pca$sdev ^ 2))
}

.norm_vec <- function(x)
  sqrt(sum(x ^ 2))

.remove_duplicates <- function(X) {
  message("for duplicated genes, entry with higher variance will be kept")
  X[["var"]] <- .rowVars(X[, 2:ncol(X)])
  X.o <-  X[order(X[, 1], -X[["var"]]), ]
  X <- X.o[!duplicated(X.o[, 1]), ]
  X[["var"]] <- NULL
  return (X)
}

.intersect.genes <-
  function(res, component)
    intersect(res[, 1], component[, 1])

# y = metagenes[[1]]
# res = res
# name = names(metagenes)[1]
.corr_matrix <- function (y, res, name) {
  res[[name]] <- NULL
  res[match(.intersect.genes(res, y), res[, 1]), name] <-
    y[match(.intersect.genes(res, y), y[, 1]), 2]
  return(res)
}

.verify.n <-
  function(n, n.genes.intersect)
    apply(n, 1, function(ic)
      any(ic > n.genes.intersect))

.fev <- function(S,A,X, i) {
  A.vec <- t(apply(A ,1, .norm_vec ))
  S.norm <- S %*% diag(array(A.vec))
  print(var(S.norm[,i]) / sum(colVars(X)) *100, "%", sep= " ")
}

# df = BRCAMatrixTP.red
#  array1 = group1
#  array2 = group2
.compute_x_test_for_many <- function (df, array1, array2, test = "t.test") {

  pvalue = NULL # Empty list for the p-values
  stat = NULL # Empty list of the t test statistics

  for(i in 1 : nrow(df)) { # For each gene :

    x = df[i,array1] # condition1 of gene number i
    y = df[i,array2] # condition2 of gene number i

    # Compute a test between the two conditions
    x = do.call(test,list(as.matrix(x), as.matrix(y)))

    # Put the current p-value in the pvalues list
    pvalue[i] = x$p.value

    # Put the current test statistic in the tstats list
    stat[i] = x$statistic

  }

  pvalue.df=data.frame(pvalue)
  stat.df=data.frame(stat)

  row.names(pvalue.df)=row.names(df)
  row.names(stat.df)=row.names(df)

  return(cbind(pvalue.df,stat.df))
}


#####-------------- WORKING PROGRESS
#' Computes correlation between ICs and an know rank (created for immune component)
#'
#' @param res S matrix of ICA, 1st column are gene names, 2:ncol(res) are ICs
#' @param M8_IMMUNE immune component or other vector to correlate with
#' @param type possible options :
#'
#'             "no_thr" - no threshold
#'
#'             "immune_threshold" - threshold >2 on immune (ranks vector)
#'
#'             "IC_threshold" threshold >3 on ICs
#'
#'             "immune_IC_threshold" both
#' @param plot if true plot are saved in working repo
#' @param name eee
#' @param S.thr eee
#' @param forCor r"r"
#'
#' @return returns r2 of correlations between all ICs and the rank vector
#'
#' @examples
#'
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- data.frame(S %*% A)
#' res <- run_fastica(X, center = TRUE, n.comp = 2)
#' data(M8_IMMUNE)
#' res.orient <- .orient_funct(res$S)
#' .scatter_correl(res.orient, M8_IMMUNE, type = "no_thr", plot = TRUE)

.scatter_correl <-
  function(res,
           M8_IMMUNE,
           type = "no_thr",
           name = "IMMUNE",
           plot = TRUE,
           S.thr = FALSE,
           forCor = FALSE) {
    #
    # res <- res_tcga.2
    # M8_IMMUNE
    # type = "no_thr"
    # plot = FALSE
    # apply if type adequate
    if (type == "immune_threshold" |
        type == "immune_IC_threshold") {
      M8_IMMUNE <-
        M8_IMMUNE %>% dplyr::filter(M8_IMMUNE[, 2] > 3 * stats::sd(M8_IMMUNE[, 2]))
    }
    if (type == "IC_threshold" |
        type == "immune_IC_threshold") {
      res.thr <- res[, 2:ncol(res)]
      res.thr [res.thr < 3] <- NA #all values less than 3
      res <- cbind(res[, 1], res.thr)
    }
    #verify interesection for each component with immune
    intersect.ics <-
      apply(res[, 2:ncol(res)], 2, function(x)
        length(intersect(res[!is.na(x), 1], M8_IMMUNE[, 1])))
    #delete component if not enough genes to correlate
    if (length(which(intersect.ics <= 30)) == (ncol(res) - 1)) {
      stop("no IC to correlate")
    }
    if (length(which(intersect.ics <= 30)) > 0) {
      res <-
        res[, -c(which(intersect.ics <= 30) + 1)]
    }
    #work only which genes that are common to IMMUNE and ICA matrix
    intersect.vec  <-
      intersect(res[, 1], M8_IMMUNE[, 1])
    #sanity check
    condition <- !(identical(as.matrix(res[match(intersect.vec, res[, 1]), 1]) ,
                             as.matrix(M8_IMMUNE[match(intersect.vec, M8_IMMUNE[, 1]), 1])))
    if (condition) {
      stop ("problem with gene names intersection")
    }
    if (S.thr &
        type == "immune_IC_threshold" | S.thr & type == "IC_threshold") {
      colnames(res)[1] <- "GENES"
      return(res)
    }
    # in case "IMMUNE" column exist delete
    res[[name]] <- NULL
    res[match(intersect.vec, res[, 1]), name] <-
      M8_IMMUNE[match(intersect.vec, M8_IMMUNE[, 1]), 2]
    if (forCor) {
      return(as.matrix(res[, 2:ncol(res)]))
    }
    #make correlatons one to one, NA are ommited
    rcorr.res <-
      Hmisc::rcorr(as.matrix(res[, 2:ncol(res)]))
    #keep only results for correlation with IMMUNE
    df <-
      data.frame(name = rcorr.res$r[, c(name)])
    colnames(df)[which(colnames(df) == "name")] <- name
    #keep only correlations >0.1
    res.fil <-
      df %>% tibble::rownames_to_column() %>% dplyr::filter(eval(parse(text = name)) > 0.1)
    # save scatter plots
    if (plot) {
      for (i in 1:(nrow(res.fil) - 1)) {
        grDevices::jpeg(paste(type, paste(as.matrix(res.fil[i, 1])), "jpg", sep = "."))
        plot(
          res[[name]],
          res[[paste(as.matrix(res.fil[i, 1]))]],
          main = paste(type, " r2 = ", res.fil[[i, 2]], sep = ""),
          ylab = paste(as.matrix(res.fil[i, 1]))
        )
        grDevices::dev.off()
      }
    }
    # return all correlations with IMMUNE
    return(df)
  }

.getCorr <- function(X, thr = 0.2) {
  res <-
    X[, "IMMUNE", drop = FALSE] %>% data.frame %>% tibble::rownames_to_column() %>% filter(IMMUNE > thr)  %>% arrange(-IMMUNE) %>% filter(row_number() > 1)
  return(res)
}


.resampling.correl <- function(X, genes, n) {
  df <- list()
  del.genes.list <- list()
  ngenes <- sum(is.na(X[["IMMUNE"]]))
  if (ngenes == 0L) {
    ngenes <- length(genes)
  }
  indexes.non.na <- which(!is.na(X[["IMMUNE"]]))
  #sanity check
  if (!magrittr::equals(length(genes), nrow(X)))
    stop("problem with genes")
  for (i in 1:n) {
    #for (i in 1:2) {
    set.seed(i + 100)
    delete.index <- sample(indexes.non.na, 0.1 * ngenes)
    #deleted genes
    del_genes <- genes[delete.index]
    del.genes.list[[i]] <- del_genes
    X.del <- X[-delete.index,]
    rcorr.res <-
      Hmisc::rcorr(as.matrix(X.del))
    df[[i]] <-  rcorr.res$r[, c("IMMUNE")]
  }
  myres <-
    list(corel.simul = data.frame(df), del.genes.list = del.genes.list)
  return(myres)
}

#X - df with n simulatons in colums and IC in rows
.median.corr <- function(X) {
  data.frame(corr = X %>%  apply(1, median)) %>% tibble::rownames_to_column() %>% arrange(-corr) %>% filter(row_number() > 1)
}

#from sampling correlations
.getOutliers <-
  function(X,
           genes.list,
           icn,
           percent_simul = 0.05,
           percent_genes = 0.05,
           plot_top = 10) {
    ic <- as.matrix(X)[icn, ]
    q <- quantile(ic, c(percent_simul, 1 - percent_simul))
    out <- which(ic > q[2] | ic < q[1])
    frq.tab <- data.frame (table(unlist(genes.list[out])))
    q.g <- quantile(frq.tab[["Freq"]], c(1 - percent_genes))
    res <-
      data.frame(frq.tab %>%  arrange(-Freq) %>% filter(Freq > q.g) %>% filter(row_number() <= plot_top))
    return(res)
  }


.export_top_genes <- function(X, ic, n, name) {
  topgenes <-
    X %>% arrange(-(eval(parse(text = ic)))) %>% select(1) %>% filter(row_number() <=
                                                                        n) #tcell
  write.table(
    topgenes,
    paste(name, "_", ic, ".txt", sep = ""),
    quote = FALSE,
    sep = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
}

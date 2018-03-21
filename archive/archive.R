#' @title Run enrichment analysis for metagenes
#'
#' @description Verifies enrichment in provided signatures for selected
#'   metagenes
#'
#' @details \code{get_enrichment} runs enrichment of an independent component
#'   (or any ranked list) in known immune cell types signatures. By default it
#'   is using \code{S} matrix from ICA \code{\link{run_fastica}} and runs
#'   enrichment using \code{\link[ACSNMineR]{enrichment}} only on components
#'   indentified as correlated with immune metagene through function
#'   \code{identify_immune_ic} and it serches in Immgen signatures
#'   \url{http://Immgen.org}.
#'
#' @param S ranks of components, dim \code{n} corresponding to genes, \code{m}
#'   corresponding to number of components
#' @param gene.names list of gene names, length needs to be equal to \code{n}
#' @param immune.ics list of components to use for enrichment test
#' @param gmt list of signatures, by default Immgen signatures
#'   \url{http://Immgen.org}
#' @param quantile takes numerical value \code{[0,1]}, default \code{= 0.99}, to
#'   idicate quantile of ranked list that should be used of enrichment, if
#'   \code{ = NULL}, \code{abs.thr} needs to be indicated
#' @param abs.thr impose threshold for ranked list, specify  \code{quantile =
#'   NULL} if you want to use \code{abs.thr}
#' @param ... other parameters you want to pass to
#'   \code{\link[ACSNMineR]{enrichment}}
#' @inheritParams ACSNMineR::enrichment
#' @return Function returns value if there is an enrichment in provided
#'   signatures: \describe{ \item{metagenes}{interpreted metagene gene ranking}
#'   \item{enrich}{full results of the enrichment analysis sorted by corrected
#'   p.value} }
#' @seealso \code{\link{identify_immune_ic}} identificng immune related
#'   components, \code{\link{run_fastica}} for running Independent Components
#'   Analysis, and \code{\link[ACSNMineR]{enrichment}} for enrichment in gmt
#'   files
#' @export
#'
#' @examples
#'
#' # run enrichment
#'

# S = BRCAMatrixTP.ica$S
# gene.names = BRCAMatrixTP.ica$names
# immune.ics = names(BRCAMatrixTP.identify)
# gmt = ImmgenHUGO
# quantile = 0.995
# abs.thr = NULL
# threshold = 0.05
get_enrichment <-
  function(S,
           gene.names,
           immune.ics,
           gmt = ImmgenHUGO,
           quantile = 0.995,
           abs.thr = NULL,
           threshold = 0.05,
           ...) {
    genes.list <- list()
    results.list <- list()
    metgenes.list <- list()
    # orient components
    S <- .orient_funct(S)
    if (length(as.matrix(gene.names)) != nrow(S))
      stop("wrong number of gene names")
    # add colnames
    gene.names <- as.data.frame(gene.names)
    colnames(S) <- paste("IC", 1:ncol(S), sep = "")
    colnames(gene.names) <- "gene.names"
    table <- data.frame(gene.names, S)[c("gene.names", immune.ics)]
    message("running enrichment tests \n can take some time...")
    # for each possible immune related component
    metgenes.list <- sapply(immune.ics, function(i) {
      # apply quantile or direct threshold
      if (!(is.null(quantile)))
        abs.thr <- quantile(table[, i], quantile)
      t <- table[which(table[, i] > abs.thr), c("gene.names", i)]
      t <- t[order(-t[,2]),]
    }, simplify = FALSE, USE.NAMES = TRUE)

    genes.list <-
      sapply(metgenes.list, function(t)
        names <- as.array(t[, 1]), simplify = FALSE, USE.NAMES = TRUE)


    results.list <-
      ACSNMineR::multisample_enrichment(
        genes.list,
        threshold = threshold,
        maps = gmt,
        cohort_threshold = FALSE,
        ...
      )

    results.list <-
      lapply(results.list, function(Example)
        suppressWarnings(if (!(is.na(Example)))
          Example[order(Example[["p.value.corrected"]]),]))
    results.list <-results.list[sapply(results.list, function(s)
      ! (is.null(s)))]
    metgenes.list <-
      metgenes.list[sapply(results.list, function(s)
        ! (is.null(s)))]
    return(list(enrich = results.list ,
                metagenes = metgenes.list))
  }


.export_top_genes <- function(X, ic, n, name) {
  topgenes <-
    X %>% arrange(-(eval(parse(text = ic)))) %>% select(1) %>% filter(row_number() <= n) #tcell
  write.table(
    topgenes,
    paste(name, "_", ic, ".txt", sep = ""),
    quote = FALSE,
    sep = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
}

#################
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
#
#X - df with n simulatons in colums and IC in rows
.median.corr <- function(X) {
  data.frame(corr = X %>%  apply(1, median)) %>% tibble::rownames_to_column() %>% arrange(-corr) %>% filter(row_number() > 1)
}
################
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


.getCorr <- function(X, thr = 0.2) {
  res <-
    X[, "IMMUNE", drop = FALSE] %>% data.frame %>% tibble::rownames_to_column() %>% filter(IMMUNE > thr)  %>% arrange(-IMMUNE) %>% filter(row_number() > 1)
  return(res)
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


#'Correlate components with known ranked lists of genes
#'
#'Components obtained, for example, with \code{\link{run_fastica}} can be
#'characterized through correlation with known ranked list (metagenes or profiles), by
#'default this function is using metagenes from Biton et al. (2015), Cell. It is
#'using \code{\link[Hmisc:rcorr]{rcorr}} function for correlations
#'
#'@param S S matrix of components
#'@param gene.names list of gene names, needs to be of the same length as nrow
#'  of \code{S}, for ICA it is recommended to run \code{\link{run_fastica}}
#'  \code{with.names = TRUE} to assure compatibility
#'@param metagenes named list of datasets, each with two columns 1st - gene
#'  names, 2nd - ranks, by default 11 metagenes from Biton et al. (2015), Cell
#'@param threshold threshold for components (columns of \code{S}) to be
#'  applied before correlation, default set to -Inf (all ranks are kept)
#'@param n.genes.intersect minimum of genes that should intersect between a component
#'and a metagene to keep the component in correlation matrix
#'@param orient.long orient by long tails, default TRUE
#'@param orient.max orient by maximal correlation, default FALSE, can be used if
#'  there is no long tails
#'@param ... additional params you can pass to \code{\link[Hmisc:rcorr]{rcorr}}
#'@inheritParams Hmisc::rcorr
#'@return a correlation matrix with correlation coefficient \code{r},
#'  p.values \code{P} and number of overlapping genes \code{n}, oriented
#'  \code{S} matrix
#'
#'@export
#'
#'@seealso \code{\link[Hmisc]{rcorr}} \code{\link{run_fastica}} \code{\link{make_list}}
#'
#'@examples
#'res_run_ica <- run_fastica (
#'  Example_ds,
#'  overdecompose = FALSE,
#'  n.comp = 5,
#'  with.names = TRUE
#')
#'correlate_metagenes(
#'    S = res_run_ica$S,
#'    gene.names = res_run_ica$names)
#'
correlate_metagenes <-
  function(S,
           gene.names,
           metagenes = Biton.list,
           threshold = -Inf,
           n.genes.intersect = 30,
           orient.long = TRUE,
           orient.max = FALSE,
           ...) {
    # orient components in the direction of the long tail
    if (orient.long &
        orient.max)
      stop("select one or none orienting method")
    if (orient.long) {
      S_or <- S <- .orient_funct(S)
      colnames(S_or) <- colnames(S) <- paste0("IC", 1:ncol(S))
    } else {
      colnames(S) <- paste0("IC", 1:ncol(S))
    }

    # verify if gene names are correct
    if (length(as.matrix(gene.names)) != nrow(S))
      stop("wrong number of gene names")
    # add colnames
    # colnames(S) <- paste("IC", 1:ncol(S), sep = "")
    gene.names <- as.data.frame(gene.names)
    colnames(gene.names) <- "gene.names"
    # all values less than "threshold"
    if (!orient.max)
      S[S < threshold] <- NA
    res <- data.frame(gene.names, S)
    for (i in 1:length(metagenes)) {
      res <- .corr_matrix(metagenes[[i]], res, names(metagenes)[i])
    }
    rcorr.res <-
      Hmisc::rcorr(as.matrix(res[, 2:ncol(res)]), ...)
    n <-
      rcorr.res$n[1:(nrow(rcorr.res$r) - length(metagenes)), names(metagenes)]
    verify.n <- .verify.n(n, n.genes.intersect)
    if (length(which(verify.n == FALSE)) > 0)
      message(paste("deleting ", names(which(!verify.n)), "\n", sep = ""))
    n <- rcorr.res$n[which(verify.n), names(metagenes)]
    r <- rcorr.res$r[which(verify.n), names(metagenes)]
    p <- rcorr.res$P[which(verify.n), names(metagenes)]

    if (orient.max) {
      S <- .orient_max(S, r)
      colnames(S) <- paste0("IC", 1:ncol(S))
      S[S < threshold] <- NA
      res <- data.frame(gene.names, S)
      for (i in 1:length(metagenes)) {
        res <- .corr_matrix(metagenes[[i]], res, names(metagenes)[i])
      }
      rcorr.res <-
        Hmisc::rcorr(as.matrix(res[, 2:ncol(res)]), ...)
      n <-
        rcorr.res$n[1:(nrow(rcorr.res$r) - length(metagenes)), names(metagenes)]
      verify.n <- .verify.n(n, n.genes.intersect)
      if (length(which(verify.n == FALSE)) > 0)
        message(paste("deleting ", names(which(!verify.n)), "\n", sep = ""))
      n <- rcorr.res$n[which(verify.n), names(metagenes)]
      r <- rcorr.res$r[which(verify.n), names(metagenes)]
      p <- rcorr.res$P[which(verify.n), names(metagenes)]
      return(list(
        S.max = S,
        n = n,
        r = r,
        P = p
      ))
    } else{
      return(list(
        S.or =  S_or,
        n = n,
        r = r,
        P = p
      ))
    }
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Assign components to a metagene through mutual reciprocity
#'
#' Attributes labels to components under condition of mutual reciprocal correlation
#'
#' @details  This function assign a component to a metagene/profile through
#' verification if the component's the maximal correlation points to a given profile and if
#' for this profile the maximal correlation points back the that component. In mathematical
#' terms, given correlations between the set of profiles/metagenes \eqn{A = {A_1,...,A_m}} and
#' \eqn{S} components matrix \eqn{S = {IC1,...,ICN}}, if
#' \deqn{Si = argmaxi(corr(Aj,S))} and \deqn{A_j = argmax_j(corr(S_i,A))}
#'
#' @param r the correlation matrix, \code{r} matrix, can be generated from
#' \code{\link{correlate_metagenes}} function
#' @param exclude_name name of the components (present in \code{r}) to
#'   be excluded from this analysis (for example immune), by default "M8_IMMUNE"
#'   is excluded
#'
#' @return returns a \code{dataf.rame} with component name in the first
#' column and assigned profile/metagene name in second column
#'
#' @seealso \code{\link{get_max_correlations}}, \code{\link{correlate_metagenes}}
#'
#' @export
#'
#' @examples
#'
#' res_run_ica <- run_fastica (
#'  Example_ds,
#'  overdecompose = FALSE,
#'  n.comp = 5,
#'  with.names = TRUE
#')
#'corr <- correlate_metagenes(
#'    S = res_run_ica$S,
#'    gene.names = res_run_ica$names)
#'
#'assign_metagenes(corr$r)
#'
#'
assign_metagenes <- function(r, exclude_name = "M8_IMMUNE") {
  r <- data.frame(r)
  if (!is.null(exclude_name)) {
    r <- r[, -c(match(exclude_name, colnames(r)))]
    message("profiles excluded")
  } else
    message("no profiles to exlude provided")
  col.as <- apply(r, 2, function(ic)
    which(ic == max(ic)))
  row.as <- apply(r, 1, function(meta)
    which(meta == max(meta)))
  l <- list()
  for (i in 1:length(col.as))
    for (j in 1:length(row.as))
      if (col.as[i] == j &
          row.as[j] == i)
        l[[i]] <- c(names(col.as)[i], names(row.as)[j])
  # creating data frame
  df <-
    data.frame(profile = unlist(sapply(l, "[[", 1)),
               component = unlist(sapply(l, "[[", 2)))
  message("DONE")
  return(df)
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Identify components related to immune signal
#'
#' @param x the correlation with immune metagene can be retrieved from
#'   \code{\link{correlate_metagenes}} output
#' @param l vector of names of assigned components
#' @param threshold lower bound for filtering correlation [0,1]
#'
#' @return
#' it returns data frame of component names and correlations passing the
#' \code{threshold}
#' @export
#'
#' @examples
#'
#' res_run_ica <- run_fastica (
#'  Example_ds,
#'  overdecompose = FALSE,
#'  n.comp = 20,
#'  with.names = TRUE
#')
#'corr <- correlate_metagenes(
#'    S = res_run_ica$S,
#'    gene.names = res_run_ica$names)
#'
#'assign <- assign_metagenes(corr$r)
#'
#'identify_immune_comp(corr$r[,"M8_IMMUNE"], assign[, "component"], threshold = 0.1)
identify_immune_comp <- function(x, l, threshold = 0.1) {
  x[which(x > threshold)][!(names(which(x > threshold)) %in% l)]
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Assign through maximal correlations
#'
#' It assigns maximal correlations between set of correlated vectors
#'
#' @param corr list of correlation matrices with correlation coefficients
#' and p-values, can be obtained from \code{\link{correlate_metagenes}} or
#' \code{\link[Hmisc:rcorr]{rcorr}}
#'
#' @return \code{data.frame} with matched column names,
#' Pearson correlation coefficient, p.value
#'
#' @export
#'
#' @seealso \code{\link[Hmisc:rcorr]{rcorr}}, \code{\link{correlate_metagenes}},
#' \code{\link{assign_metagenes}}
#'
#' @examples
#'res_run_ica <- run_fastica (
#'  Example_ds,
#'  overdecompose = FALSE,
#'  n.comp = 5,
#'  with.names = TRUE
#')
#'corr <- correlate_metagenes(
#'    S = res_run_ica$S,
#'    gene.names = res_run_ica$names)
#'
#'get_max_correlations(corr)
#'
get_max_correlations <- function(corr) {
  r <- data.frame(corr$r)
  p <- data.frame(corr$P)
  col.as <-
    apply(r, 2, function(ic)
      cbind(which(ic == max(ic)), max(ic)))
  p.val <- sapply (1:ncol(p), function(i)
    p[col.as[1, i], i])
  data.frame(
    TYPE = colnames(col.as),
    IC = row.names(r)[col.as[1, ]],
    r = col.as[2, ],
    p.val = p.val
  )
}

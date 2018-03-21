#stopifnot(require(deconica, quietly=TRUE))

#'Correlate independent components with known ranked lists of metagenes
#'
#'Independent components obtained with \code{\link{run_fastica}} can be
#'characterised through correlation with known ranked list (metagenes), by
#'default this function is using metagenes from Biton et al. (2015) Cell It is
#'using \code{\link[Hmisc]{rcorr}} function for correlations
#'
#'@param S S matrix of ICA decomposition
#'@param gene.names list of gene names, needs to be of the same length as nrow
#'  of S, it is recommended to run \code{\link{run_fastica}} \code{with.names =
#'  TRUE} to assure compatibility
#'@param metagenes named list of datasets, each with two columns 1st - gene
#'  names, 2nd - ranks, by default 11 metagenes from Biton et al....
#'@param threshold threshold for ICA components (columns of \code{S}) to be
#'  applied before correlation, default set to -Inf (all ranks are kept)
#'@param n.genes.intersect minimum of genes that should intersect between IC and
#'  metagene to keep the IC in correlation matrix
#'@param orient.long orient by long tails,  default TRUE
#'@param orient.max orient by maximal correlation, default FALSE, can be used if
#'  there is no long tails
#'@param ... additional params you can pass to \code{\link[Hmisc]{rcorr}}
#'@inheritParams Hmisc::rcorr
#'@return provides correlation matrix with correlation cofficient \code{r},
#'  p.values  \code{P} and number of overlapping genes \code{n} and oriented
#'  \code{S} matrix
#'@export
#'
#'@seealso \code{\link[Hmisc]{rcorr}} \code{\link{run_fastica}}
#'
#'@examples
#'res_run_ica <- run_fastica (
#'  Example_ds,
#'  optimal = FALSE,
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
           metagenes = data.list,
           threshold = -Inf,
           n.genes.intersect = 30,
           orient.long = TRUE,
           orient.max = FALSE,
           ...) {
    # orient components in the direction of the long tail
    if (orient.long &
        orient.max)
      stop("select one or none orienting method")
    if (orient.long)
      S_or <- S <- .orient_funct(S)
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
#' Assign Independent Component to a metagene through reciprocity
#'
#' This function assign an idependent component to a metagene through
#' verification if for an IC the max correlation points to given metagene and if
#' for this metageneif the max correlation points back the the
#'
#' @param r the correlation matrix, \code{r} matrix from out of
#'   \code{\link{correlate_metagenes}}
#' @param immune_name name of the immune component (if present in \code{r}) to
#'   be excluded from this analysis
#'
#' @return this function returns a data frame with number of IC in the first
#' column and assigned metagene name in second column
#'
#'  @export
#'
#' @examples
#'
#'
#'
assign_metagenes <- function(r, immune_name = "M8_IMMUNE") {
  r <- data.frame(r)
  if (!is.null(immune_name)) {
    r <- r[, -c(match(immune_name, colnames(r)))]
    message("immune metagene excluded")
  } else
    message("no immune metagen name provided")
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
    data.frame(metagene = unlist(sapply(l, "[[", 1)),
               IC = unlist(sapply(l, "[[", 2)))
  message("DONE")
  return(df)
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Identify components related to immune signal
#'
#' @param x the correaltion with immune metagene can be retreived from
#'   \code{\link{correlate_metagenes}} output
#' @param l vector of names of assigned components
#' @param threshold lower bound for filtering correlation [0,1]
#'
#' @return
#' it returns data frame of IC names and correlations passing the
#' \code{threshold}
#' @export
#'
#' @examples
#'
#'
identify_immune_ic <- function(x, l, threshold = 0.1) {
  x[which(x > threshold)][!(names(which(x > threshold)) %in% l)]
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#' Assign through highest correlations
#'
#' It assigns max correlations between set of correlated vectors
#'
#' @param corr list of correlation matrices with correlation coefficients
#' and p-values
#'
#' @return dataframe with column names, correlation coefficient, p.value
#'
#' @export
#'
#' @examples
#'
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

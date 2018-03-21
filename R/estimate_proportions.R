#' Generate markers from metagenes
#'
#' It extracts from set of sources (i.e. ICA S matrix) the n top genes(with
#' weights if needed) to use as marker list or markers with weights for
#' deconvolution through \code{\link{get_scores}}
#'
#' @param df output of \code{run_fastica} containing at least \code{S} and \code{names}
#'   elements
#' @param n number of top genes considered from each signature, n = 30 by
#'   default
#' @param thr max gene expression, if removal of outliers is necessary, Inf (no
#'   threshold) by default.
#' @param sel.ic ICs indentified as specific sources (i.e. immune cells), by
#'   default it takes all ICs of S matrix, can be provided as valid column names
#'   or numeric index
#' @param return should it return "gene.list" or "metagenes.ranked"
#' @param orient.long TRUE by default, if S is oriented change to FALSE
#'
#' @return function returns either list of gene markers "gene.list" for each
#'   signal (selected IC component) or list of "metagenes.ranked" which are gene
#'   nemes with ICA weights
#'
#' @export
#'
#' @examples
#'
#' #gene markers
#' #generate_markers(df, 10)
#'
generate_markers <-
  function(df,
           n = 30,
           thr= Inf,
           sel.ic = paste("IC", 1:ncol(df$S), sep = ""),
           return = "gene.list",
           orient.long = TRUE) {
    if (orient.long) {
      S_or <- .orient_funct(df$S)
    } else {
      S_or <- df$S
    }
    colnames(S_or) <- paste("IC", 1:ncol(S_or), sep = "")
    row.names(S_or) <- df$names
    S_or <- S_or[, sel.ic]
    metagenes <-
      apply(S_or, 2, function(col)
        data.frame(GENE = row.names(S_or), col))

    switch(
      return,
      gene.list = lapply(metagenes, function(x) {
        x <- x[order(-x[, 2]), ]
        x <- x[which(x[, 2] < thr), ]
        return(as.character(x[1:n, 1]))
      }),
      metagenes.ranked = lapply(metagenes, function(x) {
        x <- x[order(-x[, 2]), ]
        return(x[which(x[, 2] < thr), ][1:n, ])
      })
    )
  }

#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Generate basis matrix
#'
#' It generates a basis matrix that can be used for regression from list of
#' weighted markers
#'
#' @param df output of \code{run_fastica} containing at least \code{S} and
#'   \code{names} elements
#' @param sel.ic ICs indentified as specific sources (i.e. immune cells), by
#'   default it takes all ICs of S matrix, can be provided as valid column names
#'   or numeric index
#' @param markers list of markers that should be used for basis matrix (i.e.
#'   "gene.list" from \code{generate_markers}), can be also simple vecor or
#'   list of gene names
#' @param orient.long TRUE by default, if you modified S matrix and you don't want it
#' to be oriented select FALSE
#'
#' @return it returns a \code{data.frame} (basis matrix) that can be used for
#'   regression or visualization purposes
#'
#' @export
#'
#' @examples
#'
#' #generate_basis(df, c("IC1", "IC3", "IC5"), markers = my.gene.list )
#'
#' #visualize
#' #basis <- generate_basis(df, c("IC1", "IC3", "IC5"), markers = my.gene.list )
#' #pheatmap::pheatmap(basis)
#'
generate_basis <- function(df, sel.ic, markers, orient.long = TRUE) {
  if (orient.long) {
    S_or <- .orient_funct(df$S)
  } else {
    S_or <- df$S
  }
  colnames(S_or) <- paste("IC", 1:ncol(df$S), sep = "")
  row.names(S_or) <- df$names
  genes <- unique(as.character(unlist(markers)))
  S_sel <- S_or[genes, sel.ic]
  return(S_sel)
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Get aboundance scores
#'
#'
#' It calulates aboundance scores through a mean of marker genes
#'
#' @param df gene matrix with samples in columns and genes in rows with named
#'   rows
#' @param markers.list list of genes or list of genes with weights
#' @param summary can be any type of mean i.e. \code{mean}, \code{gm_mean}
#'   (geometric mean), \code{harmonic_mean}, \code{weighted.mean}. For wieghted
#'   mean weights are needed along with gene names
#' @param ... optional parameters for the mean function
#'
#' @return Function returns numerical value for each column (sample) of provided data frame
#'
#' @export
#'
#' @examples
#'
#' #get_scores (gene_expression, my_markers.list, summary = "mean", na.rm = TRUE)

get_scores <-
  function(df,
           markers.list,
           summary = "mean",
           ...) {
    ifelse(summary == "weighted.mean",
           type <- "metagenes",
           type <- "gene.list")
    switch(
      type,
      gene.list = sapply(markers.list, function(markers) {
        if (ncol(markers) > 1L)
          markers <- markers[, 1]
        fun <- match.fun(summary)
        common.markers <- intersect(markers, row.names(df))
        if (length(common.markers) < 0.5 * length(markers.list$markers)) {
          warning(paste("not enough markers for"), markers, sep = " ")
          return(NA)
        } else {
          df.m <- df[common.markers, ]
          apply(df.m, 2, fun, ...)
        }
      }),
      metagenes = sapply(markers.list, function(metagene) {
        markers <- metagene[, 1]
        common.markers <- intersect(markers, row.names(df))
        if (length(common.markers) < 0.5 * length(markers.list$markers)) {
          warning(paste("not enough markers for"), markers, sep = " ")
          return(NA)
        } else {
          df.m <- df[common.markers, ]
          apply(df.m, 2, stats::weighted.mean, w = metagene[, 2])
        }
      })
    )
  }

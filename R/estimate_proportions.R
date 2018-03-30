#' Generate markers from components
#'
#' It extracts from set of components (i.e. ICA S matrix) the n top genes (with
#' weights if needed) to use as marker list or markers with weights for
#' estimation of abundance through \code{\link{get_scores}}
#'
#' @param df list (usually output of \code{run_fastica}) containing at least \code{S}
#' and \code{names} elements
#' @param n number of top genes considered from each signature, n = 30 by
#'   default
#' @param thr max gene expression, if removal of outliers is necessary, Inf (no
#'   threshold) by default.
#' @param sel.comp components of interest (i.e. identified as specific to some profiles/metagenes
#' (i.e. immune cells)), by default it takes all columns of S matrix, can be provided as valid
#' column names or numeric index
#' @param return return \code{gene.list} or \code{gene.ranked}
#' @param orient.long \code{TRUE} by default, if \code{S} is oriented change to \code{FALSE}
#'
#' @return function returns either list of gene markers \code{gene.list} for each
#'   component or list of \code{gene.ranked} which are gene names with weights
#'
#' @export
#'
#' @seealso \code{run_fastica}, \code{\link{get_scores}}
#'
#' @examples
#'set.seed(123)
#'  res_run_ica <- run_fastica (
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
#'immune <- identify_immune_comp(corr$r[,"M8_IMMUNE"], assign[, "component"], threshold = 0.1)
#'
#'generate_markers(df = res_run_ica,n = 10,sel.comp= names(immune))
#'generate_markers(df = res_run_ica,n = 10,sel.comp= names(immune), return= "gene.ranked")
generate_markers <-
  function(df,
           n = 30,
           thr = Inf,
           sel.comp = paste("IC", 1:ncol(df$S), sep = ""),
           return = "gene.list",
           orient.long = TRUE) {
    if (orient.long) {
      S_or <- .orient_funct(df$S)
    } else {
      S_or <- df$S
    }
    colnames(S_or) <- paste("IC", 1:ncol(S_or), sep = "")
    row.names(S_or) <- df$names
    S_or <- S_or[, sel.comp]
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
      gene.ranked = lapply(metagenes, function(x) {
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
#' @param sel.comp components identified as specific sources (i.e. immune cells), by
#'   default it takes all components of S matrix, can be provided as valid column names
#'   or numeric index
#' @param markers list of markers that should be used for basis matrix (i.e.
#'   "gene.list" from \code{generate_markers}), can be also simple vector or
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
#'set.seed(123)
#'  res_run_ica <- run_fastica (
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
#'immune <- identify_immune_comp(corr$r[,"M8_IMMUNE"], assign[, "component"], threshold = 0.1)
#'
#'markers <- generate_markers(df = res_run_ica,n = 10,sel.comp= names(immune), return= "gene.list")
#'basis <- generate_basis(df = res_run_ica,sel.comp= names(immune),markers= markers )
#'pheatmap::pheatmap(basis )
generate_basis <- function(df, sel.comp, markers, orient.long = TRUE) {
  if (orient.long) {
    S_or <- .orient_funct(df$S)
  } else {
    S_or <- df$S
  }
  colnames(S_or) <- paste("IC", 1:ncol(df$S), sep = "")
  row.names(S_or) <- df$names
  genes <- unique(as.character(unlist(markers)))
  S_sel <- S_or[genes, sel.comp]
  return(S_sel)
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Get abundance scores
#'
#'
#' It calculates abundance scores through a mean of marker genes
#'
#' @param df gene matrix with samples in columns and genes in rows with named
#'   rows
#' @param markers.list list of genes or list of genes with weights
#' @param summary can be any type of mean i.e. \code{mean}, \code{gm_mean}
#'   (geometric mean), \code{harmonic_mean}, \code{weighted.mean}. For weighted
#'   mean weights are needed along with gene names
#' @param ... optional parameters for the mean function
#'
#' @return Function returns numerical value for each column (sample) of provided data frame
#'
#' @export
#'
#' @examples
#'set.seed(123)
#'  res_run_ica <- run_fastica (
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
#'immune <- identify_immune_comp(corr$r[,"M8_IMMUNE"], assign[, "component"], threshold = 0.1)
#'counts.abs <- (2^res_run_ica$log.counts)-1
#'row.names(counts.abs) <- res_run_ica$names
#'
#'markers <- generate_markers(df = res_run_ica,n = 10,
#'                             sel.comp= names(immune),
#'                             return= "gene.list")
#'get_scores (counts.abs, markers, summary = "mean", na.rm = TRUE)
#'
#'markers <- generate_markers(df = res_run_ica,n  = 10,
#'                            sel.comp= names(immune),
#'                            return= "genes.ranked")
#'get_scores (counts.abs, markers, summary = "weighted.mean", na.rm = TRUE)
get_scores <-
  function(df,
           markers.list,
           summary = "mean",
           ...) {
    type <- NULL
    if (summary == "weighted.mean") {
      type = "metagenes"
    } else {
      type = "gene.list"
    }
    switch(
      type,
      gene.list = sapply(markers.list, function(markers) {
        if (!is.null(ncol(markers)))
          markers <- markers[, 1]
        fun <- match.fun(summary)
        common.markers <- intersect(markers, row.names(df))
        if (length(common.markers) < 0.5 * length(markers.list$markers)) {
          warning(paste("not enough markers for"), markers, sep = " ")
          return(NA)
        } else {
          df.m <- df[common.markers,]
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
          df.m <- df[common.markers,]
          apply(df.m, 2, stats::weighted.mean, w = metagene[, 2])
        }
      })
    )
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Make list of weighted markers
#'
#' Transforms a data frame with multiple columns into a named list of weighted markers
#' with gene names in the first column and values in the second column.
#'
#' @param df \code{data.frame} to be transformed with gene names in the \code{row.names}
#'
#' @return
#' named \code{list} of \code{data.frame}s with gene names in the first column and values in the second column.
#' @export
#'
#' @examples
#' X <- as.data.frame(matrix(runif(10000), 50, 10))
#' row.names(X) <- paste("A",1:nrow(X), sep="")
#' make_list(X)
#'
#'
make_list <-function(df){
  apply(data.frame(df), 2, function(col)
    data.frame(GENE = row.names(data.frame(df)), col))
}

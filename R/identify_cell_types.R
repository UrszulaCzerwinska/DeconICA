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
get_enrichment <-
  function(S,
           gene.names,
           immune.ics,
           gmt = get(utils::data("ImmgeneHUGO",  envir = environment())),
           quantile = 0.99,
           abs.thr = NULL,
           threshold = 0.05,
           ...) {
    # empty list for outputs
    enrich <- list()
    metagenes <- list()
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
    for (i in 2:(length(immune.ics) + 1)) {
      # apply quantile or direct threshold
      if (!(is.null(quantile)))
        abs.thr <- quantile(table[, i], quantile)
      t <- table[which(table[, i] > abs.thr), c(1, i)]
      names <- as.matrix(t[, 1])
      Example <-
        ACSNMineR::enrichment(names, threshold = threshold, maps = gmt, ...)
      message(paste(immune.ics[i - 1]), "...DONE", sep = "")
      # save only if some enrichment detected
      if (length(Example) > 1) {
        # retrun enrichment results ordered by corrected p.value
        enrich[[immune.ics[i - 1]]] <-
          Example[order(Example[["p.value.corrected"]]), ]
        # return actual metagene with weights
        metagenes[[immune.ics[i - 1]]] <- t
      }
    }
    return(list(enrich = enrich,
                metagenes = metagenes))
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Attribute cell type to a component
#'
#' @param enrich enrichment reults from \code{\link{get_enrichment}}
#' @param n look into \code{n} top results
#'
#' @return it returns list of data frames for each result of enrichment list
#'   from  \code{\link{get_enrichment}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CIT.cell.types.immgen <- identify_cell_types(S = res.pipeline.BRCACIT$ica$S,res.pipeline.BRCACIT$ica$names, immune.ics=res.pipeline.BRCACIT$immune.ics)
#' cell_voting_immgene (CIT.cell.types.immgen$enrich, 10)
#' }

cell_voting_immgene <-
  function(enrich, n = 10) {
    sapply(names(enrich), function(j) {
      t <- table(enrich[[j]][1:n, 1])
      data.frame(cell.type = names(t[order(-t)]),
                 vote = paste(round(t[order(-t)] / n * 100, 2), " %", sep =
                                ""))

    }, simplify = FALSE, USE.NAMES = TRUE)
  }


# #
# table(CIT.cell.types.immgen$enrich$IC25[1:10,1:8])/10*100
#
#
# CIT.cell.types.immgen$enrich$IC41 %>% arrange(p.value.corrected)
# CIT.cell.types.immgen <- get_enrichment(S = res.pipeline.BRCACIT$ica$S,res.pipeline.BRCACIT$ica$names, immune.ics="IC51")
# t <- table(CIT.cell.types.immgen$enrich$IC41[1:10,1])
# t[order(-t)]
# BRCACIT.cell.types.immgen <- get_enrichment(S = res.pipeline.BRCACIT$ica$S,res.pipeline.BRCACIT$ica$names, immune.ics=res.pipeline.BRCACIT$immune.ics, threshold = 1)
# cell_voting_immgene (BRCACIT.cell.types.immgen$enrich, 10)
# BRCAWAN.cell.types.immgen <- get_enrichment(S = res.pipeline.BRCAWAN$ica$S,res.pipeline.BRCAWAN$ica$names, immune.ics=res.pipeline.BRCAWAN$immune.ics)
# cell_voting_immgene (BRCAWAN.cell.types.immgen$enrich, 10)
# METABRIC.cell.types.immgen <- get_enrichment(S = res.pipeline.METABRIC$ica$S,res.pipeline.METABRIC$ica$names, immune.ics="IC26")
# cell_voting_immgene (METABRIC.cell.types.immgen$enrich, 10)
#
# BRCACIT.cell.types.immgen$enrich$IC25
# BRCACIT.cell.types.immgen$enrich$IC5

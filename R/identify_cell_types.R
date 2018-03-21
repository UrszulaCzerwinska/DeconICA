#' @title Enrichment analysis
#'
#' @description Computes an enrichment score (fisher exact test) in provided signatures for selected
#'   metagenes
#'
#' @details \code{gene_enrichment_test} runs enrichment of an independent component
#'   (or any ranked list) in known immune cell types signatures. By default it
#'   is using \code{S} matrix from ICA \code{\link{run_fastica}}
#'   \code{\link{fisher.test}} only on components indentified as correlated with
#'   immune metagene through function \code{identify_immune_ic} and it searches
#'   in Immgen signatures \url{http://Immgen.org}.
#'
#' @param S ranks of components, dim \code{n} corresponding to genes, \code{m}
#'   corresponding to number of components, use oriented matrix
#' @param gene.names character vector of gene names, length needs to be equal to
#'   \code{n}
#' @param immune.ics vector of charatcter names of components to use for
#'   enrichment test
#' @param gmt data.frame obtained from gmt file with a function
#'   \code{\link[ACSNMineR:format_from_gmt]{format_from_gmt}}, by default Immgen signatures
#'   \url{http://Immgen.org}
#'@inheritParams ACSNMineR::enrichment
#'@inheritParams stats::fisher.test
#' @param p.adjust.method correction method
#' @param n number of top genes that wlll be used to test signature
#' @param n.consider number of genes from the positive end to be considered
#' @param p.value.threshold maximal p-value (corrected if correction is enabled) that will be displayed
#' @param max_module_size maximum moudule size from gmt file to be considered in enrichment
#' @param orient.long TRUE by default, in case you applied transformation to your S metagenes, select FALSE
#'
#' @return Function returns value if there is an enrichment in provided
#'   signatures: \describe{ \item{metagenes}{interpreted metagene gene ranking}
#'   \item{enrichment}{full results of the enrichment analysis sorted by corrected
#'   p.value} \item{genes.list}{list of genes used for enrichment} }
#'
#' @seealso \code{\link{identify_immune_ic}} identifying immune related
#'   components, \code{\link{run_fastica}} for running Independent Components
#'   Analysis, and \code{\link[ACSNMineR]{enrichment}} for enrichment in gmt
#'   files
#' @export
#'
#' @examples
#'
#'
#'
#'
#'
gene_enrichment_test <-
  function(S,
           gene.names,
           immune.ics,
           gmt = ImmgenHUGO,
           alternative = c("greater", "lower"),
           p.adjust.method = c("holm",
                               "hochberg",
                               "hommel",
                               "bonferroni",
                               "BH",
                               "BY",
                               "fdr",
                               "none"),
           n = 100,
           n.consider = 500,
           min_module_size = 5,
           max_module_size = 500,
           p.value.threshold = 0.05,
           orient.long = TRUE) {
    #extract genes ignoring empty strings
    non.empty.gmt <- sapply(1:nrow(gmt), function(row) {
      gmt[row, which(gmt[row, 3:ncol(gmt)] != "")]
    })
    if (sum(sapply(non.empty.gmt, length) < min_module_size) > 0) {
      message(paste(
        sum(sapply(non.empty.gmt, length) < min_module_size),
        " modules lower than threshold: ",
        min_module_size,
        sep = ""
      ))
      non.empty.gmt <-
        non.empty.gmt[-which(sapply(non.empty.gmt, length) < min_module_size)]
    }
    if (sum(sapply(non.empty.gmt, length) > max_module_size) > 0) {
      message(paste(
        sum(sapply(non.empty.gmt, length) > max_module_size),
        " modules higher than threshold: ",
        max_module_size,
        sep = ""
      ))
      non.empty.gmt <-
        non.empty.gmt[-which(sapply(non.empty.gmt, length) > max_module_size)]
    }
    #count unique genes in gmt
    genes.universe <- unique(unlist(non.empty.gmt))
    #count overlap between provided genes and the gmt
    present <- gene.names %in% genes.universe
    if (sum(present) > 0.3 * length(genes.universe)) {
      message(
        paste(
          "Number of genes present in the gmt from provided list: ",
          sum(present),
          "/",
          length(gene.names),
          sep = ""
        )
      )
    } else {
      stop(
        paste(
          "Not enough overlap between provided list and gmt signatures: ",
          sum(present),
          "/",
          length(gene.names),
          sep = ""
        )
      )
    }
    if (length(as.matrix(gene.names)) != nrow(S))
      stop("wrong number of gene names")
    #ortient if needed
    if (orient.long)
      S <- .orient_funct(S)
    # add colnames
    gene.names <- as.data.frame(gene.names)
    colnames(S) <- paste("IC", 1:ncol(S), sep = "")
    colnames(gene.names) <- "gene.names"
    table <- data.frame(gene.names, S)[c("gene.names", immune.ics)]
    table.uni <- table[present, ]

    message("saving metagenes")
    #the list of IC and genes (no filter)
    metagenes.list <- sapply(immune.ics, function(i) {
      # apply quantile or direct threshold

      table <- table[order(-table[, i]), ]
      t <- table[, c("gene.names", i)]

    }, simplify = FALSE, USE.NAMES = TRUE)

    message("extracting top genes")
    message("")
    # the list of genes for enrichment respecing threshold n
    genes.list <-
      sapply(immune.ics, function(i) {
        t <- table.uni [order(-table.uni [, i]), ]
        t <- t[, c("gene.names", i)]
        names <- as.array(t[, 1])
        }, simplify = FALSE, USE.NAMES = TRUE)

    #number of ros of result

    enrichment <- sapply(immune.ics, function(n.ic) {
      module <- nrow(gmt)
      df <-
        as.data.frame(matrix("NA", module * 9, nrow = module, ncol = 9) ,
                      stringsAsFactors = FALSE)
      colnames(df) <-
        c(
          "module",
          "module_size",
          "nb_genes_in_module",
          "genes_in_module",
          "universe_size",
          "nb_genes_in_universe",
          "p.value",
          "p.value.corrected",
          "test"
        )

      df[["universe_size"]] <- rep(length(genes.universe), module)
      df[["test"]] <- rep(alternative, module)
      df[["nb_genes_in_universe"]] <-
        rep(sum(as.character(gene.names[, 1]) %in% genes.universe), module)
      message(paste("running enrichment for: ", n.ic, sep = ""))
      # for each line of gmt
      for (n.module in seq(1, module, 1)) {
        df[n.module, "module"] <- as.character(gmt[n.module, "V1"])
        df[n.module, "module_size"] <-
          as.character(gmt[n.module, "V2"])
        sig <-
          gmt[n.module, which(gmt[n.module, 3:ncol(gmt)] != "")]
        ica <- as.character(genes.list[[n.ic]])[1:n.consider]

        a <- length(which(which(ica %in% sig) <= n))
        c <- n - a
        b <- length(which(which((ica %in% sig)) > n))
        d <- length(ica) - n - b
        df[n.module, "nb_genes_in_module"] <- as.character(a)
        df[n.module, "genes_in_module"] <-
          paste(ica[which(which(ica %in% sig) <= n)],  collapse = " ")
        M <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        f.p.v <-
          stats::fisher.test(M, alternative = alternative)$p.value
        #ica.matlab.BRCAWAN.enrichment$enrich[[1]]
        df[n.module, "p.value"] <- f.p.v
      }

      message("correcting p.values")
      #multiple testing correction
      df[["p.value.corrected"]] <-
        stats::p.adjust(df[["p.value"]], method = p.adjust.method)
      df <- df[order(df[["p.value.corrected"]]), ]

      message("applying p.value.threshold")
      df <-
        df[which(df[["p.value.corrected"]] < p.value.threshold), ]

    }, simplify = FALSE, USE.NAMES = TRUE)

    enrichment <- enrichment[sapply(enrichment, function(s)
      nrow(s) > 0)]
    metagenes.list <-
      metagenes.list[sapply(enrichment, function(s)
        ! (is.null(s)))]
    genes.list <-
      genes.list[sapply(enrichment, function(s)
        ! (is.null(s)))]
    message("")
    message("DONE")
    return(list(
      metagenes = metagenes.list,
      genes.list = genes.list,
      enrichment = enrichment
    ))
  }

#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Attribute cell type to a component
#'
#' Works only with immgen signatures
#'
#' @param enrich enrichment reults from \code{\link{gene_enrichment_test}}
#' @param n look into \code{n} top results
#'
#' @return it returns list of data frames for each result of enrichment list
#'   from  \code{\link{gene_enrichment_test}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CIT.cell.types.immgen <- identify_cell_types(S = res.pipeline.BRCACIT$ica$S,res.pipeline.BRCACIT$ica$names, immune.ics=res.pipeline.BRCACIT$immune.ics)
#' cell_voting_immgene (CIT.cell.types.immgen$enrich, 10)
#' }


cell_voting_immgen <-
  function(enrich, n = 10) {
    sapply(names(enrich), function(j) {
      if (n > nrow(enrich[[j]]))
        n <- nrow(enrich[[j]])
      t <- table(enrich[[j]][1:n, 1])
      weight <- sapply(names(t), function(name) {
        switch (
          name,
          "alpha.beta.T.cells" =  68,
          "B.cells" = 39,
          "CD.positive" = 8,
          "Fetal.Liver" = 2,
          "gamma.delta.T.cells" = 18,
          "Myeloid.Cells" = 70,
          "NK.cells" = 11,
          "Plasmacytoid" = 3,
          "Stem.Cells" = 11,
          "Stromal.Cells" = 11
        )
      })

      adjusted.val <- t / weight

      data.frame(cell.type = names(adjusted.val[order(-adjusted.val)]),
                 vote = paste(round(
                   adjusted.val[order(-adjusted.val)] / sum(adjusted.val) * 100, 2
                 ), " %", sep = ""))
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

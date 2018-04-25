# go to ./data-raw repo
setwd("inst/data-raw/")
# list all .gmt files
all.raw.gmt <- list.files(".", pattern = "*.gmt$")

# Use funtion form ACSNMineR package to import gmt file
list2env(gmt.list <-
           lapply(
             setNames(all.raw.gmt , make.names(gsub(
               "*.gmt$", "", all.raw.gmt
             ))),
             ACSNMineR::format_from_gmt
           ),
         envir = .GlobalEnv)

#
# --------------------------------------------------
# simplify names and export table of frequencies
#
ImmgenHUGO[, 1] <-
  sapply(strsplit(ImmgenHUGO[, 1], "..",  fixed = TRUE), "[[", 1)

#
# --------------------------------------------------
# export frequencies table
#
# t <-
#   data.frame(table(sapply(
#     strsplit(ImmgeneHUGO[, 1], "..",  fixed = TRUE), "[[", 1
#   )))
# colnames(t) <- c("cell.type", "frequency")
# tab <- ggpubr::ggtexttable(t , rows = NULL, theme = ggpubr::ttheme("classic"))
# ggsave( "./ImmGeneFreqTable.png", tab)
#
# --------------------------------------------------
# export all objects to the ./data
#

devtools::use_data(ImmgenHUGO,
                   TIMER_cellTypes,
                   overwrite = TRUE)




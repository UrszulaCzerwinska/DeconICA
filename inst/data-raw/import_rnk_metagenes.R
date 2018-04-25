# go to ./data-raw repo
setwd("./inst/data-raw/")
# list all .rnk files
all.raw.rnk <- list.files(".", pattern = "*.rnk$")
# user can import specific file with the code
# system.file("data-raw", file, package = "deconica")
# --------------------------------------------------
#
# import rnk files to the environment
list2env(
  rnk.list <-
    lapply(
      setNames(all.raw.rnk, make.names(gsub(
        "*.rnk$", "", all.raw.rnk
      ))),
      read.table,
      sep = "\t",
      strip.white = TRUE,
      stringsAsFactors = FALSE
    ),
  envir = .GlobalEnv
)

Biton.list <- rnk.list
#
# --------------------------------------------------
# export all objects to the ./data
#
devtools::use_data(Biton.list, internal = FALSE, overwrite= TRUE)

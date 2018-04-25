setwd("./data-raw/")
BEK_log_cen <-
  read.csv("BRCABEK.txt",
           sep = "",
           stringsAsFactors = FALSE)
Example_ds <- BEK_log_cen[,c(1, sample(2:ncol(BEK_log_cen),60))]
devtools::use_data(Example_ds, overwrite = TRUE)



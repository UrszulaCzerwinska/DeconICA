LM22 <- read.table("./inst/data-raw/LM22.txt", sep= "\t", header= TRUE, row.names = 1)
LM22.list <- make_list(LM22)
devtools::use_data(LM22.list)

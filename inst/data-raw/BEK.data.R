#setwd("/Users/ulala/Documents/CURIE/R_code/DeconICA/inst/data-raw/")
source("./change_to_hugo_official.R")
library(Biobase)
library(GEOquery)
#library(limma)
getwd()
# load series and platform data from GEO
gset <- getGEO("GSE23720", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

row.names(exprs(gset))


biocLite("hgu133plus2.db")
library("AnnotationDbi")
library("hgu133plus2.db")    ##for Human
names <- select(hgu133plus2.db, row.names(exprs(gset)), c("SYMBOL")) ##  This is just a trying example
names$SYMBOL <- change_to_hugo_official(names$SYMBOL)
Means <- rowMeans(exprs(gset))
names_uniques <- names[!duplicated(names$PROBEID),]

t1 = data.frame(head(Means, 100))
t2 = head(names_uniques,100)

t1= data.frame(Means)
t2 = names_uniques
t1$PROBEID <- row.names(t1)

probe_symbol_mean <- merge(t1,t2, by="PROBEID")


probe_symbol_mean_sorted <- probe_symbol_mean[order(probe_symbol_mean$SYMBOL, -abs(probe_symbol_mean$Means) ), ] #sort by id and reverse of abs(value)
probe_symbol_mean_sorted_u <- probe_symbol_mean_sorted[ !duplicated(probe_symbol_mean_sorted$SYMBOL), ]

ex <- data.frame(PROBEID= row.names(exprs(gset)),exprs(gset))
dim(ex)
df_hugo <- merge(na.omit(probe_symbol_mean_sorted_u), ex, by="PROBEID")
dim(df_hugo)
dim(probe_symbol_mean_sorted_u)
which(is.na(df_hugo$SYMBOL))
row.names(df_hugo) <- df_hugo$SYMBOL

df_hugo$PROBEID<-NULL
df_hugo$Means<-  NULL
df_hugo$SYMBOL<-  NULL


which(is.na(df_hugo))

max(df_hugo)

df_hugo_unlog <- (2^df_hugo)-1
row.names(df_hugo_unlog)


BEK <-df_hugo_unlog

###### run ICA in MATLAB

BEK_ica_overdecompose <- run_fastica (
  BEK,
  isLog = FALSE,
  overdecompose = TRUE,
  with.names = FALSE,
  gene.names = row.names(BEK),
  R = FALSE
)
devtools::use_data(BEK_ica_overdecompose, overwrite = TRUE)

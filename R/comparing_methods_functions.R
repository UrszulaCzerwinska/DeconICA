########################### WORKING
best.correlations <- function(x,y, mean = TRUE, ...){
  cor.fill <-corr.scores.plot(x,y, ...)$filtered
  assign <- assign_metagenes(cor.fill, immune_name = NULL)
  res <- data.frame(assign, value = apply(assign,1, function(row) cor.fill [row[2], row[1]]))
  if(mean) {
    return(mean(res[,"value"]) )
  } else {
    return(res)
  }
}

make_list <-function(df){
  apply(data.frame(df), 2, function(col)
    data.frame(GENE = row.names(data.frame(df)), col))
}

increase_sparseness <- function(marker.list, df, n) {
  lapply(marker.list, function(sig)
  {
    sp <- apply(df[sig, ], 1, function(row)
      NMF::sparseness(row))
    names(sp[order(-sp)][1:n])
  })
}




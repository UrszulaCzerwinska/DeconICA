#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

change_to_hugo_official <- function(vector.of.names) {
  rownames.test <- vector.of.names
  rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
  rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
  rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
  rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
  require(org.Hs.eg.db)
  require(RSQLite)

  dbCon <- org.Hs.eg_dbconn()
  # write your SQL query
  
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  # execute the query on the database
  aliasSymbol <- dbGetQuery(dbCon, sqlQuery)

  hg.rows.na<-which(!( rownames.test %in% aliasSymbol[,5]  ))
  correspondance <- lapply(rownames.test[hg.rows.na], function (x) if( x %in% aliasSymbol[,2] == TRUE) {  x = aliasSymbol[which(aliasSymbol[,2] == x),5][1] } else {x})
  hugo <- unlist(correspondance)
  
  
  HUGO <-   rownames.test
  HUGO[hg.rows.na] <- hugo
  
  return(HUGO)
  
}

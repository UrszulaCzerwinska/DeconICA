MetICA_consistent<-function(source_list,limit){
  
  # This function evaluate the stability of MetICA component
  # source_list is a list object, in which each element is the centrotype or estimated source (sample in rows)
  # limit is the minimal spearman's correlation to accept the two sources are correlated, between 0 and 1
  
  nb_source=ncol(source_list[[1]]) # Nb of sources
  score_source=rep(1,nb_source)
  for (i in 1:(length(source_list)-1)){
    cor_matrix=cor(source_list[[i]],source_list[[i+1]],method='spearman')
    max_cor=apply(abs(cor_matrix),1,max)
    uncor_source=which(max_cor<limit)
    for (j in 1:nb_source){
      if (score_source[j]==1 & j %in%uncor_source){score_source[j]=0}}
  }
  # It returns a vector having the same length as evaluated components
  # 1 present the source is present in every simulation, otherwise 0
  
  return(score_source)}
MetICA_cluster_generator<-function(S,type_correlation=c('pearson','spearman'),max_cluster){
  
  # This function finds out clusters from estimated sources and outputs 3 files necessary for CCA(Matlab) visualization
  # S: estimated source matrice from MetICA_source_generator, n estimations = n coloums
  # type_correlation: type of correlation, 'pearson' or 'spearman' 
  # max_cluster: maximal number of partitions tested
  
  type_correlation <- match.arg(type_correlation)
  if (max_cluster<2)
  {message("'max_cluster' is too small: reset to ", 2)
   max_cluster=2}
  
  colnames(S)=paste0('IC',1:ncol(S))
  write.table(S,file='source_list.txt',sep='\t',row.names=F) 
  print('First file for CCA exported: source_list.txt') 
  
  R=cor(S,method=type_correlation) # Correlation matrix 
  dist_R=as.dist(1-abs(R))  # Disimilariy Matrix in dist format
  dist_R_matrix=data.matrix(dist_R) # Disimilariy Matrix in matrix format
  system.time(write.table(dist_R_matrix,file='distance.txt',sep='\t',row.names=F,col.names=F))
  print('Second file for CCA exported: distance.txt') 
  
  clusterObj <- hclust(dist_R, method="average") # Hierarchical clustering
  cluster_summary=c() 
  for (nb_cluster in 2:max_cluster){  
    cluster <- cutree(clusterObj,nb_cluster) 
    cluster_summary=cbind(cluster_summary,cluster)} 
  write.table(cluster_summary,file='cluster_labels.txt',col.names=F,row.names=F,sep='\t')
  print('Third file for CCA exported: cluster_labels.txt') 
  # Partition results: each column=sample label for a number of partitions given
  
  print('Cluster generation finished\n')
  
  return(list(S=S,D=dist_R_matrix,C=cluster_summary))
}

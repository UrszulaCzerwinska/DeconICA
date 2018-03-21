MetICA_cluster_center<-function(S,D,Nb){
  
  # This function finds out the centrotype of each cluster
  # S: estimated sources, generated from MetICA_cluster_generator
  # D: dissimilarity matrix, generated from MetICA_cluster_generator
  # Nb: desired number of clusters, decided by the evaluation of CCA in matlab
  
  if (ncol(S)!=ncol(D) || nrow(D)!=ncol(D)){stop('Matrix dimension not correct!')}  
  if (Nb<2)
  {message("'Nb' is too small: reset to ", 2)
   Nb=2}
  
  clusterObj <- hclust(as.dist(D),method="average") # Hierarchical clustering
  cluster <- cutree(clusterObj,Nb)  # Label each estimate into one cluster
  cluster_center<- c() # Matrix contains cluster centers (sources)
  center_ID<-c()  # The ID of the corresponding estimate, from this ID we could know which fastICA run produces this centrotype
  for (p in 1:Nb){
    cl=which(cluster==p) # Indices of estimates that belong to cluster p
    Si=S[,cl] # Estimated sources belong to cluster p
    dist_R_matrix_cluster=D[cl,cl] # Distance matrix for estimates belonging to this cluster
    mini_dis=which.min(apply(dist_R_matrix_cluster,1,sum))  # Which estimates has minimal distance to other points
    cluster_center=cbind(cluster_center,Si[,mini_dis])
    center_ID=c(center_ID,as.double(strsplit(names(mini_dis),'IC')[[1]][2]))}
  
  print('Cluster center calculation finished\n')
  
  return(list(center=cluster_center,center_ID=center_ID))
}

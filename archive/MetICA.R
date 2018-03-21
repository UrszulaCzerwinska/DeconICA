## Loading all necessary packages and files before running the program.

packages <- c("fastICA", "MASS", "e1071")

for (p in 1:3){
  x=packages[p]
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
  } 
}

library(fastICA)
library(MASS)
library(e1071)

source('MetICA_fastICA.R')
source('MetICA_simulated_generator.R')
source('MetICA_source_generator.R')
source('MetICA_cluster_generator.R')
source('MetICA_cluster_center.R')
source('MetICA_bootstrap.R')
source('MetICA_consistent.R')

t=as.numeric(Sys.time())
set.seed((t - floor(t))*1e8->seed) 


##############
new_data=data.matrix(read.table('Yeast-Experimental.txt',sep='\t',dec='.',header=T,check.names=FALSE))
new_data=new_data[2:nrow(new_data),2:ncol(new_data)]
row.names(new_data)=read.table('Yeast-Experimental.txt',sep='\t',dec='.',header=T)[2:(nrow(new_data)+1),1] 
new_data_centered=scale(new_data,scale=F)
dim(new_data_centered)
###
M1=MetICA_source_generator(new_data_centered,0.9999,'gaussian',100)
dim(M1)

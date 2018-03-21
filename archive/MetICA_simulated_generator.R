MetICA_simulated_generator<-function(X,I,PC_scores,PC_loadings,max_iter){
  
  # This function produce simulated data from experimental data X
  # X: original data matrix, must be n (nb of observations) * p (number of features), either centered or not
  # I: level of background noises
  # PC_scores, PC_loadings: PC scores and loadings of X used for generation
  # max_iter: repetitions for background noise generation, the final noise used is the average of repetitions
  
  p=ncol(X)
  C=cov(X)
  background=matrix(0,45,p)
  for (i in 1:max_iter){
    print(paste0('Iteration:',i))
    background=background+mvrnorm(45,rep(0,p),C)}
  X_output=PC_scores%*%t(PC_loadings)+background*I/max_iter
  
  print('Production of simulated data finished\n')
  
  return(X_output)
}
MetICA_source_generator<-function(X,n.comp,w.distribution,max_iter, type){
 X = as.matrix(X)
  # This function computes estimated sources from randomly-initialized fastICA algorithm
  # X: data matrix, must be n (nb of observations) * p (number of features), either centered or not
  # tot_var: minimal pourcentage of variance kept for ICA
  # w.distribution: type of distribution used for generation of random initial demixing matrix W, 'gaussian', 'uniform' or 'beta'
  # max_iter:number of simulations, set the highest possible based on computer memory

  ### Input max_iter, at least 50

  ### Iterations
  # Type of fastICA for the initial half of simulations
  W_sum=c()  # Matrix storing all simulated mixing matrices
  W0_sum=c() # Matrix storing the initial demixing matrices
  A_sum=c() # Matrix storing the loading matrices

  for (i in 1:max_iter){
    print(paste0('Iteration:',i))
    wines.ica <- fastICA2(X, n.comp, w.distribution=w.distribution, alg.typ = type,
                          fun = "logcosh", alpha = 1, maxit = 300, tol = 1e-04)
    W_dmix=t(wines.ica$K%*%wines.ica$W) # Calculation of demixing matrix W
    W_sum=rbind(W_sum,W_dmix) # Storage of demixing matrix
    W0_sum=rbind(W0_sum,wines.ica$W0) # Storage of initial demixing matrix
    A_sum=rbind(A_sum,wines.ica$A) # Storage of loading matrix
  }

  source_list=X%*%t(W_sum) # Matrix storing estimated sources from all runs

  ### Algorithm output: combined source matrix, demixing matrix, initial inputs, number of ICs

  print('Source generation finished,number of components:\n')

  return(list(S=source_list,W=W_sum,W0=W0_sum,A=A_sum,IC=wines.ica$IC))}

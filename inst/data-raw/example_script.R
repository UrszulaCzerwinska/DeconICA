install.packages("devtools")
library(devtools)
instal_github("UrszulaCzerwinska/DeconICA")




 S <- matrix(runif(10000), 10, 2)
 A <- matrix(sample(-3:3, 16, replace = TRUE),2,8, byrow = TRUE)
 X <- data.frame(S %*% A)
 run_fastica(X, row.center = TRUE, n.comp = 3, overdecompose = FALSE, R = FALSE, name = "test")

 print("DONE")
#### if errors try
 # matlabpath = " "  #provide your path ie "/bioinfo/opt/build/Matlab-R2013a/bin" - no slash at the end no "matlab" mention

 # S <- matrix(runif(10000), 10, 2)
 # A <- matrix(sample(-3:3, 16, replace = TRUE),2,8, byrow = TRUE)
 # X <- data.frame(S %*% A)
 # run_fastica(X, row.center = TRUE, n.comp = 3, overdecompose = FALSE, R = FALSE, name = "test", matlbpth = matlabpath)
 #

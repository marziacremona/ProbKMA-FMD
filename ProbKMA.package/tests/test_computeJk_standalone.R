set.seed(123)  # Set the desired seed value in R for obtaining the same results
library(microbenchmark)
source("../ProbKMA-FMD_functions.r")

#alpha = 0.5
#use0 = TRUE
#use1 = TRUE


#alpha = 0.0
#use0 = TRUE
#use1 = FALSE


alpha = 1.0
use0 = FALSE
use1 = TRUE

w = 1
m = 2

my_output_d1 <- compute_Jk_rcpp(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1)
my_output_efficient_d1 <- compute_Jk_rcpp_efficient(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1)
prof_output_d1 <- .compute_Jk(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1)

risultati <- microbenchmark(compute_Jk_rcpp(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1),compute_Jk_rcpp_efficient(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1))
print(risultati)

risultati2 <- microbenchmark(compute_Jk_rcpp(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1),.compute_Jk(V_new[[2]],S_k[[2]],P_k[[2]],Y,alpha,w,m,use0,use1))
print(risultati2)

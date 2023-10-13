# new test for find_min_diss_rcpp
library(microbenchmark)
# caricare data.Rdata

y = list(output_prof$Y0$c1, output_prof$Y1$c1)
v = output_prof$V0

set.seed(123)  # Set the desired seed value in R for obtaining the same results

alpha = 0.5

use0 = TRUE
use1 = TRUE

w = 1 
#c_k = 2
d = 1

# test with distance H1
output1 = .find_min_diss(y,v,alpha,w,c_k,d,use0,use1)
output1
output2 =  find_diss_normal(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2)
output2
output3 = find_diss(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2)
output3

# test with distance L2

alpha = 0.0
use0 = TRUE
use1 = FALSE

output1 = .find_min_diss(y,v,alpha,w,c_k,d,use0,use1)
output1
output2 =  find_diss_normal(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2)
output2
output3 = find_diss(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2)
output3



# test with distance L2 in derivative

alpha = 1.0
use0 = FALSE
use1 = TRUE

output1 = .find_min_diss(y,v,alpha,w,c_k,d,use0,use1)
output1
output2 =  find_diss_normal(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2)
output2
output3 = find_diss(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2)
output3


# execution time comparison 
risultati <- microbenchmark(.find_min_diss(y,v,alpha,w,c_k,d,use0,use1), find_diss(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2))
print(risultati)
  
tempo_funzione1 <- system.time(.find_min_diss(y,v,alpha,w,c_k,d,use0,use1))
tempo_funzione2 <- system.time( find_diss(y,v,w, alpha, c_k, d, use0, use1, .domain, .select_domain, .diss_d0_d1_L2))
tempo_funzione1 
tempo_funzione2  

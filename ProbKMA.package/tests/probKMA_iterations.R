# Set working directory to the folder of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # mi metto nel path corrente 

source("../ProbKMA-FMD_functions.r")
#source("../ProbKMA-FMD_functions_modified_test_elongation.r")
#source("../ProbKMA-FMD_functions_modified_block_elongate.r")
set.seed(123)  # Set the desired seed value in R for obtaining the same results
load(paste0('len200_sd0.1.RData')) 

# probkma iterations diss='d0_d1_L2' 
diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0 
trials_elong = 201 
c_max = 71 
K = 2
c = 61
my_output = probKMA(Y0,Y1,standardize=FALSE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                 diss=diss,alpha=alpha,w=1,m=2,
                 iter_max=10,stop_criterion='max',quantile=NULL,tol=1e-8,
                 iter4elong=2,tol4elong=1e-3,max_elong=0.5,trials_elong=10,deltaJk_elong=0.05,max_gap=max_gap,
                 iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                 return_options=TRUE,return_init=TRUE,worker_number=NULL)

output_elongate_block_rcpp_d0d1 <- my_output
save(output_elongate_block_rcpp_d0d1,file = "output_elongate_block_rcpp_d0d1.RData")

setequal(my_output,output_elongate_block_rcpp_d0d1)

# probkma iterations diss = 'd0_L2' 
diss = 'd0_L2' 
alpha = 0.0
max_gap = 0 
iter4elong = 1 
trials_elong = 201 
c_max = 71 
K = 2
c = 61
output_d0 = probKMA(Y0,Y1,standardize=FALSE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                 diss=diss,alpha=alpha,w=1,m=2,
                 iter_max=10,stop_criterion='max',quantile=NULL,tol=1e-8,
                 iter4elong=2,tol4elong=1e-3,max_elong=0.5,trials_elong=10,deltaJk_elong=0.05,max_gap=max_gap,
                 iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                 return_options=TRUE,return_init=TRUE,worker_number=NULL)

output_elongate_block_rcpp_d0 <- output_d0
save(output_elongate_block_rcpp_d0,file = "output_elongation_rcpp_d0.RData")

setequal(output_d0,output_elongate_block_rcpp_d0)

# probkma iterations diss = 'd1_L2' 

diss = 'd1_L2' 
alpha = 1.0
max_gap = 0 
iter4elong = 1 
trials_elong = 201 
c_max = 71 
K = 2
c = 61
output_d1 = probKMA(Y0,Y1,standardize=FALSE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                 diss=diss,alpha=alpha,w=1,m=2,
                 iter_max=10,stop_criterion='max',quantile=NULL,tol=1e-8,
                 iter4elong=2,tol4elong=1e-3,max_elong=0.5,trials_elong=10,deltaJk_elong=0.05,max_gap=max_gap,
                 iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                 return_options=TRUE,return_init=TRUE,worker_number=NULL)

output_elongate_block_rcpp_d1 <- output_d1
save(output_elongate_block_rcpp_d1,file = "c.RData")

setequal(output_d1,output_elongate_block_rcpp_d1)

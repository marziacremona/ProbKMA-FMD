setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
source("../ProbKMA-FMD_functions.r")
set.seed(123)  # Set the desired seed value in R for obtaining the same results
load(paste0('substitution_rates.RData')) 
load(paste0('indel_rates.RData')) 

Y0 = mapply(function(sub, indel){
  cbind(sub, indel)
}, substitution_smoothed, indel_smoothed, SIMPLIFY = FALSE)
Y1 = mapply(function(sub, indel){
  cbind(sub, indel)
}, substitution_smoothed_derivative, indel_smoothed_derivative, SIMPLIFY = FALSE)


# load DataDebugMatDataElongate.RData

result <- elongate_motifs(V_new,
                          V_dom,
                          S_k,
                          P_k,
                          Y,
                          w, 
                          m, 
                          use0,
                          use1,
                          alpha,
                          c,
                          c_max, 
                          max_elong, 
                          deltaJk_elong,
                          trials_elong,
                          D,
                          K,
                          max_gap,
                          .compute_motif,
                          .compute_Jk,
                          .domain,
                          .select_domain,
                          .diss_d0_d1_L2)

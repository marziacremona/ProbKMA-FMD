# Set working directory to the folder of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../../ProbKMA-FMD_functions.r")
source("../../ProbKMA-FMD_functions_modified.r")

#############################
### SIMULATION SCENARIO 1 ###

load(paste0('../Data/indel_rates.RData'))
load(paste0('../Data/substitution_rates.RData'))

Y0 = mapply(function(sub, indel){
  
  cbind(sub, indel)
  
}, substitution_smoothed, indel_smoothed, SIMPLIFY = FALSE)

Y1 = mapply(function(sub, indel){
  
  cbind(sub, indel)
  
}, substitution_smoothed_derivative, indel_smoothed_derivative, SIMPLIFY = FALSE)

#######################
###     RUN FMD     ###
#######################

# use Sobolev-like distance d_0.5
diss = 'd0_d1_L2' 
alpha = 0.5

max_gap = 0 # no gaps allowed
iter4elong = 1 # perform elongation
trials_elong = 201 # try all possible elongations
c_max = 61 # maximum motif length 61

### run probKMA multiple times (2x2x3=12 times)
K = c(2) # number of clusters to try
c = c(51, 41) # minimum motif lengths to try
n_init = 3 # number of random initializations to try
 
files = list.files('./test_motif_search')
if('test.RData' %in% files){
  # candidate motifs already present, load them
  load('./test_motif_search/test.RData')
}else{
# find candidate motifs
  find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,
                                                        name = './test_motif_search/test', names_var = 'x(t)',
                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                               return_options = TRUE, return_init = TRUE,
                                                                               diss = diss, alpha = alpha),
                                                        plot = FALSE, worker_number = NULL)
  save(find_candidate_motifs_results, file = './test_motif_search/test.RData')
}
### filter candidate motifs based on silhouette average and size
silhouette_average = Reduce(rbind, Reduce(rbind, find_candidate_motifs_results$silhouette_average_sd))[ , 1] # retrieve silhouette average for all candidate motifs
filter_candidate_motifs_results = filter_candidate_motifs(find_candidate_motifs_results,
                                                          sil_threshold = quantile(silhouette_average, 0.9),
                                                          size_threshold = 5)

### cluster candidate motifs based on their distance and select radii
cluster_candidate_motifs_results = cluster_candidate_motifs(filter_candidate_motifs_results,
                                                            motif_overlap = 0.6)

save(cluster_candidate_motifs_results,file ='./test_motif_search/cluster_candidate_motifs_results.RData')





# Set working directory to the folder of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../ProbKMA-FMD_functions.r")




#############################
### SIMULATION SCENARIO 1 ###
#############################

# 20 curves (N = 20) with 2 different motifs (nmotifs = 2) of length 60 (len_motifs = 60, corresponding to 61 evaluation points)
# curves 1-6 contain one occurrence of motif 1
# curves 7-12 contain one occurrence of motif 2
# curves 13-14 contain one occurrence of motif 1 and one occurrence of motif 2
# curves 15-16 contain two occurrences of motif 1
# curves 17-18 contain two occurrences of motif 2
# curves 19-20 contain no motif

#################################################################
#                                                               #
# curves of length 200 (corresponding to 201 evaluation points) #
#                                                               #
# level of noise (sigma) 0.1                                    #
#                                                               #
#################################################################

load(paste0('len200_sd0.1.RData')) 

# Y0: list of 20 vectors (univariate curves) with the evaluation of the 20 curves in 201 points
# Y1: list of 20 vectors (univariate curves) with the evaluation of the 20 curve derivatives in 201 points


#######################
###     RUN FMD     ###
#######################

# use Sobolev-like distance d_0.5
diss = 'd0_d1_L2' 
alpha = 0.5

max_gap = 0 # no gaps allowed
iter4elong = 1 # perform elongation
trials_elong = 201 # try all possible elongations
c_max = 71 # maximum motif length 70



### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(61, 51, 41) # minimum motif lengths to try
n_init = 10 # number of random initializations to try

###
# NOTE: rename "results" folder to re-run everything (TIME CONSUMING, TAKES ABOUT 30 MINUTES ON 3 CORES)!
###

files = list.files('./results')
if('len200_sd0.1_candidate.RData' %in% files){
  # candidate motifs already present, load them
  load('./results/len200_sd0.1_candidate.RData')
}else{
  # Y0: list of N vectors, for univariate curves y_i(x), or
  #     list of N matrices with d columns, for d-dimensional curves y_i(x),
  #     with the evaluation of curves (all curves should be evaluated on a uniform grid).
  #     When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
  # Y1: list of N vectors, for univariate derivative curves y'_i(x), or
  #     list of N matrices with d columns, for d-dimensional derivatibe curves y'_i(x),
  #     with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
  #     When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions.
  #     Must be provided when diss='d1_L2' or diss='d0_d1_L2'.

  # find candidate motifs
  find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,
                                                        name = './results/len200_sd0.1', names_var = 'x(t)',
                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                               return_options = TRUE, return_init = TRUE,
                                                                               diss = diss, alpha = alpha),
                                                        plot = TRUE, worker_number = NULL)
  save(find_candidate_motifs_results, file = './results/len200_sd0.1_candidate.RData')
}


### filter candidate motifs based on silhouette average and size
silhouette_average = Reduce(rbind, Reduce(rbind, find_candidate_motifs_results$silhouette_average_sd))[ , 1] # retrieve silhouette average for all candidate motifs
filter_candidate_motifs_results = filter_candidate_motifs(find_candidate_motifs_results,
                                                          sil_threshold = quantile(silhouette_average, 0.9),
                                                          size_threshold = 5)

### cluster candidate motifs based on their distance and select radii
cluster_candidate_motifs_results = cluster_candidate_motifs(filter_candidate_motifs_results,
                                                            motif_overlap = 0.6)

### plot cluster candidate motifs results
pdf('./results/len200_sd0.1_clustering_candidate_motifs.pdf', height = 12, width = 9)
cluster_candidate_motifs_plot(cluster_candidate_motifs_results, ask = FALSE)
dev.off()

### search selected motifs
motifs_search_results = motifs_search(cluster_candidate_motifs_results,
                                      use_real_occurrences = FALSE, length_diff = +Inf)
        
### plot FMD results (NB: no threshold of frequencies of motif found!)
pdf('./results/len200_sd0.1_results.pdf', height = 7, width = 17)
motifs_search_plot(motifs_search_results, ylab = 'x(t)', freq_threshold = 1)
dev.off()

save(find_candidate_motifs_results, silhouette_average, filter_candidate_motifs_results,
     cluster_candidate_motifs_results, motifs_search_results,
     file='./results/len200_sd0.1.RData')

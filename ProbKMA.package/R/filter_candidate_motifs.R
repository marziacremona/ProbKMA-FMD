
filter_candidate_motifs <- function(find_candidate_motifs_results,sil_threshold=0.5,size_threshold=2,
                                    K=find_candidate_motifs_results$K,c=find_candidate_motifs_results$c){
  # Filter the candidate motifs on the basis of a threshold on the average silhouette index
  # and a threshold on the size of the curves in the motif. 
  # find_candidate_motifs_results: output of find_candidate_motif function.
  # sil_threshold: threshold on the average silhouette index.
  # size_threshold: threshold on the size of the motif (number of curves in the cluster).
  # K: vector with numbers of motifs that must be considered (should be a subset of find_candidate_motifs_results$K).
  # c: vector with minimum motifs lengths that must be considered (should be a subset of find_candidate_motifs_results$c).
  
  ### check input #############################################################################################
  # check sil_threshold
  if(length(sil_threshold)!=1)
    stop('sil_threshold not valid.')
  if((sil_threshold<(-1))|(sil_threshold>1))
    stop('sil_threshold should be a number bewteen -1 and 1.')
  # check size_threshold
  if(length(size_threshold)!=1)
    stop('size_threshold not valid.')
  if(size_threshold%%1!=0)
    stop('size_threshold should be an integer.')
  if(size_threshold<1)
    stop('size_threshold should be at least 1.')
  # check K
  if(TRUE %in% !(K %in% find_candidate_motifs_results$K)){
    warning('K is not a subset of find_candidate_motifs_results$K. Using default K.')
    K=find_candidate_motifs_results$K
  }
  # check c
  if(TRUE %in% !(c %in% find_candidate_motifs_results$c)){
    warning('c is not a subset of find_candidate_motifs_results$c. Using default c.')
    c=find_candidate_motifs_results$c
  }
  
  ### filter motifs ###########################################################################################  
  n_init=find_candidate_motifs_results$n_init
  name=find_candidate_motifs_results$name
  motifs=lapply(K,
                function(K){
                  lapply(c,
                         function(c){
                           lapply(seq_len(n_init),
                                  function(i){
                                    load(paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
                                    motifs=.probKMA_silhouette_filter(probKMA_results,silhouette,sil_threshold,size_threshold)
                                    return(motifs)
                                  })
                         })
                })
  motifs=unlist(unlist(motifs,recursive=FALSE),recursive=FALSE)
  if(is.null(unlist(motifs)))
    stop("No motif present after filtering. Please re-run the function with less stringent parameters.")
  V0_clean=unlist(lapply(motifs,function(motifs) motifs$V0_clean),recursive=FALSE)
  V_clean_length=unlist(lapply(V0_clean,length))
  index=order(V_clean_length,decreasing=TRUE) # order from the longest to the shortest
  V0_clean=V0_clean[index]
  V1_clean=unlist(lapply(motifs,function(motifs) motifs$V1_clean),recursive=FALSE)[index]
  D_clean=Reduce(cbind,lapply(motifs,function(motifs) motifs$D_clean))[,index]
  P_clean=Reduce(cbind,lapply(motifs,function(motifs) motifs$P_clean))[,index]
  c=Reduce(c,lapply(motifs,function(motifs) motifs$c))[index]
  K=Reduce(c,lapply(motifs,function(motifs) motifs$K))[index]
  
  ### output ##################################################################################################
  load(paste0(name,"_K",K[1],"_c",c[1],'/random',1,'.RData'))
  return(list(V0_clean=V0_clean,V1_clean=V1_clean,D_clean=as.matrix(D_clean),P_clean=as.matrix(P_clean),c=c,K=K,
              Y0=probKMA_results$Y0,Y1=probKMA_results$Y1,
              diss=probKMA_results$diss,alpha=probKMA_results$alpha,w=probKMA_results$w,max_gap=probKMA_results$max_gap))
}
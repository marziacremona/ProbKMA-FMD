
.probKMA_silhouette_filter <- function(probKMA_results,silhouette,sil_threshold=0.5,size_threshold=2){
  # Filter the motifs found by probKMA on the basis of a threshold on the average silhouette index
  # probKMA_results: output of probKMA function (with return_options=TRUE)
  # silhouette: output of probKMA_silhouette function
  # sil_threshold: threshold on the average silhouette index
  # size_threshold: threshold on the size of the motif (number of curves in the cluster)
  
  index_sil=which(silhouette$silhouette_average>=sil_threshold)
  index_size=which(colSums(probKMA_results$P_clean)>=size_threshold)
  index=intersect(index_sil,index_size)
  if(length(index)==0)
    return(NULL)
  return(list(V0_clean=probKMA_results$V0_clean[index],V1_clean=probKMA_results$V1_clean[index],
              D=probKMA_results$D[,index],D_clean=probKMA_results$D_clean[,index],P=probKMA_results$P[,index],P_clean=probKMA_results$P_clean[,index],
              c=probKMA_results$c[index],K=rep(probKMA_results$K,length(index))))
}
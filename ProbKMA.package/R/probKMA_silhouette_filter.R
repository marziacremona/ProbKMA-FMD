#' @title .probKMA_silhouette_filter
#'
#' @description Filter the motifs found by probKMA on the basis of a threshold on the average silhouette index
#'
#' @param probKMA_results output of probKMA function (with return_options=TRUE)
#' @param silhouette output of probKMA_silhouette function
#' @param sil_threshold threshold on the average silhouette index
#' @param size_threshold threshold on the size of the motif (number of curves in the cluster)
#' @return A list containing all the elements in probKMA_results correspondent to the filtered motifs;
#' @return \item{V0_clean}{ filtered motifs}
#' @return \item{V1_clean}{ derivative of the filtered motifs}
#' @return \item{D}{ filtered dissimilarity matrix}
#' @return \item{S_clean}{ shift warping matrix after cleaning motifs}
#' @return \item{D}{ dissimilarity matrix}
#' @return \item{D_clean}{ filtered dissimilarity matrix after cleaning motifs}
#' @return \item{P}{ filtered membership matrix}
#' @return \item{P_clean}{ filtered membership matrix after cleaning motifs}
#' @return \item{c}{ filtered minimum motif lengths}
#' @return \item{K}{ vector containing the number of motifs repeated a number of filtered motifs times}
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
.probKMA_silhouette_filter <- function(probKMA_results,silhouette,sil_threshold=0.5,size_threshold=2){
 
  index_sil=which(silhouette$silhouette_average>=sil_threshold)
  index_size=which(colSums(probKMA_results$P_clean)>=size_threshold)
  index=intersect(index_sil,index_size)
  if(length(index)==0)
    return(NULL)
  return(list(V0_clean=probKMA_results$V0_clean[index],V1_clean=probKMA_results$V1_clean[index],
              D=probKMA_results$D[,index],D_clean=probKMA_results$D_clean[,index],P=probKMA_results$P[,index],P_clean=probKMA_results$P_clean[,index],
              c=probKMA_results$c[index],K=rep(probKMA_results$K,length(index))))
}
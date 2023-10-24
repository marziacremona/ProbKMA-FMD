#' @title cluster_candidate_motifs
#'
#' @description Determine a global radius Rall (based on distances between all motifs and all curves), and cluster candidate
#' motifs based on their distance, requiring groups to be more than 2*Rall apart. Then for each group m=1,...,M
#' determine a group-specific radius Rm (based on distances between motifs of the same group and all curves).
#'
#' @param filter_candidate_motifs_results output of filter_candidate_motifs function.
#' @param motif_overlap minimum overlap required between candidate motifs in their distance computation,
#' in % of the shortest motif.
#' @param k_knn number of neighbors to be used in k-nearest neighbors, in order to determine Rall and Rm
#'when the two groups (curves with/without motif) are not separated.
#' @param votes_knn_Rall threshold on the percentage of votes for class 1 (curve has the motif) in
#' k-nearest neighbors, in order to determine Rall.
#' @param votes_knn_Rm threshold on the percentage of votes for class 1 (curve has the motif) in
#' k-nearest neighbors, in order to determine Rm.
#' @param worker_number number of CPU cores to be used for parallelization (default number of CPU cores).
#' If worker_number=1, the function is run sequentially.
#' @return A list containing some its input: k_knn, votes_knn_Rall and votes_knn_Rm;
#' the outputs of filter_candidate_patterns function;
#' @return \item{VV_D}{ matrix of distances between candidates motifs after filtering}
#' @return \item{VV_S}{ matrix of shift deformations of candidates motifs after filtering}
#' @return \item{R_all}{ global radius used to cut the dendrogram}
#' @return \item{hclust_res}{list of result hierarchical clustering based on the distances candidate motifs obtained after filtering}
#' @author Marzia Angela CREMONA & Francesca CHIAROMONTE
#' @export
cluster_candidate_motifs <- function(filter_candidate_motifs_results,motif_overlap=0.6,
                                     k_knn=3,votes_knn_Rall=0.5,votes_knn_Rm=0.5,worker_number=NULL){
 
  ### set parallel jobs ###################################################################################
  core_number <- parallel::detectCores()
  # check worker number
  if(!is.null(worker_number)){
    if(!is.numeric(worker_number)){
      warning('Worker number not valid. Selecting default number.')
      worker_number=NULL
    }else{
      if((worker_number%%1!=0)|(worker_number<1)|(worker_number>core_number)){
        warning('Worker number not valid. Selecting default number.')
        worker_number=NULL
      }
    }
  }
  if(is.null(worker_number))
    worker_number <- core_number-1
  rm(core_number)
  if(worker_number>1){
    cl_search=parallel::makeCluster(worker_number,timeout=60*60*24*30)
    on.exit(parallel::stopCluster(cl_search))
  }else{
    cl_search=NULL
  }
  
  ### prepare input data ##################################################################################
  if(filter_candidate_motifs_results$diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(filter_candidate_motifs_results$Y0,function(y0) list(y0=y0,y1=NULL))
    V=lapply(filter_candidate_motifs_results$V0_clean,function(v0) list(v0=v0,v1=NULL))
  }
  if(filter_candidate_motifs_results$diss=='d1_L2'){
    alpha=1
    use0=FALSE
    use1=TRUE
    Y=lapply(filter_candidate_motifs_results$Y1,function(y1) list(y0=NULL,y1=y1))
    V=lapply(filter_candidate_motifs_results$V1_clean,function(v1) list(v0=NULL,v1=v1))
  }
  if(filter_candidate_motifs_results$diss=='d0_d1_L2'){
    alpha=filter_candidate_motifs_results$alpha
    use0=TRUE
    use1=TRUE
    Y=mapply(function(y0,y1) list(y0=y0,y1=y1),filter_candidate_motifs_results$Y0,filter_candidate_motifs_results$Y1,SIMPLIFY=FALSE)
    V=mapply(function(v0,v1) list(v0=v0,v1=v1),filter_candidate_motifs_results$V0_clean,filter_candidate_motifs_results$V1_clean,SIMPLIFY=FALSE)
  }
  w=filter_candidate_motifs_results$w
  max_gap=filter_candidate_motifs_results$max_gap
  d=ncol(filter_candidate_motifs_results$Y0[[1]])
  N=nrow(filter_candidate_motifs_results$D)
  K=ncol(filter_candidate_motifs_results$D)
  V_dom=lapply(filter_candidate_motifs_results$V0,function(v) rowSums(!is.na(v))!=0)
  V_length=unlist(lapply(V_dom,length))
  
  
  ### compute distances between motifs ####################################################################
  if(length(V)==1){
    VV_D=as.matrix(0)
    VV_S=as.matrix(1)
  }else{
    VV=combn(V,2,simplify=FALSE)
    VV=array(unlist(VV,recursive=FALSE),dim=c(2,length(VV)))
    VV_lengths=as.matrix(combn(V_length,2))
    VV_motif_overlap=floor(apply(VV_lengths,2,min)*motif_overlap)
    SD=mapply(.find_min_diss,VV[1,],VV[2,],VV_motif_overlap,
              MoreArgs=list(alpha=alpha,w=w,d=d,use0=use0,use1=use1),SIMPLIFY=TRUE)
    VV_D=matrix(0,nrow=length(V),ncol=length(V))
    VV_D[lower.tri(VV_D)]=SD[2,]
    VV_D=VV_D+t(VV_D) # matrix of distances
    VV_S=matrix(0,nrow=length(V),ncol=length(V))
    VV_S[lower.tri(VV_S)]=SD[1,]
    VV_S=VV_S+t(VV_S)+diag(1,nrow=length(V)) # matrix of shifts
  }
  
  ### determine a global radius Rall ######################################################################
  dataR=data.frame(P=as.vector(filter_candidate_motifs_results$P_clean),D=as.vector(filter_candidate_motifs_results$D_clean))
  if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
    R_all=max(dataR$D[dataR$P==1])
  }else{
    D_new=seq(0,max(filter_candidate_motifs_results$D_clean),length.out=10000)
    pred_knn=class::knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=k_knn,prob=TRUE)
    R_all=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=votes_knn_Rall)[1]]
  }
  
  ### cluster motifs and determine a group-specific radius Rm ##############################################
  if(length(V)==1){
    hclust_res=NULL
    R_m=R_all
  }else{
    VV_dist=as.dist(VV_D)
    hclust_res=hclust(VV_dist,method='average') # hierarchical clustering based on motif-motif distances
    V_hclust=cutree(hclust_res,h=2*R_all) # cut at high 2*R_all
    n_hclust=max(V_hclust)
    
    R_m=rep(NA,n_hclust)
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_length_i=V_length[index_i]
      V_P_i=as.matrix(filter_candidate_motifs_results$P_clean[,index_i])
      V_D_i=as.matrix(filter_candidate_motifs_results$D_clean[,index_i])
      
      dataR=data.frame(P=as.vector(V_P_i),D=as.vector(V_D_i))
      if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
        R_m[i_hclust]=max(dataR$D[dataR$P==1])
      }else{
        D_new=seq(0,max(V_D_i),length.out=10000)
        pred_knn=class::knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=k_knn,prob=TRUE)
        R_m[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=votes_knn_Rm)[1]]
      }
    }
  }
  
  ### output ###############################################################################################
  return(c(list(VV_D=VV_D,VV_S=VV_S,k_knn=k_knn,votes_knn_Rall=votes_knn_Rall,R_all=R_all,hclust_res=hclust_res,votes_knn_Rm=votes_knn_Rm,R_m=R_m),filter_candidate_motifs_results)) 
}
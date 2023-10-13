library(microbenchmark)
# load output_prof_elongation_new_data

source("../ProbKMA-FMD_functions.r")
probKMA_results <- output_prof_elongation_new_data

output <- probKMA_silhouette_rcpp(probKMA_results,TRUE,.domain,.select_domain,.diss_d0_d1_L2)

probKMA_silhouette <- function(probKMA_results,align=TRUE){
  # Compute the adapted silhouette index on the results of probKMA.
  # probKMA_results: output of probKMA function (with return_options=TRUE).
  # align: if TRUE, try all possible alignements between pieces of curves (corresponding to the same or to different motifs).
  # plot: if TRUE, the silhouette plot is drawn.

  ### compute silhouette #####################################################################################
  if(probKMA_results$diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(probKMA_results$Y0,function(y0) list(y0=y0,y1=NULL))
  }
  if(probKMA_results$diss=='d1_L2'){
    alpha=1
    use0=FALSE
    use1=TRUE
    Y=lapply(probKMA_results$Y1,function(y1) list(y0=NULL,y1=y1))
  }
  if(probKMA_results$diss=='d0_d1_L2'){
    alpha=probKMA_results$alpha
    use0=TRUE
    use1=TRUE
    Y=mapply(function(y0,y1) list(y0=y0,y1=y1),probKMA_results$Y0,probKMA_results$Y1,SIMPLIFY=FALSE)
  }
  w=probKMA_results$w
  d=ncol(probKMA_results$Y0[[1]])
  N=nrow(probKMA_results$P_clean)
  K=ncol(probKMA_results$P_clean)
  V_dom=lapply(probKMA_results$V0,function(v) rowSums(!is.na(v))!=0)
  V_length=unlist(lapply(V_dom,length))

  # extract pieces of curves that fall in the different motifs
  curves_in_motifs=apply(probKMA_results$P_clean,2,function(P_clean_k) which(P_clean_k==1))
  if(!is.null(ncol(curves_in_motifs))){
    curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs)))
  }
  
  curves_in_motifs_number=colSums(probKMA_results$P_clean)
  S_clean_k=mapply(function(s_k,curves_in_motif) s_k[curves_in_motif],
                   split(probKMA_results$S_clean,rep(seq_len(K),each=N)),curves_in_motifs,SIMPLIFY=FALSE)
  

  # compute distances between pieces of curves
  Y_in_motifs=unlist(mapply(function(curves_in_motif,s_k,v_dom){
    Y_in_motif=mapply(function(y,s){
      if(use0){
        d=ncol(y$y0) # dimension of curves
        y_len=nrow(y$y0)
        index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
        y$y0=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                   matrix(y$y0[index[index<=y_len],],ncol=d),
                   matrix(NA,nrow=sum(index>y_len),ncol=d))
        y$y0[!v_dom,]=NA
      }
      if(use1){
        d=ncol(y$y1) # dimension of curves
        y_len=nrow(y$y1)
        index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
        y$y1=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                   matrix(y$y1[index[index<=y_len],],ncol=d),
                   matrix(NA,nrow=sum(index>y_len),ncol=d))
        y$y1[!v_dom,]=NA
      }
      return(y)
    },Y[curves_in_motif],s_k,SIMPLIFY=FALSE)
  },curves_in_motifs,S_clean_k,V_dom,SIMPLIFY=FALSE),recursive=FALSE)
  Y_motifs=rep.int(1:K,curves_in_motifs_number)
  YY=combn(Y_in_motifs,2,simplify=FALSE)
  YY=array(unlist(YY,recursive=FALSE),dim=c(2,length(YY)))
  YY_lengths=combn(V_length[Y_motifs],2)
  swap=YY_lengths[1,]<YY_lengths[2,]
  YY[,swap]=rbind(YY[2,swap],YY[1,swap])
  if(align){
    # find distance between the two pieces of curves
    # no alignment for pieces corresponding to motifs with the same length
    # alignment for pieces corresponding to motifs with different lengths (requiring one piece inside the other)
    equal_length=(YY_lengths[1,]==YY_lengths[2,])
    SD=mapply(.find_diss,YY[1,],YY[2,],equal_length,
              MoreArgs=list(alpha=alpha,w=w,d,use0,use1),SIMPLIFY=TRUE)
  }else{
    # find minimum distance between the two pieces of curves, allowing alignment
    # minimum overlap required: minimum motif length
    cc_motifs=apply(combn(probKMA_results$c[Y_motifs],2),2,min)
    SD=mapply(.find_min_diss,YY[1,],YY[2,],cc_motifs,
              MoreArgs=list(alpha=alpha,w=w,d,use0,use1),SIMPLIFY=TRUE)
  }
  YY_D=matrix(0,nrow=length(Y_motifs),ncol=length(Y_motifs))
  YY_D[lower.tri(YY_D)]=SD[2,]
  YY_D=YY_D+t(YY_D)

  # compute intra-cluster distances
  intra=Reduce(cbind,lapply(Y_motifs,function(motif) Y_motifs==motif))
  diag(intra)=FALSE
  a=colSums(intra*YY_D)/(curves_in_motifs_number[Y_motifs]-1)

  # compute inter-cluster distances
  b_k=Reduce(rbind,lapply(1:(K-1),
                          function(k){
                            inter=Reduce(cbind,lapply(Y_motifs,
                                                      function(motif) Y_motifs==ifelse(motif+k>K,(motif+k)%%K,motif+k)))
                            b_k=colSums(inter*YY_D)/(curves_in_motifs_number[ifelse(Y_motifs+1>K,(Y_motifs+1)%%K,Y_motifs+1)])
                            return(b_k)}))
  if(is.matrix(b_k)){
    b=apply(b_k,2,min)
  }else{
    b=b_k
  }
  # compute silhouette
  silhouette=(b-a)/pmax(a,b)
  silhouette[is.nan(silhouette)]=0

  # compute average silhouette per cluster
  silhouette_average=rep(NA,K)
  for(k in seq_len(K)){
    silhouette_k=silhouette[Y_motifs==k]
    curves_in_motifs[[k]]=curves_in_motifs[[k]][order(silhouette_k,decreasing=TRUE)]
    silhouette[Y_motifs==k]=sort(silhouette_k,decreasing=TRUE)
    silhouette_average[k]=mean(silhouette_k)
  }


  return(list(silhouette=silhouette,motifs=Y_motifs,curves=unlist(curves_in_motifs),
              silhouette_average=silhouette_average))
}

output_prof <- probKMA_silhouette(probKMA_results,TRUE)
##setequal(output_prof,output)

results <- microbenchmark(probKMA_silhouette_rcpp(probKMA_results,FALSE,.domain,.select_domain,.diss_d0_d1_L2),probKMA_silhouette(probKMA_results,FALSE))
print(results)

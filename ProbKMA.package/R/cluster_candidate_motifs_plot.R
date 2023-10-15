
#' @export
cluster_candidate_motifs_plot <- function(cluster_candidate_motifs_results,ylab='',
                                          R_all=cluster_candidate_motifs_results$R_all,R_m=NULL,ask=TRUE){
  # Plot the results of cluster_candidate_motifs. 
  # cluster_candidate_motifs_results: output of cluster_candidate_motifs function.
  # R_all: global radius, used to cut the dendrogram (requiring groups to be more than 2*Rall apart).
  # R_m: vector with group-specific radii. The length of the vector must match the number of clusters
  #      obtained cutting the dendrogram at height 2*Rall. If NULL, Rm is determined in each group (based on distances between motifs of the same group and all curves).
  # ask: if TRUE the user is prompted before a new plot is drawn.
  
  ### prepare input data ##################################################################################
  if(cluster_candidate_motifs_results$diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(cluster_candidate_motifs_results$Y0,function(y0) list(y0=y0,y1=NULL))
    V=lapply(cluster_candidate_motifs_results$V0_clean,function(v0) list(v0=v0,v1=NULL))
  }
  if(cluster_candidate_motifs_results$diss=='d1_L2'){
    alpha=1
    use0=FALSE
    use1=TRUE
    Y=lapply(cluster_candidate_motifs_results$Y1,function(y1) list(y0=NULL,y1=y1))
    V=lapply(cluster_candidate_motifs_results$V1_clean,function(v1) list(v0=NULL,v1=v1))
  }
  if(cluster_candidate_motifs_results$diss=='d0_d1_L2'){
    alpha=cluster_candidate_motifs_results$alpha
    use0=TRUE
    use1=TRUE
    Y=mapply(function(y0,y1) list(y0=y0,y1=y1),cluster_candidate_motifs_results$Y0,cluster_candidate_motifs_results$Y1,SIMPLIFY=FALSE)
    V=mapply(function(v0,v1) list(v0=v0,v1=v1),cluster_candidate_motifs_results$V0_clean,cluster_candidate_motifs_results$V1_clean,SIMPLIFY=FALSE)
  }
  w=cluster_candidate_motifs_results$w
  max_gap=cluster_candidate_motifs_results$max_gap
  d=ncol(cluster_candidate_motifs_results$Y0[[1]])
  N=nrow(cluster_candidate_motifs_results$D)
  K=ncol(cluster_candidate_motifs_results$D)
  V_dom=lapply(cluster_candidate_motifs_results$V0,function(v) rowSums(!is.na(v))!=0)
  V_length=unlist(lapply(V_dom,length))
  
  ### plot global distances and global radius Rall ########################################################
  par(ask=ask,mfrow=c(2,1))
  dist_breaks=diff(hist(cluster_candidate_motifs_results$D_clean[cluster_candidate_motifs_results$D_clean<quantile(cluster_candidate_motifs_results$D_clean,0.90)],plot=FALSE)$breaks[1:2])
  breaks=seq(0,ceiling(max(cluster_candidate_motifs_results$D_clean)/10)*10+dist_breaks,by=dist_breaks)
  xlim=c(0,quantile(cluster_candidate_motifs_results$D_clean,0.95))
  
  # histogram of distances between all motifs and all curves
  hist_res=hist(cluster_candidate_motifs_results$D_clean,breaks=breaks,
                probability=TRUE,xlim=xlim,xlab='Distances motif-curve',main='Distances motif-curve - All motifs')
  D_clean_density=density(cluster_candidate_motifs_results$D_clean,bw=hist_res$breaks[2])
  lines(D_clean_density$x,D_clean_density$y,lwd=2)
  abline(v=R_all,lwd=2)
  text(x=R_all,y=max(hist_res$density),labels=expression(R[all]),pos=4)
  
  # histogram of distances for all motifs, separately for curves with the motifs and for curves without the motifs
  hist_res=hist(cluster_candidate_motifs_results$D_clean[cluster_candidate_motifs_results$P_clean==1],breaks=breaks,
                probability=TRUE,xlim=xlim,xlab='Distances motif-curve',main='Distances motif-curve - All motifs')
  D_clean_density1=density(cluster_candidate_motifs_results$D_clean[cluster_candidate_motifs_results$P_clean==1],bw=breaks[2])
  lines(D_clean_density1$x,D_clean_density1$y,col='red',lwd=2)
  hist(cluster_candidate_motifs_results$D_clean[cluster_candidate_motifs_results$P_clean==0],breaks=breaks,probability=TRUE,add=TRUE)
  D_clean_density0=density(cluster_candidate_motifs_results$D_clean[cluster_candidate_motifs_results$P_clean==0],bw=breaks[2])
  lines(D_clean_density0$x,D_clean_density0$y,col='blue',lwd=2)
  abline(v=R_all,lwd=2)
  text(x=R_all,y=max(hist_res$density),labels=expression(R[all]),pos=4)
  legend('topright',legend=c('Curves with motif','Curves without motif'),lwd=2,col=c('red','blue'))
  
  ### cut dendrogram and check group-specific radius Rm ##################################################
  if(is.null(cluster_candidate_motifs_results$hclust_res)){
    V_hclust=1
  }else{
    V_hclust=cutree(cluster_candidate_motifs_results$hclust_res,h=2*R_all) # cut at high 2*R_all
  }
  n_hclust=max(V_hclust)
  if(!is.null(cluster_candidate_motifs_results$hclust_res)){
    par(mfrow=c(1,1))
    dendr=as.dendrogram(cluster_candidate_motifs_results$hclust_res,hang=1)
    labels_cex(dendr)=0.8
    plot(dendr,ylab='Distance motif-motif',main='Dendrogram of motifs')
    abline(h=2*R_all,col='black',lwd=2)
    text(x=1,y=2*R_all,labels=expression('2R'[all]),pos=3,offset=0.5)
    if(n_hclust>1)
      rect.dendrogram(dendr,k=n_hclust)
  }
  
  # re-compute or check group-specific radius Rm, and plot
  if(is.null(R_m))
    R_m=rep(NA,n_hclust)
  if(length(R_m)!=n_hclust)
    stop(paste0('The length of the vector R_m must match the number of clusters: ',n_hclust))
  
  for(i_hclust in seq_len(n_hclust)){
    index_i=which(V_hclust==i_hclust) # motifs in cluster i
    V_length_i=V_length[index_i]
    V_P_i=as.matrix(cluster_candidate_motifs_results$P_clean[,index_i])
    V_D_i=as.matrix(cluster_candidate_motifs_results$D_clean[,index_i])
    
    if(is.na(R_m[i_hclust])){
      dataR=data.frame(P=as.vector(V_P_i),D=as.vector(V_D_i))
      if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
        R_m[i_hclust]=max(dataR$D[dataR$P==1])
      }else{
        D_new=seq(0,max(V_D_i),length.out=10000)
        pred_knn=knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=cluster_candidate_motifs_results$k_knn,prob=TRUE)
        R_m[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=cluster_candidate_motifs_results$votes_knn_Rm)[1]]
      }
    }
    
    # plot histograms and motifs
    par(mfrow=c(2,1))
    
    # histogram of distances between motifs in cluster i and all curves
    hist_res=hist(V_D_i,breaks=breaks,probability=TRUE,xlim=xlim,xlab='Distances motif-curve',main=paste0('Distances motif-curve - Cluster ',i_hclust))
    D_clean_density=density(V_D_i,bw=hist_res$breaks[2])
    lines(D_clean_density$x,D_clean_density$y,lwd=2)
    abline(v=R_m[i_hclust],lwd=2)
    text(x=R_m[i_hclust],y=max(hist_res$density),labels=bquote(R[.(i_hclust)]),pos=4)
    
    # histogram of distances for motifs in cluster i, separately for curves with the motifs and for curves without the motifs
    hist_res=hist(V_D_i[V_P_i==1],breaks=breaks,probability=TRUE,xlim=xlim,xlab='Distances motif-curve',main=paste0('Distances motif-curve - Cluster ',i_hclust))
    D_clean_density1=density(V_D_i[V_P_i==1],bw=breaks[2])
    lines(D_clean_density1$x,D_clean_density1$y,col='red',lwd=2)
    hist(V_D_i[V_P_i==0],breaks=breaks,probability=TRUE,add=TRUE)
    D_clean_density0=density(V_D_i[V_P_i==0],bw=breaks[2])
    lines(D_clean_density0$x,D_clean_density0$y,col='blue',lwd=2)
    abline(v=R_m[i_hclust],lwd=2)
    text(x=R_m[i_hclust],y=max(hist_res$density),labels=bquote(R[.(i_hclust)]),pos=4)
    legend('topright',legend=c('Curves with motif','Curves without motif'),lwd=2,col=c('red','blue'))
    
    par(mfrow=c(2+d,1))
    
    # plot aligned motifs
    VV_S_i=as.matrix(cluster_candidate_motifs_results$VV_S[index_i,index_i])
    V0_i=cluster_candidate_motifs_results$V0_clean[index_i]
    for(jd in seq_len(d)){
      V0_i_jd=lapply(V0_i,function(v0_i) as.matrix(v0_i[,jd]))
      plot(1,1,ylim=range(unlist(V0_i_jd)),xlim=c(min(VV_S_i[1,]),max(VV_S_i[1,]+V_length_i-1)),ylab=ylab[jd],
           xlab='Position',main=paste0('Motifs in Cluster ',i_hclust),type='n',cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
      for(j in seq_along(V_length_i))
        lines(seq(VV_S_i[1,j],VV_S_i[1,j]+V_length_i[j]-1),V0_i_jd[[j]],col=j)
    }
    
    # plot alignment
    plot(1,1,ylim=c(1,length(V_length_i)),xlim=c(min(VV_S_i[1,]),max(VV_S_i[1,]+V_length_i-1)),
         ylab='Motifs',xlab='Position',main=paste0('Alignment motifs in Cluster ',i_hclust),type='n',yaxt='n',cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
    for(j in seq_along(index_i)){
      lines(c(VV_S_i[1,j],VV_S_i[1,j]+V_length_i[j]-1),c(j,j),lwd=2,col=j)
      text(VV_S_i[1,j],j,labels=index_i[j],pos=2)
    }
    
    # plot
    V_mean_diss_approx_i=apply(V_D_i,2,function(x) mean(x[x<=R_m[i_hclust]])) # approximate average distance (actual radius)
    V_frequencies_approx_i=colSums(V_D_i<=R_m[i_hclust]) # approximate frequency 
    V_frequencies_approx_i_scatter=V_frequencies_approx_i+rnorm(length(V_frequencies_approx_i),sd=0.1)
    plot(V_mean_diss_approx_i,V_frequencies_approx_i_scatter,cex=V_length_i/min(V_length_i),col=seq_along(V_length_i),
         xlab='Approx average distance',ylab='Approx frequency',main=paste0('Motif in Cluster ',i_hclust),yaxt='n',cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
    axis(2,at=seq(min(V_frequencies_approx_i),max(V_frequencies_approx_i)),cex.axis=1.5,cex.lab=1.5)
    text(V_mean_diss_approx_i,V_frequencies_approx_i_scatter,labels=index_i,pos=4)
  }
  return()
}
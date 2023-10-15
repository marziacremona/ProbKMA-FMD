
probKMA_silhouette <- function(probKMA_results,align=FALSE,plot=TRUE){
  # Compute the adapted silhouette index on the results of probKMA.
  # probKMA_results: output of probKMA function (with return_options=TRUE).
  # align: if TRUE, try all possible alignements between pieces of curves (corresponding to the same or to different motifs).
  # plot: if TRUE, the silhouette plot is drawn.
 
  result <- .probKMA_silhouette_rcpp(probKMA_results,
                                     domain,select_domain,
                                     diss_d0_d1_L2,align) 
  
  
  ### plot silhouette ########################################################################################
  if(plot){
    K=ncol(probKMA_results$P_clean)
    n=length(result$silhouette)
    sil=rev(result$silhouette)
    y=barplot(sil,space=c(0,rev(diff(result$Y_motifs))),xlab='Silhouette index',names='',
              xlim=c(min(0,min(sil)-0.05),1.2),horiz=TRUE,las=1,mgp=c(2.5,1,0),col='gray')
    text(ifelse(sil>=0,sil+0.03,sil-0.03),ifelse(sil>0,y,y+0.2),labels=rev(unlist(result$curves_in_motifs)),cex=0.5)
    title(main='Silhouette plot')
    title(sub=paste("Average silhouette width:",round(mean(result$silhouette_average),2)),adj=0)
    mtext(paste("n =",n),adj=0)
    mtext(substitute(K~~"motifs",list(K=K)),adj=1)
    mtext(expression(paste(j," : ",n[j]," | avg ",s[i])),adj=1.04,line=-1.2)
    y=rev(y)
    for(k in seq_len(K)){
      y_k=mean(y[result$Y_motifs==k])
      text(1,y_k,
           paste0(k,": ",result$curves_in_motifs_number[k]," | ",format(result$silhouette_average[k],digits=1,nsmall=2)),xpd=NA,adj=0.1)
    }
  }
  
  return(list(silhouette=silhouette,motifs=Y_motifs,curves=unlist(curves_in_motifs),
              silhouette_average=silhouette_average))
}
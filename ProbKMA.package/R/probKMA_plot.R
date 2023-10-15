
probKMA_plot <- function(probKMA_results,ylab='',cleaned=FALSE){
  # Plot the results of probKMA.
  # probKMA_results: output of probKMA function.
  # ylab: a vector of length d, with the titles for the y axis for each dimension.
  # cleaned: if TRUE, plot the cleaned motifs.
  
  d=ncol(probKMA_results$Y0[[1]])
  N=nrow(probKMA_results$P)
  K=ncol(probKMA_results$P)
  V_dom=lapply(probKMA_results$V0,function(v) rowSums(!is.na(v))!=0)
  S_k=split(probKMA_results$S,rep(seq_len(K),each=N))
  P_k=split(probKMA_results$P,rep(seq_len(K),each=N))
  
  ### plot motifs with matched curves ########################################################################
  if(cleaned){
    S_clean_k=split(probKMA_results$S_clean,rep(seq_len(K),each=N))
    P_clean_k=split(probKMA_results$P_clean,rep(seq_len(K),each=N))
    if(is.null(probKMA_results$V1[[1]])){
      mapply(function(v,v_dom,s_k,p_clean_k,k){
        keep=which(p_clean_k==1)
        Y_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          probKMA_results$Y0[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
    }else{
      mapply(function(v0,v1,v_dom,s_k,p_clean_k,k){
        keep=which(p_clean_k==1)
        Y0_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          probKMA_results$Y0[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          probKMA_results$Y1[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y0_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v0[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y1_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j],' derivative'))
                 points(v1[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0_clean,probKMA_results$V1_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
    }
  }else{
    if(is.null(probKMA_results$V1[[1]])){
      mapply(function(v,v_dom,s_k,p_k,k){
        Y_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          probKMA_results$Y0,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0,V_dom,S_k,P_k,seq_len(K))
    }else{
      mapply(function(v0,v1,v_dom,s_k,p_k,k){
        Y0_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          probKMA_results$Y0,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          probKMA_results$Y1,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y0_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v0[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y1_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j],' derivative'))
                 points(v1[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0,probKMA_results$V1,V_dom,S_k,P_k,seq_len(K))
    }
  }
  layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
  
  ### plot motifs ############################################################################################
  if(cleaned){
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             motif_length=unlist(lapply(probKMA_results$V0_clean,nrow))
             plot(probKMA_results$V0_clean[[1]][,j],type='l',col=2,lwd=7,lty=1,main=ylab[j],xlim=c(1,max(motif_length)),
                  ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V0_clean),na.rm=TRUE),max(unlist(probKMA_results$V0_clean),na.rm=TRUE)))
             mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                    probKMA_results$V0_clean[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             plot.new()
             legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
             return()})
  }else{
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             motif_length=unlist(lapply(probKMA_results$V0,nrow))
             plot(probKMA_results$V0[[1]][,j],type='l',col=2,lwd=7,lty=1,main=ylab[j],xlim=c(1,max(motif_length)),
                  ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V0),na.rm=TRUE),max(unlist(probKMA_results$V0),na.rm=TRUE)))
             mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                    probKMA_results$V0[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             plot.new()
             legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
             return()})
  }
  if(!is.null(probKMA_results$V1[[1]])){
    layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
    if(cleaned){
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               motif_length=unlist(lapply(probKMA_results$V1_clean,nrow))
               plot(probKMA_results$V1_clean[[1]][,j],type='l',col=2,lwd=7,lty=1,main=paste0(ylab[j],' derivative'),xlim=c(1,max(motif_length)),
                    ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V1_clean),na.rm=TRUE),max(unlist(probKMA_results$V1_clean),na.rm=TRUE)))
               mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                      probKMA_results$V1_clean[-1],seq_len(K-1))
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
               return()})
    }else{
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               motif_length=unlist(lapply(probKMA_results$V1,nrow))
               plot(probKMA_results$V1[[1]][,j],type='l',col=2,lwd=7,lty=1,main=paste0(ylab[j],' derivative'),xlim=c(1,max(motif_length)),
                    ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V1),na.rm=TRUE),max(unlist(probKMA_results$V1),na.rm=TRUE)))
               mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                      probKMA_results$V1[-1],seq_len(K-1))
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
               return()})
    }
  }
  
  ### plot memberships #######################################################################################
  par(mfrow=c(K,1),mar=c(3,4,4,2)+0.1)
  if(cleaned){
    mapply(function(p_k,p_clean_k,k){
      col=rep('lightgray',N)
      col[p_clean_k==1]='gray35'
      barplot(p_k,names.arg=seq_len(N),col=col,las=2,ylim=c(0,1),ylab='memberships',main=paste('Motif',k))
    },P_k,P_clean_k,seq_len(K))
  }else{
    mapply(function(p_k,k){
      barplot(p_k,names.arg=seq_len(N),col='gray',las=2,ylim=c(0,1),ylab='memberships',main=paste('Motif',k))
    },P_k,seq_len(K))
  }
  
  ### plot objective function and Bhattacharyya distance #####################################################
  par(mfrow=c(1,1))
  plot(seq_len(probKMA_results$iter),probKMA_results$J_iter,type='l',xlab='iter',ylab='objective function',main='Objective function Jm')
  plot(seq_len(probKMA_results$iter),probKMA_results$BC_dist_iter,type='l',xlab='iter',ylab='distance between memberships',main='Bhattacharyya distance between memberships')
  
  return()
}
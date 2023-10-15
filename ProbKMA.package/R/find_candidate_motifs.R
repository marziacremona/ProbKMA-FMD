
find_candidate_motifs <- function(Y0,Y1=NULL,K,c,n_init=10,name='results',names_var='',
                                  probKMA_options=NULL,silhouette_align=FALSE,plot=TRUE,worker_number=NULL){
  # Run multiple times probKMA function with different K,c and initializations, 
  # with the aim to find a set of candidate motifs.
  # If the folder name_KK_cc is already present and n result files are already present, 
  # load them and continue with the n_init-n runs.
  # Y0: list of N vectors, for univariate curves y_i(x), or
  #     list of N matrices with d columns, for d-dimensional curves y_i(x), 
  #     with the evaluation of curves (all curves should be evaluated on a uniform grid).
  #     When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
  # Y1: list of N vectors, for univariate derivative curves y'_i(x), or
  #     list of N matrices with d columns, for d-dimensional derivatibe curves y'_i(x),
  #     with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
  #     When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions. 
  #     Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
  # K: vector with numbers of motifs that must be tested.
  # c: vector with minimum motifs lengths that must be tested.
  # n_init: number of random initialization for each combination of K and c.
  # name: name of the folders when the results are saved.
  # names_var: vector of length d, with names of the variables in the different dimensions.
  # probKMA_options: list with options for probKMA (see the help of probKMA).
  # plot: if TRUE, summary plots are drawn. 
  # worker_number: number of CPU cores to be used for parallelization (default number of CPU cores -1). If worker_number=1, the function is run sequentially. 
  
  ### check input #############################################################################################
  # check required input
  if(missing(K))
    stop('K must be specified')
  if(missing(c))
    stop('c must be specified')
  # check K
  if(!is.vector(K))
    stop('K should be a vector.')
  # check c
  if(!is.vector(c))
    stop('c should be a vector.')
  # check n_init
  if(n_init%%1!=0)
    stop('Number of initializations n_init should be an integer.')
  if(n_init<1)
    stop('Number of initializations n_init should be at least 1.')
  # check name
  if(!is.character(name))
    stop('name should be a string.')
  if(grepl(' ',name)){
    warning('name contains spaces. Changing spacing in underscores.')
    name=gsub(" ","_",name,fixed=TRUE)
  }
  # check Y0
  if(missing(Y0))
    stop('Y0 must be specified.')
  if(!is.null(Y0))
    if(!is.list(Y0))
      stop('Y0 should be a list of vectors or matrices.')
  if((FALSE %in% lapply(Y0,is.matrix))&&(FALSE %in% lapply(Y0,is.vector)))
    stop('Y0 should be a list of vectors or matrices.')
  N=length(Y0) # number of curves
  if(N<5)
    stop('More curves y_i(x) needed.')
  Y0=lapply(Y0,as.matrix)
  
  # set return_options=TRUE
  if(!is.null(probKMA_options$return_options)){
    if(probKMA_options$return_options==FALSE){
      warning('Setting return_option=TRUE')
      probKMA_options$return_options=TRUE
    }
  }
  
  ### set parallel jobs #############################################################################
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
  len_mean=mean(unlist(lapply(Y0,nrow)))
  c_mean=mean(c)
  if(((len_mean/c_mean>10)&(length(K)*length(c)*n_init<40))|(length(K)*length(c)*n_init==1)){
    if(is.null(probKMA_options$worker_number))
      probKMA_options$worker_number=worker_number
    cl_find=NULL
  }else{
    probKMA_options$worker_number=1
    if(worker_number>1){
      cl_find=parallel::makeCluster(worker_number,timeout=60*60*24*30)
      parallel::clusterExport(cl_find,c('name','names_var','Y0','Y1','probKMA_options',
                              'probKMA','probKMA_plot','probKMA_silhouette','compute_motif', # delete these in the package
                              'mapply_custom','diss_d0_d1_L2','domain','select_domain','find_min_diss','compute_Jk'),envir=environment()) 
      parallel::clusterCall(cl_find,function()library(parallel,combinat))
      on.exit(parallel::stopCluster(cl_find))
    }else{
      cl_find=NULL
    }
  }
  
  ### run probKMA ##########################################################################################
  i_c_K=expand.grid(seq_len(n_init),c,K)
  results=mapply_custom(cl_find,function(K,c,i){
    dir.create(paste0(name,"_K",K,"_c",c),showWarnings=FALSE)
    files=list.files(paste0(name,"_K",K,"_c",c))
    message("K",K,"_c",c,'_random',i)
    if(paste0('random',i,'.RData') %in% files){
      load(paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
      return(list(probKMA_results=probKMA_results,
                  time=time,silhouette=silhouette))
    }else{
      iter=iter_max=1
      while(iter==iter_max){
        start=proc.time()
        probKMA_results=do.call(probKMA,c(list(Y0=Y0,Y1=Y1,K=K,c=c),probKMA_options))
        end=proc.time()
        time=end-start
        iter=probKMA_results$iter
        iter_max=probKMA_results$iter_max
        if(iter==iter_max)
          warning('Maximum number of iteration reached. Re-starting.')
      }
      pdf(paste0(name,"_K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
      probKMA_plot(probKMA_results,ylab=names_var,cleaned=FALSE)
      dev.off()
      pdf(paste0(name,"_K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
      probKMA_plot(probKMA_results,ylab=names_var,cleaned=TRUE)
      dev.off()
      pdf(paste0(name,"_K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
      silhouette=probKMA_silhouette(probKMA_results,align=silhouette_align,plot=TRUE)
      dev.off()
      save(probKMA_results,time,silhouette,
           file=paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
      return(list(probKMA_results=probKMA_results,
                  time=time,silhouette=silhouette))
    }
  },i_c_K[,3],i_c_K[,2],i_c_K[,1],SIMPLIFY=FALSE)
  results=split(results,list(factor(i_c_K[,2],c),factor(i_c_K[,3],K)))
  results=split(results,rep(K,each=length(c)))
  
  ### plot silhouette average #################################################################################
  silhouette_average_sd=lapply(results,
                               function(results){
                                 silhouette_average=lapply(results,
                                                           function(results){
                                                             silhouette_average=numeric(n_init)
                                                             silhouette_sd=numeric(n_init)
                                                             for(i in seq_len(n_init)){
                                                               silhouette_average[i]=mean(results[[i]]$silhouette$silhouette_average)
                                                               silhouette_sd[i]=sd(results[[i]]$silhouette$silhouette_average)
                                                             }
                                                             return(cbind(silhouette_average,silhouette_sd))
                                                           })
                               })
  if(plot){
    pdf(paste0(name,'_silhouette.pdf'),width=7,height=5)
    for(i in seq_along(K)){
      silhouette_average_plot=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) average_sd[order(average_sd[,1],decreasing=TRUE),1])),ncol=length(c))
      silhouette_sd_plot=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) average_sd[order(average_sd[,1],decreasing=TRUE),2])),ncol=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      matplot(matrix(rep(1:n_init,length(c))+rep(seq(-0.1,0.1,length.out=length(c)),each=n_init),ncol=length(c)),
              silhouette_average_plot,type='b',pch=16,lty=1,col=1+seq_along(c),ylim=c(0,1),xaxt='n',
              xlab='',ylab='Silhouette average',main=paste0('K=',K[i]))
      shift=seq(-0.1,0.1,length.out=length(c))
      for(ii in seq_along(c)){
        segments(1:n_init+shift[ii],silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
                 1:n_init+shift[ii],silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],col=ii+1)
        segments(1:n_init+shift[ii]-0.1,silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],
                 1:n_init+shift[ii]+0.1,silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],col=ii+1)
        segments(1:n_init+shift[ii]-0.1,silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
                 1:n_init+shift[ii]+0.1,silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],col=ii+1)
        text(1:n_init+shift[ii],silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
             silhouette_order[,ii],col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),lty=1)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
    }
    dev.off()
  }
  
  ### plot processing time ####################################################################################
  times=lapply(results,
               function(results){
                 times=lapply(results,
                              function(results)
                                times=unlist(lapply(results,function(results) results$time[3])))
               })
  if(plot){
    pdf(paste0(name,'_times.pdf'),width=7,height=5)
    y_max=max(unlist(times))*1.2
    times_plot=vector('list',length(K))
    for(i in seq_along(K)){
      times_plot[[i]]=Reduce(cbind,
                             lapply(seq_along(c),
                                    function(j)
                                      as.matrix(times[[i]][[j]][order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)]))
      )
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      matplot(matrix(rep(1:n_init,length(c))+rep(seq(-0.1,0.1,length.out=length(c)),each=n_init),ncol=length(c)),
              times_plot[[i]],type='b',pch=16,lty=1,col=1+seq_along(c),ylim=c(0,y_max),xaxt='n',
              xlab='',ylab='Time',main=paste0('K=',K[i]))
      shift=seq(-0.1,0.1,length.out=length(c))
      for(ii in seq_along(c)){
        text(1:n_init+shift[ii],times_plot[[i]][,ii],silhouette_order[,ii],col=ii+1,pos=1,cex=0.7)
      }
      legend('topleft',legend=paste0('c=',c),col=1+seq_along(c),lty=1)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
    }
    for(i in seq_along(K)){
      boxplot(times_plot[[i]],col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
              xlab='',ylab='Times',main=paste0('K=',K[i]))
    }
    dev.off()
  }
  
  ### plot dissimilarities ####################################################################################
  if(plot){
    D=lapply(results,
             function(results){
               D=lapply(results,
                        function(results){
                          D=lapply(results,
                                   function(results){
                                     D=as.vector(results$probKMA_results$D)
                                   })
                        })
             })
    pdf(paste0(name,'_dissimilarities.pdf'),width=7,height=5)
    y_max=max(unlist(D))
    for(i in seq_along(K)){
      D_plot=matrix(unlist(D[[i]]),ncol=length(c))
      boxplot(D_plot,col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
              xlab='',ylab='Dissimilarity',main=paste0('K=',K[i]))
    }
    for(i in seq_along(K)){
      for(j in seq_along(c)){
        silhouette_order=order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)
        D_plot=D[[i]][[j]][silhouette_order]
        boxplot(D_plot,col=j+1,ylim=c(0,y_max),names=silhouette_order,
                xaxt='n',ylab='Dissimilarity',main=paste0('K=',K[i],'   c=',c[j]))
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
      }
    }
    dev.off()
    
    D_clean=lapply(results,
                   function(results){
                     D_clean=lapply(results,
                                    function(results){
                                      D_clean=lapply(results,
                                                     function(results){
                                                       D_clean=as.vector(results$probKMA_results$D_clean)
                                                     })
                                    })
                   })
    pdf(paste0(name,'_dissimilarities_clean.pdf'),width=7,height=5)
    y_max=max(unlist(D_clean))
    for(i in seq_along(K)){
      D_plot=matrix(unlist(D_clean[[i]]),ncol=length(c))
      boxplot(D_plot,col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
              xlab='',ylab='Dissimilarity',main=paste0('K=',K[i]))
    }
    for(i in seq_along(K)){
      for(j in seq_along(c)){
        silhouette_order=order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)
        D_plot=D_clean[[i]][[j]][order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)]
        boxplot(D_plot,col=j+1,ylim=c(0,y_max),names=silhouette_order,
                xaxt='n',ylab='Dissimilarity',main=paste0('K=',K[i],'   c=',c[j]))
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
      }
    }
    dev.off()
  }
  
  ### plot motif lengths ######################################################################################
  if(plot){
    motif_length=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0,nrow))
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    pdf(paste0(name,'_lengths.pdf'),width=7,height=5)
    motif_length_plot=lapply(motif_length,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths',main=paste0('K=',K[i]))
      abline(h=c,col=1+seq_along(c))
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_length[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
    
    motif_clean_length=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0_clean,nrow))
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    pdf(paste0(name,'_lengths_clean.pdf'),width=7,height=5)
    motif_length_plot=lapply(motif_clean_length,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths',main=paste0('K=',K[i]))
      abline(h=c,col=1+seq_along(c))
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_clean_length[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
    
    pdf(paste0(name,'_lengths_perc.pdf'),width=7,height=5)
    motif_length_perc=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0,nrow))/results$probKMA_results$c*100
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    motif_length_plot=lapply(motif_length_perc,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths (% of minimum length)',main=paste0('K=',K[i]))
      abline(h=100)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_length_perc[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
    
    pdf(paste0(name,'_lengths_clean_perc.pdf'),width=7,height=5)
    motif_clean_length_perc=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0_clean,nrow))/results$probKMA_results$c*100
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    motif_length_plot=lapply(motif_clean_length_perc,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=10^-1),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths (% of minimum length)',main=paste0('K=',K[i]))
      abline(h=100)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_clean_length_perc[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
  }
  
  ### output ##################################################################################################
  return(list(name=name,K=K,c=c,n_init=n_init,silhouette_average_sd=silhouette_average_sd,times=times))
}
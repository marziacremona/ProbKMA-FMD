# ProbKMA: Probabilistic k-mean with local alignment
# ProbKMA-FMD: ProbKMA-based Functional Motif Discovery

library(combinat)
library(parallel)
library(class)
library(dendextend)


# Use mapply if paralization is not needed(cluster is not available)
# Use clusterMap if paralization is needed
.mapply_custom <- function(cl,FUN,...,MoreArgs=NULL,SIMPLIFY=TRUE,USE.NAMES=TRUE){
  if(is.null(cl)){
    mapply(FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }else{
    clusterMap(cl,FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }
}
.diss_d0_d1_L2 <- function(y,v,w,alpha){
  # Dissimilarity index for multidimensional curves (dimension=d).
  # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
  # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
  
  .diss_L2 <- function(y,v,w){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # L2 distance with normalization on common support.
    # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    
    sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
  }
  
  if(alpha==0){
    return(.diss_L2(y[[1]],v[[1]],w))
  }else if(alpha==1){
    return(.diss_L2(y[[2]],v[[2]],w))
  }else{
    return((1-alpha)*.diss_L2(y[[1]],v[[1]],w)+alpha*.diss_L2(y[[2]],v[[2]],w))
  }
}
.domain <- function(v,use0){# determine the domain (or support) of a curve or its derivative based on the boolean argument use0.
  if(use0){
    #is.na returns a matrix whose elements are TRUE or FALSE wheter it is Na or not Na 
    rowSums(!is.na(v[[1]]))!=0 
  }else{
    rowSums(!is.na(v[[2]]))!=0
  } # return true if there exists at least one element different from Na in the row.
}
# v[[1]] accesses to the center cluster(matrix)
.select_domain <- function(v,v_dom,use0,use1){ # select the portion of the domain of the curve based on .domain
  if(use0)
    v[[1]]=as.matrix(v[[1]][v_dom,])
  if(use1)
    v[[2]]=as.matrix(v[[2]][v_dom,])
  return(v)
}
 # the find_min_diss function performs a search for the best shift (warping) value between 
 # two multidimensional curves, y and v, by evaluating dissimilarity scores for different shifts. 
 # It then returns the shift that minimizes the dissimilarity and the corresponding dissimilarity value.
.find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1){
  # y: list of two elements y0=y(x), y1=y'(x), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x), matrices with d columns.
  # alpha: weight coefficient between d0.L2 and d1.L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  
  v_dom=.domain(v,use0) # domain of the v curve
  v=.select_domain(v,v_dom,use0,use1) # select a portion of the curve
  v_len=length(v_dom)
  y_len=unlist(lapply(y,nrow))[1] # number of rows of the curve 
  s_rep=(1-(v_len-c_k)):(y_len-v_len+1+(v_len-c_k)) #sequence of integers

# The function is designed to create a list y_rep containing different subsets of the original y curve, 
# each shifted based on the current shift index i.
  y_rep=lapply(s_rep,  # the function apply function(i) for each element i in the sequence s_rep, y_rep è una lista fatta di y_rep_i
               function(i){ 
                 index=i-1+seq_len(v_len) #This line of code calculates the indices used to extract subsets of y based on the current shift index i
                 y_rep_i=list(y0=NULL,y1=NULL)# initialization of an empty list 
                 if(use0){
                   y_rep_i$y0=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[1]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }# Construct a matrix y0 combaning three parts: a matrix of NA values with the same number
                # of rows as the number of negative or zero values in index, a submatrix of the curve y[[1]]
                # of the list y using positive and non-zero values from index, and another matrix of NA  
                # values with the same number of rows as the number of positive values exceeding y_len in index. 
                 if(use1){
                   y_rep_i$y1=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[2]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 y_rep_i=.select_domain(y_rep_i,v_dom,use0,use1)
                 return(y_rep_i)
               })
                
  length_inter=unlist(lapply(y_rep,
                             function(y_rep_i){
                               if(use0)
                                 return(sum((!is.na(y_rep_i$y0[,1]))))
                               return(sum((!is.na(y_rep_i$y1[,1]))))
                             }))#This line of code calculates the length of the intersection between the modified subsets of y (obtained from y_rep) and the v curve.
  valid=length_inter>=c_k # This line of code generates a logical vector valid, where each element is TRUE if the length of the intersection between the modified subset of y and v is greater than or equal to c_k
  if(sum(valid)==0){ # meaning there is no intersection of sufficient length
    valid[length_inter==max(length_inter)]=TRUE # the code ensures that at least the subset with the longest intersection is marked as valid
  }
  s_rep=s_rep[valid] 
  y_rep=y_rep[valid]
  # This line calculates the dissimilarity values for each element in y_rep using the function .diss_d0_d1_L2. 
  # The function .diss_d0_d1_L2 takes as input y_rep, v, w, and alpha
  d_rep=unlist(mapply(.diss_d0_d1_L2,y_rep,MoreArgs=list(v,w,alpha))) # mapply(FUN,... (the list for which we want to apply FUN multiple times), the other args of the FUN)
  return(c(s_rep[which.min(d_rep)],min(d_rep)))# returns a vector with two elements: the shift value that minimizes the dissimilarity, and the minimum dissimilarity value.
}

.find_diss <- function(y,v,alpha,w,aligned,d,use0,use1){
  # Find dissimilarity between multidimensional curves (dimension=d), without alignment unless their lengths are different.
  # Return shift and dissimilarity.
  # To be used by probKMA_silhouette function.
  # y: list of two elements y0=y(x), y1=y'(x), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x), matrices with d columns.
  # alpha: weight coefficient between d0.L2 and d1.L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # aligned: if TRUE, curves are already aligned. If FALSE, the shortest curve is aligned inside the longest.
  
  v_dom=.domain(v,use0)
  v=.select_domain(v,v_dom,use0,use1)
  v_len=length(v_dom)
  y_len=unlist(lapply(y,nrow))[1]
  if(aligned){
    s_rep=1
  }else{
    s_rep=1:(y_len-v_len+1)
  }
  y_rep=lapply(s_rep,
               function(i){
                 index=i-1+seq_len(v_len)
                 y_rep_i=list(y0=NULL,y1=NULL)
                 if(use0){
                   y_rep_i$y0=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[1]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 if(use1){
                   y_rep_i$y1=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[2]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 y_rep_i=.select_domain(y_rep_i,v_dom,use0,use1)
                 return(y_rep_i)
               })
  d_rep=unlist(mapply(.diss_d0_d1_L2,y_rep,MoreArgs=list(v,w,alpha)))
  return(c(s_rep[which.min(d_rep)],min(d_rep)))
}


.compute_motif <- function(v_dom,s_k,p_k,Y,m,use0,use1){ #compute the new center of the cluster, ha bisogno dei dati e delle probabilità di assegnazione e degli shift al passo precedente 
  # Compute the new motif v_new.
  # v_dom: TRUE for x in supp(v). #vettore di bool = TRUE se e solo se la corrispondente riga sta nel supporto
  # s_k: shift vector for motif k.
  # p_k: membership vector for motif k.
  # Y: list of N lists of two elements, Y0=y_i(x), Y1=y'_i(x), matrices with d columns, for d-dimensional curves.
  # Y contiene tutte le curve del dataset, quindi lista di N curve d-dimensionali (N matrici)

  .domain <- function(v,use0){
    if(use0){
      rowSums(!is.na(v[[1]]))!=0
    }else{
      rowSums(!is.na(v[[2]]))!=0
    }
  }
  .select_domain <- function(v,v_dom,use0,use1){
    if(use0)
      v[[1]]=as.matrix(v[[1]][v_dom,]) #select the portion of the curve in the domain 
    if(use1)
      v[[2]]=as.matrix(v[[2]][v_dom,])
    return(v)
  }
  # The function performs a weighted averaging of the elements in Y_inters_k to compute the new motif
  .compute_v_new <- function(Y_inters_k,Y_inters_supp,v_dom,v_len,p_k,d,m){
    # Y_inters_k: modified subsets of curves 
    v_new=matrix(NA,nrow=v_len,ncol=d)
    if(length(Y_inters_k)==1){ #there is only one modified subset of curves 
      v_new[v_dom,]=Y_inters_k[[1]]
      return(v_new)
    }
    Y_inters_supp=Reduce(rbind,Y_inters_supp) #concatena le matrici che formano Y_inters_supp in una unica matrice Y_inters_supp
    coeff_k=p_k[p_k>0]^m/(rowSums(Y_inters_supp)) # NB: divide for the length of the interval, not for the squared length!
    coeff_x=colSums(Y_inters_supp*coeff_k)#represents the weighted sum across the columns of Y_inters_supp
    coeff_x[colSums(Y_inters_supp)==0]=NA #prevents division by zero in the subsequent calculation.
    v_new[v_dom,]=Reduce('+',mapply('*',Y_inters_k,coeff_k,SIMPLIFY=FALSE))/coeff_x # For each element of Y_inters_k i multiply element-wise coefficienti_k and then i sum the result matrices
    return(v_new)
  }
  
  if(sum(p_k)==0){
    stop('Motif with no members! Degenerate cluster!')
  }
  v_len=length(v_dom)
  d=unlist(lapply(Y[[1]],ncol))[1] # dimension of curves
  # It applies a function using mapply to the elements of the Y list that have a membership 
  # probability greater than zero (Y[p_k>0]).
  Y_inters_k=mapply(function(y,s_k_i,d,use0,use1){# It processes each curve in the list, y, that has a nonzero membership probability to the cluster, k, and a corresponding shift value, s_k_i. It generates a modified subset of the curve based on the shift value
                      y_len=unlist(lapply(y,nrow))[1]
                      y_inters_k=list(y0=NULL,y1=NULL)
                      index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i)) # It calculates the indices used to extract the relevant part of the curve based on the shift value
                      y_inters_k$y0=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                                            matrix(y$y0[index[index<=y_len],],ncol=d),
                                            matrix(NA,nrow=sum(index>y_len),ncol=d)) # shfted portion of the curve
                      if(use1)
                        y_inters_k$y1=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                                            matrix(y$y1[index[index<=y_len],],ncol=d),
                                            matrix(NA,nrow=sum(index>y_len),ncol=d))
                      return(.select_domain(y_inters_k,v_dom,use0,use1))#remove any rows that are not part of the support, as indicated by the boolean vector v_dom
                    },Y[p_k>0],s_k[p_k>0],MoreArgs=list(d,use0,use1),SIMPLIFY=FALSE)  
  Y_inters_supp=lapply(Y_inters_k,.domain,use0) #For each shifted curve returns a vector of ones or zero, 1 iff the row does not contain NA(which elements of the curves belong to the support)
  Y_inters_k=mapply(function(y_inters_k,y_inters_supp,use0,use1){#setting the non-support elements to zero.
                      if(use0)
                        y_inters_k$y0[!y_inters_supp]=0
                      if(use1)
                        y_inters_k$y1[!y_inters_supp]=0
                      return(y_inters_k)},
                    Y_inters_k,Y_inters_supp,SIMPLIFY=FALSE,MoreArgs=list(use0,use1))
  v_new=list(v0=NULL,v1=NULL)
  if(use0)
    v_new$v0=.compute_v_new(lapply(Y_inters_k,function(y_inters_k) y_inters_k$y0),
                            Y_inters_supp,v_dom,v_len,p_k,d,m)
  if(use1)
    v_new$v1=.compute_v_new(lapply(Y_inters_k,function(y_inters_k) y_inters_k$y1),
                            Y_inters_supp,v_dom,v_len,p_k,d,m)
  # check if there are NA at the corners (it can happen after cleaning), and remove it
  range_v_new=range(which(.domain(v_new,use0)))
  v_dom_new=rep(FALSE,v_len)
  v_dom_new[(range_v_new[1]):(range_v_new[2])]=TRUE
  v_new=.select_domain(v_new,v_dom_new,use0,use1)
  if(range_v_new[1]>1){
    return(c(v_new,list(shift=range_v_new[1]-1)))
  }
  return(v_new)
}


.compute_Jk <- function(v,s_k,p_k,Y,alpha,w,m,c_k=NULL,keep_k=NULL,use0,use1){
  # Compute the objective function J for the motif k.
  # v: list of two elements, v0=v(x), v1=v'(x), matrices with d columns.
  # s_k: shift vector for motif k.
  # p_k: membership vector for motif k.
  # Y: list of N lists of two elements, Y0=y_i(x), Y1=y'_i(x), matrices with d columns, for d-dimensional curves.
  # alpha: weight coefficient between d0_L2 and d1_L2.
  # m>1: weighting exponent in least-squares functional.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  # keep_k: check c_k only when keep=TRUE for y_shifted.
  
  .diss_d0_d1_L2 <- function(y,v,w,alpha){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
    # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
    # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
    
    .diss_L2 <- function(y,v,w){
      # Dissimilarity index for multidimensional curves (dimension=d).
      # L2 distance with normalization on common support.
      # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
      # w: weights for the dissimilarity index in the different dimensions (w>0).
      
      sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
    }
    
    if(alpha==0){
      return(.diss_L2(y[[1]],v[[1]],w))
    }else if(alpha==1){
      return(.diss_L2(y[[2]],v[[2]],w))
    }else{
      return((1-alpha)*.diss_L2(y[[1]],v[[1]],w)+alpha*.diss_L2(y[[2]],v[[2]],w))
    }
  }
  .domain <- function(v,use0){
    if(use0){
      rowSums(!is.na(v[[1]]))!=0
    }else{
      rowSums(!is.na(v[[2]]))!=0
    }
  }
  .select_domain <- function(v,v_dom,use0,use1){
    if(use0)
      v[[1]]=as.matrix(v[[1]][v_dom,])
    if(use1)
      v[[2]]=as.matrix(v[[2]][v_dom,])
    return(v)
  }
  
  v_dom=.domain(v,use0)
  v_len=length(v_dom)
  v=.select_domain(v,v_dom,use0,use1)
  d=unlist(lapply(Y[[1]],ncol))[1]
  Y_inters_k=mapply(function(y,s_k_i,d){ # generate the shifted subsets of curves Y_inters_k. For each curve y in Y, it computes the shifted version y_inters_k according to its shift s_k_i, creating a list of such shifted curves.
                      y_len=unlist(lapply(y,nrow))[1]
                      y_inters_k=list(y0=NULL,y1=NULL)
                      index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
                      if(use0)
                        y_inters_k$y0=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                                            matrix(y$y0[index[index<=y_len],],ncol=d),
                                            matrix(NA,nrow=sum(index>y_len),ncol=d))
                      if(use1)
                        y_inters_k$y1=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                                            matrix(y$y1[index[index<=y_len],],ncol=d),
                                            matrix(NA,nrow=sum(index>y_len),ncol=d))
                      return(.select_domain(y_inters_k,v_dom,use0,use1))
                    },Y,s_k,MoreArgs=list(d),SIMPLIFY=FALSE) 
  if(!is.null(keep_k)){# If any of the supports have a length less than c_k, the function returns NA
    supp_inters_length=unlist(lapply(Y_inters_k[keep_k],function(y_inters_k) sum(.domain(y_inters_k,use0))))
    if(TRUE %in% (supp_inters_length<c_k))
      return(NA)
  }
  dist=unlist(mapply(.diss_d0_d1_L2,Y_inters_k,MoreArgs=list(v,w,alpha)))
  return(sum(dist*(p_k^m),na.rm=TRUE))
}


probKMA <- function(Y0,Y1=NULL,standardize=FALSE,K,c,c_max=Inf,P0=NULL,S0=NULL,
                       diss='d0_L2',alpha=NULL,w=1,m=2,
                       iter_max=1000,stop_criterion='max',quantile=NULL,tol=1e-8,
                       iter4elong=10,tol4elong=1e-3,max_elong=0.5,trials_elong=10,deltaJk_elong=0.05,max_gap=0.2,
                       iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                       return_options=TRUE,return_init=TRUE,worker_number=NULL){
  # Probabilistic k-mean with local alignment to find candidate motifs.
  # Y0: list of N vectors, for univariate curves y_i(x), or
  #     list of N matrices with d columns, for d-dimensional curves y_i(x),
  #     with the evaluation of curves (all curves should be evaluated on a uniform grid).
  #     When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
  # Y1: list of N vectors, for univariate derivative curves y'_i(x), or
  #     list of N matrices with d columns, for d-dimensional derivatibe curves y'_i(x),
  #     with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
  #     When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions.
  #     Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
  # standardize: if TRUE, each dimension is standardized (Z score on all the regions together).
  # K: number of motifs.
  # c: minimum motif lengths. Can be an integer (or a vector of K integers).
  # c_max: maximum motif lengths. Can be an integer (or a vector of K integers).
  # P0: initial membership matrix, with N row and K column (if NULL, a random P0 is choosen).
  # S0: initial shift warping matrix, with N row and K column (if NULL, a random S0 is choosen).
  # diss: dissimilarity. Possible choices are 'd0_L2', 'd1_L2', 'd0_d1_L2'. 
  # alpha: when diss='d0_d1_L2', weight coefficient between d0_L2 and d1_L2 (alpha=0 means d0_L2, alpha=1 means d1_L2).
  # w: vector of weights for the dissimilarity index in the different dimensions (w>0).
  # m>1: weighting exponent in least-squares functional.
  # iter_max: the maximum number of iterations allowed.
  # stop_criterion: criterion to stop iterate, based on the Bhattacharyya distance between memberships in subsequent iterations. 
  #                 Possible choices are: 'max' for the maximum of distances in the different motifs;
  #                                       'mean' for the average of distances in the different motifs;
  #                                       'quantile' for the quantile of distances in the different motifs (in this case, quantile must be provided).
  # quantile: probability in [0,1] to be used if stop.criterion='quantile'.
  # tol: method tolerance (method stops if the stop criterion <tol).
  # iter4elong: motifs elongation is performed every iter4elong iterations (if iter4elong>iter.max, no elongation is done).
  # tol4elong: tolerance on the Bhattacharyya distance (with the choice made in stop.criterion) for performing motifs elongation.
  # max_elong: maximum elongation allowed in a single iteration, as percentage of the motif length.
  # trials_elong: number of (equispaced) elongation trials at each side of the motif in a single iteration.
  # deltaJk_elong: maximum relative objective function increasing allowed in each motif elongation (for gaps and each side).
  # max_gap: maximum gap allowed in each alignment (percentage of motif length).
  # iter4clean: motif cleaning is performed every iter4clean iterations (if iter4clean>iter_max, no cleaning is done).
  # tol4clean: tolerance on the Bhattacharyya distance (with the choice made in stop_criterion) for performing motif cleaning.
  # quantile4clean: dissimilarity quantile to be used in motif cleaning.
  # return_options: if TRUE, the options K,c,diss,w,m are returned by the function.
  # return_init: if TRUE, P0 and S0 are returned by the function.
  # worker_number: number of CPU cores to be used for parallelization (default number of CPU cores -1). If worker_number=1, the function is run sequentially. 
  
  ### set parallel jobs #############################################################################
  start=proc.time()
  core_number <- detectCores()
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
    cl_probKMA=makeCluster(worker_number,timeout=60*60*24*30)# creates a cluster of worker nodes to enable parallel processing. It end up at most in 30 days
    clusterExport(cl_probKMA,c('.diss_d0_d1_L2','.domain','.select_domain'))# export functions in each node
    on.exit(stopCluster(cl_probKMA))
  }else{
    cl_probKMA=NULL
  }
  end=proc.time()
  #message('set parallel jobs: ',round((end-start)[3],2))
  
  ### check input ####################################################################################
  start=proc.time()
  # check dissimilarity
  if(!(diss %in% c('d0_L2','d1_L2','d0_d1_L2')))
    stop('Dissimilarity not valid. Possible choices are: \'d0_L2\', \'d1_L2\' and \'d0_d1_L2\'.')
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
  d=unique(unlist(lapply(Y0,ncol),use.names=FALSE)) # dimension of curves
  if(length(d)>1)
    stop('Curves have not all the same dimension.')
  Y0_NA=lapply(Y0,function(y,d) sum(rowSums(!is.na(y)) %in% c(0,d))==nrow(y),d)
  if(FALSE %in% Y0_NA){
    warning('y_j(x)=NA only in some dimensions, for some x. Putting NA in all dimensions.')
    Y0=lapply(Y0,
              function(y,d){
                y[rowSums(!is.na(y))!=d,]=NA
                return(y)},
              d) # Set NA for dimensions that are not completely on the support 
  }
  rm(Y0_NA)
  if(standardize){
    Y0_tot=Reduce(rbind,Y0)# Stack all the matrices in Y0 row-wise
    Y0_mean=colMeans(Y0_tot,na.rm=TRUE)
    Y0_sd=apply(Y0_tot,2,sd,na.rm=TRUE)# applied sd to each column(margin=2)
    rm(Y0_tot)
    Y0=lapply(Y0,function(y0,m,s) t((t(y0)-m)/s*100),m=Y0_mean,s=Y0_sd) # values in percentage
  }
  if(diss=='d0_d1_L2'){
    # check required input
    if(!is.numeric(alpha)){
      warning('alpha not valid. Setting alpha=0.5.')
      alpha=0.5
    }
    if((alpha<0)|(alpha>1)){
      warning('alpha not valid. Setting alpha=0.5.')
      alpha=0.5
    } else if(alpha==0){
      diss='d0_L2'
    } else if(alpha==1){
      diss='d1_L2'
    }
  }
  if(diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(Y0,function(y0) list(y0=y0,y1=NULL))
  }else{
    # check Y1
    if(is.null(Y1))
      stop(paste0('Y1 must be specified with diss=\'',diss,'\'.'))
    if(!is.list(Y1))
      stop('Y1 should be a list of vectors or matrices.')
    if((FALSE %in% lapply(Y1,is.matrix))&&(FALSE %in% lapply(Y1,is.vector)))
      stop('Y1 should be a list of vectors or matrices.')
    if(N!=length(Y1))
      stop('The number of curves in Y0 and derivatives Y1 should be the same.')
    Y1=lapply(Y1,as.matrix)
    d1=unique(unlist(lapply(Y1,ncol),use.names=FALSE)) # dimension of derivatives
    if(length(d1)>1)
      stop('Derivatives have not all the same dimension.')
    if(d!=d1)
      stop('The dimension of curves in Y0 and derivatives Y1 should be the same')
    Y1_NA=lapply(Y1,function(y,d) sum(rowSums(!is.na(y)) %in% c(0,d))==nrow(y),d)
    if(FALSE %in% Y1_NA){
      warning('y\'_j(x)=NA only in some dimensions, for some x. Putting NA in all dimensions.')
      Y1=lapply(Y1,
                function(y,d){
                  y[rowSums(!is.na(y))!=d,]=NA
                  return(y)},
                d)
    }
    rm(Y1_NA)
    Y0_NA=lapply(Y0,function(y) rowSums(is.na(y))!=0)
    Y1_NA=lapply(Y1,function(y) rowSums(is.na(y))!=0)
    diff_NA=mapply(function(y0,y1) y0!=y1,Y0_NA,Y1_NA)
    index_diff_NA=which(unlist(lapply(diff_NA,sum))!=0) # total number of differences between Y0 and Y1
    if(length(index_diff_NA)>0){
      warning('y_j(x) and y\'_j(x) are not both defined, for some x. Putting NA in that case.')
      same_NA=mapply(function(y0,y1,diff_NA){
                       y0[diff_NA]=NA
                       y1[diff_NA]=NA
                       return(list(y0,y1))
                     },Y0[index_diff_NA],Y1[index_diff_NA],diff_NA[index_diff_NA])
      Y0[index_diff_NA]=same_NA[1,]
      Y1[index_diff_NA]=same_NA[2,]
    }
    rm(Y0_NA,Y1_NA,diff_NA)
    if(standardize){
      Y1=lapply(Y1,function(y1,s) t(t(y1)/s*100),s=Y0_sd)
    }
    if(diss=='d1_L2'){
      alpha=1
      use0=FALSE
      use1=TRUE
      Y=lapply(Y1,function(y1) list(y0=NULL,y1=y1))
    }
    if(diss=='d0_d1_L2'){
      use0=TRUE
      use1=TRUE
      Y=mapply(function(y0,y1) list(y0=y0,y1=y1),Y0,Y1,SIMPLIFY=FALSE)
    }
  }
  # check required input
  if(missing(K))
    stop('K must be specified')
  if(missing(c))
    stop('c must be specified')
  # check K
  if(length(K)!=1)
    stop('Number of motifs K not valid.')
  if(K%%1!=0)
    stop('Number of motifs K should be an integer.')
  if(K<1)
    stop('Number of motifs K should be at least 1.')
  # check c
  if(!(length(c) %in% c(1,K)))
    stop('Minimum motif length c not valid.')
  if(sum(c%%1!=0))
    stop('Minimum motif lengths should be integers.')
  if(sum(c<=1))
    stop('Minimum motif lengths should be at least 2.')
  c=rep(c,length.out=K)
  # find all intervals contained in supp(y_i) and their lengths
  Y_intervals=lapply(Y0,
                     function(y,d){
                       y_not_NA=(rowSums(!is.na(y))==d) # It checks if the number of non-NA values in each row is equal to d.
                       intervals=which((y_not_NA[2:length(y_not_NA)]-y_not_NA[1:(length(y_not_NA)-1)])==1)+1 #It compute the intervals by finding the indices where the value changes from FALSE to TRUE and adds 1 to convert it into the right index.
                       if(y_not_NA[1])
                         intervals=c(1,intervals)# If the first element of y_not_NA is TRUE, it means the first row of the data y contains non-NA values. In this case, we add 1 to the beginning of the intervals vector to account for the first interval
                       intervals_lengths=rle(y_not_NA)$lengths[rle(y_not_NA)$values==TRUE]# lengths of the runs(interval) with TRUE values(consecutive) 
                       intervals_all=data.frame(start=intervals,
                                                end=intervals+intervals_lengths-1,
                                                length=intervals_lengths)
                       return(intervals_all)
                     },d) #These intervals represent consecutive regions in each element of Y0 have exactly d non-NA values.
  # check minimum motif length compatibility
  min_length=unlist(lapply(Y_intervals,
                           function(y_intervals,c){
                             return(sum(y_intervals$length>=c))},
                           max(c)))#this part calculates the number of intervals (motifs) whose length is greater than or equal to the specified motif length c
  if(0 %in% min_length)
    stop('Minimum motifs length is incompatible with supplied curves. Choose a smaller minimum motifs length.')
  # check c_max
  if(!(length(c_max) %in% c(1,K)))
    stop('Maximum motif length c_max not valid.')
  c_max=rep(c_max,length.out=K)
  for(k in 1:K){
    if(c_max[k]!=Inf){
      if(c_max[k]%%1!=0)
        stop('Maximum motif lengths should be integers.')
      if(c_max[k]<=1)
        stop('Maximum motif length should be at least 2.')
      # check maximum motif length compatibility
      max_length=unlist(lapply(Y_intervals,
                               function(y_intervals,c_max){
                                 return(sum(floor((1+max_gap)*y_intervals$length)>=c_max))},
                               c_max[k]))
      if(0 %in% max_length){
        warning('Maximum motif length is incompatible with supplied curves. Choosing default maximum motif length for motif ',k,'.')
        c_max[k]=Inf
      }
    }
    if(c_max[k]==Inf){
      c_max[k]=floor((1+max_gap)*min(unlist(lapply(Y_intervals,function(y_intervals) max(y_intervals$length)))))
    }
  }
  # check P0
  if(!is.null(P0)){
    if(sum(dim(P0)!=c(N,K))!=0){
      warning('Membership matrix dimensions not valid. Choosing random initial membership matrix.')
      P0=NULL
    }else{
      if(sum((P0<0)|(P0>1))){
        warning('Memberships should be non-negative and <=1. Choosing random initial membership matrix.')
        P0=NULL
        }else{
          if(sum(rowSums(P0)!=1)){
            warning('Memberships of each curve should sum to 1. Choosing random initial membership matrix.')
            P0=NULL
            }else{
              if(sum(colSums(P0)==0)){
                warning('Sum of memberships of each cluster should be positive. Choosing random initial membership matrix.')
                P0=NULL
              }
            }
        }
    }
  }
  # check S0
  if(!is.null(S0)){
    if(sum(dim(S0)!=c(N,K))!=0){
      warning('Shift warping matrix dimensions not valid. Choosing random initial shift warping matrix.')
      S0=NULL
    }else{
      Y_segments=mapply(function(y,S_i,c){# The outer lapply returns S_i for each i from 1...N. Then ot checks if there exists any NA values starting from the shifted curve y[s]
                          return(mapply(function(s,c) NA %in% y[s+seq_len(c)-1],S_i,c))},
                        Y0,lapply(seq_len(N),function(i) S0[i,]),MoreArgs=list(c),SIMPLIFY=TRUE)
      if(sum(Y_segments)){
        warning('Shift warping matrix not valid. Choosing random initial shift warping matrix.')
        S0=NULL
      }
    }#check whether there are any NA values in the segments of the curves corresponding to the centroid shifts in the S0 matrix.
  }
  # check weigths
  if(!is.vector(w)|!is.numeric(w))
    stop('Weights w not valid.')
  if(!(length(w) %in% c(1,d)))
    stop('Weights w not valid.')
  if(TRUE %in% (w<=0))
    stop('Weights w not valid.')
  w=rep(w,length.out=d)
  # check weighting exponent
  if(!is.numeric(m)|(length(m)>1))
    stop('Weighting exponent m not valid.')
  if(m<=1)
    stop('Weighting exponent m should be >1')
  # check maximum number of iterations
  if(!is.numeric(iter_max)|(length(iter_max)!=1))
    stop('Maximum number of iterations iter_max not valid.')
  if(((iter_max%%1)!=0)|(iter_max<=0))
    stop('Maximum number of iterations iter_max not valid.')
  # check stop criterion
  if(!(stop_criterion %in% c('max','mean','quantile')))
    stop('Stop criterion not valid. Possible choices are: \'max\', \'mean\', \'quantile\'.')
  if(stop_criterion=='quantile')
    if((!is.numeric(quantile))||(length(quantile)!=1)||(quantile<0)||(quantile>1))
      stop('quantile should be a number in [0,1].')
  # check tolerance
  if((!is.numeric(tol))||(length(tol)!=1)||(tol<=0))
    stop('Tolerance should be a positive number.')
  # check elongation input
  if((!is.numeric(iter4elong))||(length(iter4elong)!=1))
    stop('iter4elong not valid.')
  if(((iter4elong%%1)!=0)|(iter4elong<=0))
    stop('iter4elong not valid.')
  if((!is.numeric(tol4elong))||(length(tol4elong)!=1)||(tol4elong<=0))
    stop('tol4elong should be a positive number.')
  if((!is.numeric(max_elong))||(length(max_elong)!=1)||(max_elong<=0))
    stop('max_elong should be a positive number.')
  if((!is.numeric(trials_elong))||(length(trials_elong)!=1))
    stop('trials_elong not valid.')
  if(((trials_elong%%1)!=0)|(trials_elong<=0))
    stop('trials_elong not valid.')
  if((!is.numeric(deltaJk_elong))||(length(deltaJk_elong)!=1)||(deltaJk_elong<=0))
    stop('deltaJk_elong should be a positive number.')
  if((!is.numeric(max_gap))||(length(max_gap)!=1)||(max_gap<0))
    stop('max_gap should be a non-negative number.')
  # check cleaning input
  if((!is.numeric(iter4clean))||(length(iter4clean)!=1))
    stop('iter4clean not valid.')
  if(((iter4clean%%1)!=0)|(iter4clean<=0))
    stop('iter4clean not valid.')
  if((!is.numeric(tol4clean))||(length(tol4clean)!=1)||(tol4clean<=0))
    stop('tol4clean should be a positive number.')
  if((!is.numeric(quantile4clean))||(length(quantile4clean)!=1)||(N*K*quantile4clean<K)||(quantile4clean>1))
    stop('quantile4clean not valid.')
  end=proc.time()
  #message('check input: ',round((end-start)[3],2))
  
  ### initialize #############################################################################################
  start=proc.time()
  if(is.null(P0)){
    # create random membership matrix, with N rows and K columns
    P0=matrix(c(runif(N*(K-1)),rep(1,N)),nrow=N,ncol=K)#The first K-1 columns of P0 are filled with random values between 0 and 1.The last column of P0 is set to 1 for all rows --> each curve is assigned to at least one cluster or motif.
  
    P0=as.matrix(apply(P0,1,sort))# For each row i choose first the lower element and i compute the first Row of P0. then i go on with the second lower value for each row --> P0 becomes a KxN matrix( the last row is full of ones) --> a sort of comulative function for each column.
    if(K>1)
      P0=cbind(P0[1,],Reduce('cbind',lapply(2:K,function(k) P0[k,]-P0[k-1,]))) #I compute "atom" probability through differences in consecutiv4 rows and stacking it column-wise.Each row sums one and each column sums > 0
  }
  colnames(P0)=NULL
  P=P0
  if(is.null(S0)){
    # create shift warping matrix, with N rows and k columns
    S0=matrix(unlist(lapply(Y_intervals, # list of data-frame 
                            function(y_intervals,c,K){
                              s0=rep(NA,K)
                              for(k in 1:K){
                                y_intervals_k=y_intervals[y_intervals$length>=c[k],]# filtering intevals greater or equal to c[k]
                                y_starts=unlist(apply(y_intervals_k,1,
                                                      function(y_interval)
                                                        return(y_interval[1]:(y_interval[2]-c[k]+1))),# return start_index:(end-c[k]+1)
                                                use.names=FALSE)
                                s0[k]=sample(y_starts,1) #randomly selects one of these valid starting positions for each motif.
                              }
                              return(s0)
                            },c,K)),
              nrow=N,ncol=K,byrow=TRUE)
  }
  S=S0

  # create empty motifs
  V=lapply(c,
           function(c_k,d){
             v=list(v0=NULL,v1=NULL)
             if(use0)
               v$v0=matrix(0,nrow=c_k,ncol=d)
             if(use1)
               v$v1=matrix(0,nrow=c_k,ncol=d)
             return(v)
           },d)
  end=proc.time()
  #message('initialize: ',round((end-start)[3],2))
  
  ### iterate #############################################################################################
  iter=0
  J_iter=c()
  BC_dist_iter=c()
  BC_dist=+Inf
  while((iter<iter_max)&(BC_dist>tol)){
    iter=iter+1
    #message('iter ',iter)
    
    ##### clean motifs ###################################################################################
    start=proc.time()
    P_old=P
    if((iter>1)&&(!(iter%%iter4clean))&&(BC_dist<tol4clean)){
      keep=D<quantile(D,quantile4clean) #logical matrix.If a motif's dissimilarity values are below the specified quantile, it is considered significant or relevant and retained(N x K matrix)
      empty_k=which(colSums(keep)==0) #these are motifs that are not considered significant or relevant based on the dissimilarity values
      if(length(empty_k)>0){ # If there exist such a motif completely not significant 
        for(k in empty_k)
          keep[which.min(D[,k]),k]=TRUE #If a motif is not enterly significant for any curves i set TRUE for the curves with minimum dissimilarity ensuring that it's not entirely discarded.
      }
      P[keep]=1
      P[!keep]=0
    }
    end=proc.time()
    #message('  clean: ',round((end-start)[3],2))
    
    ##### compute motifs ###################################################################################
    start=proc.time()
    S_k=split(S,rep(seq_len(K),each=N))# S is split into a list of matrices S_k, where each element of the list corresponds to one cluster K (group) of curves
    P_k=split(P,rep(seq_len(K),each=N))
    V_dom=lapply(V,.domain,use0)
    V_new=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE) # list with computed motif for each cluster of curves 
    changed_s=which(unlist(lapply(V_new,length))>2) #This line identifies the clusters (groups) of curves where the computed motifs have more than two elements
    for(k in changed_s){
      S[,k]=S[,k]+V_new[[k]]$shift # Update S for the k-th cluster by adding the shift values computed for the motif in the V_new list.
      V_new[[k]]=V_new[[k]][c('v0','v1')]#This step updates the V_new list for the k-th cluster by keeping only the relevant information for the computed motif (it retains the elements named 'v0' and 'v1')
    }
    S_k=split(S,rep(seq_len(K),each=N))
    V_dom=lapply(V_new,.domain,use0)
    end=proc.time()
    #message('  compute motifs: ',round((end-start)[3],2))
    
    ##### elongate motifs #################################################################################
    start=proc.time()
    if((iter>1)&&(!(iter%%iter4elong))&&(BC_dist<tol4elong)){
      # fill
      with_gaps=which(unlist(lapply(V_dom,function(v_dom) sum(!v_dom)!=0)))#identifies the indices of clusters of curves that are out of the domain
      if(length(with_gaps)>0){
        V_dom_filled=lapply(V_dom[with_gaps],function(v_dom) rep_len(TRUE,length(v_dom)))# For clusters with gaps, this line creates a new V_dom_filled list where the gaps are filled with TRUE
        V_filled=mapply(.compute_motif,V_dom_filled,S_k[with_gaps],P_k[with_gaps],MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)#calculate the motif based on the filled data
        Jk_before=mapply(.compute_Jk,
                         V_new[with_gaps],S_k[with_gaps],P_k[with_gaps],
                         MoreArgs=list(Y=Y,alpha=alpha,w=w,m=m,use0=use0,use1=use1))# compute the metric J for cluster k with gaps
        Jk_after=mapply(.compute_Jk,
                        V_filled,S_k[with_gaps],P_k[with_gaps],
                        MoreArgs=list(Y=Y,alpha=alpha,w=w,m=m,use0=use0,use1=use1))#The Jk values after elongation (using filled data) are computed for the clusters with gaps
        fill=(Jk_after-Jk_before)/Jk_before<deltaJk_elong #The fill logical vector is created to identify clusters where the relative change in Jk values after elongation is less than a specified threshold deltaJk_elong
        V_dom[with_gaps[fill]]=V_dom_filled[fill]
        V_new[with_gaps[fill]]=V_filled[fill]
      }
      # elongate
      len_dom=unlist(lapply(V_dom,length))# vector of lenghts (for each motif)
      len_max_elong=mapply(min,floor(len_dom*max_elong),(c_max-len_dom))# taking the minimum value between...
      len_elong=lapply(len_max_elong, # return a list where for each cluster, all the potential lengths for elongation are computed
                       function(len_max_elong){
                         if(len_max_elong<=trials_elong){
                           len_elong=seq_len(len_max_elong)#sequence of lengths ranging from 1 to len_max_elong
                         }else{
                           len_elong=round(seq(1,len_max_elong,length.out=trials_elong))#sequence of approximately trials_elong equally spaced values between 1 and len_max_elong
                         }
                         return(len_elong)
                       })
      # left and right elongation
      keep=D<quantile(D,0.25)# The variable keep is assigned a logical matrix where the elements are TRUE if the corresponding element in the dissimilarity matrix D is less than the 25th percentile of D.
      empty_k=which(colSums(keep)==0)#Indices of clusters that have all dissimilarity values above the 25th percentile are identified.
      if(length(empty_k)>0){
        for(k in empty_k)
          keep[which.min(D[,k]),k]=TRUE # ensure that all clusters have at least one motif that is considered significant 
      }
      res_left_right=mapply(function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
                              if(length(len_elong_k)==0){ # If there are no valid elongation lengths for the cluster, it returns the original motif and domain vectors
                                return(list(v_new=v_new_k,
                                            v_dom=v_dom_k,
                                            s_k=s_k))
                              }
                              #create a matrix s_k_elong_left_right that contains different shifted warping matrices obtained by applying different elongation lengths to the original s_k matrix.
                              s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k) s_k-len_elong_k),(length(len_elong_k)+1):1)[-1]#for each element in the list the provided function subtracts it from s_k, resulting in a shifted version of s_k(matrix).At the end we have a vector of lists of shifted warping matrices.(-1 removes the null case with no elongation)
                                                                   function(len_elong_k_left) # modify previous matrices.It addresses the case where one side (left or right) has more potential for elongation than the other side.
                                                                     lapply(c(0,len_elong_k[len_elong_k<=(max(len_elong_k)-len_elong_k_left)]),
                                                                            function(len_elong_k_right) 
                                                                              c(rep_len(TRUE,len_elong_k_left),v_dom_k,rep_len(TRUE,len_elong_k_right)))),
                                                            recursive=FALSE)[-1] # non ho ben capito 
                              v_elong_left_right=mapply(.compute_motif,v_dom_elong_left_right,s_k_elong_left_right,
                                                        MoreArgs=list(p_k,Y,m,use0,use1),SIMPLIFY=FALSE)# compute new motif(return a list)
                              start_with_NA=unlist(lapply(v_elong_left_right,length))>2#This line creates a logical vector indicating which elongated motifs have more than two elements
                              v_elong_left_right=v_elong_left_right[!start_with_NA]
                              s_k_elong_left_right=s_k_elong_left_right[!start_with_NA]
                              Jk_before=.compute_Jk(v_new_k,s_k,p_k,Y,alpha,w,m,use0=use0,use1=use1)
                              c_k_after=floor(unlist(lapply(lapply(v_elong_left_right,.domain,use0),length))*(1-max_gap))#This line calculates a new potential motif length c_k_after for each elongated motif after applying a floor function and reducing it by a factor of max_gap
                              c_k_after[c_k_after<c]=c
                              Jk_after=unlist(mapply(.compute_Jk,v_elong_left_right,s_k_elong_left_right,c_k_after,
                                                     MoreArgs=list(p_k=p_k,Y=Y,alpha=alpha,w=w,m=m,keep_k=keep_k,use0=use0,use1=use1)))#This line calculates the Jk metric for each elongated motif after applying the changes in c_k_after
                              best_elong=which.min((Jk_after-Jk_before)/Jk_before)#This line identifies the index of the motif with the smallest relative change in the Jk metric after elongation.
                              if(length(best_elong)>0){#If there is a potential elongation that improves the Jk metric
                                elongate=((Jk_after-Jk_before)/Jk_before)[best_elong]<deltaJk_elong
                              }else{
                                elongate=FALSE
                              }
                              if(elongate){#If elongate is TRUE, the code returns the elongated motif, domain vector, and shifted warping matrix
                                return(list(v_new=v_elong_left_right[[best_elong]],
                                            v_dom=v_dom_elong_left_right[[best_elong]],
                                            s_k=s_k_elong_left_right[[best_elong]]))
                              }else{
                                return(list(v_new=v_new_k,
                                            v_dom=v_dom_k,
                                            s_k=s_k))
                              }
                            },V_new,V_dom,S_k,P_k,len_elong,split(keep,rep(1:K,each=N)),c)
      V_new=res_left_right[1,]
      V_dom=res_left_right[2,]
      S_k=res_left_right[3,]
      rm(res_left_right)
      S=matrix(unlist(S_k),ncol=K)
    }
    end=proc.time()
    #message('  elongate: ',round((end-start)[3],2))
    
    ##### find shift warping minimizing dissimilarities ###################################################
    start=proc.time()
    c_k=floor(unlist(lapply(V_new,function(v_new) unlist(lapply(v_new,nrow))[1]))*(1-max_gap))#For each cluster, this line calculates the desired motif length (c_k) based on the first motif in the V_new list
    c_k[c_k<c]=c
    c_k=rep(c_k,each=length(Y))
    YV=expand.grid(Y,V_new)#This line creates a data frame YV that contains all combinations of the original curves Y and the motifs in V_new
    SD=.mapply_custom(cl_probKMA,.find_min_diss,YV[,1],YV[,2],c_k,
                      MoreArgs=list(alpha=alpha,w=w,d=d,use0=use0,use1=use1),SIMPLIFY=TRUE)#find the shift warpings that minimize dissimilarities.The function .find_min_diss is applied to each combination of curve Y and motif in V_new (contained in YV), aiming to determine the optimal shift warping that minimizes dissimilarities.
    S_new=matrix(SD[1,],ncol=K) #each column corresponds to a cluster (motif).The shift warpings and dissimilarities obtained from the optimization process are stored in these matrices.
    D_new=matrix(SD[2,],ncol=K)
    end=proc.time()
    #message('  find shift: ',round((end-start)[3],2))
    
    ##### compute memberships #############################################################################
    start=proc.time()
    P_new=matrix(0,nrow=N,ncol=K) #Initializes a membership matrix P_new with dimensions N (number of curves) by K (number of motifs).
    D0=apply(D_new,2,'%in%',0) #if dist(y_i,v_k)=0 for some k, set p(i,k)=1 and p(i,h)=0 for h!=k
    for(i in which(rowSums(D0)>1)){ #This loop handles the case where a curve has a dissimilarity of 0 from two or more motifs. It randomly selects one of the motifs with dissimilarity 0 and sets the membership to 1 for that motif while setting membership to 0 for others
      warning(paste0('Curve ',i,' has dissimilarity 0 from two different motifs. Using only one of them...'))
      D0_select=sample(which(D0[i,]),1)
      D0[i,-D0_select]=FALSE
    }
    D0_index=rowSums(D0)==1 # identifies curves with dissimilarities exactly equal to 0 from a single motif
    P_new[D0_index,]=1*D0[D0_index,]
    # if dist(y_i,v_k)>0 for all k
    Dm=as.matrix(D_new[!D0_index,]^(1/(m-1)))
    P_new[!D0_index,]=1/(Dm*rowSums(1/Dm))
    for(k in which(colSums(P_new)==0)){# check degenerate clusters (zero membership)
      warning(paste0('Motif ',k,' is degenerate (zero membership). Selecting a new center...'))
      P_new[which.min(D_new[,k]),k]=1 # set the membership of the curve with the minimum dissimilarity to that motif to 1
    }
    end=proc.time()
    #message('  compute memberships: ',round((end-start)[3],2))
    
    ##### evaluate objective function #####################################################################
    J_iter[iter]=sum(D_new*(P_new^m),na.rm=TRUE)

    ##### compute Bhattacharyya distance between P_old and P_new ##########################################
    BC_dist_k=-log(rowSums(sqrt(P_old*P_new)))
    if(stop_criterion=='max')
      BC_dist=max(BC_dist_k)
    if(stop_criterion=='mean')
      BC_dist=mean(BC_dist_k)
    if(stop_criterion=='quantile')
      BC_dist=quantile(BC_dist_k,prob,type=1)#type = 1: Linear interpolation method. The quantile value is computed as a linear interpolation between the two nearest data points
    
    BC_dist_iter[iter]=BC_dist    

    ##### update ##########################################################################################
    V=V_new
    P=P_new
    S=S_new
    D=D_new
  }

  ### prepare output ####################################################################################
  start=proc.time()
  if(iter==iter_max){
    warning('maximum number of iterations reached, method stops')
  }
  # compute motifs
  S_k=split(S,rep(seq_len(K),each=N)) #  shifts and memberships matrices for each cluster(Nx1 matrix)
  P_k=split(P,rep(seq_len(K),each=N))
  if(!use0){
    use0=TRUE
    Y=mapply(function(y,y0) list(y0=y0,y1=y$y1),Y,Y0,SIMPLIFY=FALSE)
  }
  V=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
  # compute cleaned motifs
  keep=D<quantile(D,quantile4clean) #matrix keep, which is TRUE where the dissimilarity values in D are below the specified quantile (quantile4clean).
  empty_k=which(colSums(keep)==0)
  if(length(empty_k)>0){
    for(k in empty_k)
      keep[which.min(D[,k]),k]=TRUE
  }
  P_clean=P
  P_clean[keep]=1
  P_clean[!keep]=0
  P_k=split(P_clean,rep(seq_len(K),each=N))
  S_clean=S
  V_clean=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
  changed_s=which(unlist(lapply(V_clean,length))>2) #Identify clusters (changed_s) where the computed cleaned motifs have more than two elements(well computed).
  for(k in changed_s){
    S_clean[,k]=S_clean[,k]+V_clean[[k]]$shift
    V_clean[[k]]=V_clean[[k]][c('v0','v1')]# retains only v0 and v1
  }
  S_k=split(S_clean,rep(seq_len(K),each=N))
  V_dom=lapply(V_clean,.domain,use0)
  # compute dissimilarities from cleaned motifs
  D_clean=mapply(function(s_k,v_dom,v_clean,Y){
                   Y_in_motifs=mapply(function(y,s){#For each combination of s_k, v_dom, and v_clean, the inner mapply function calculates Y_in_motifs by shifting the original data Y according to s_k and adjusting the domain based on v_dom.
                                        if(use0){
                                          d=ncol(y$y0) # dimension of curves
                                          y_len=nrow(y$y0)
                                          index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
                                          y$y0=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                                                     matrix(y$y0[index[index<=y_len],],ncol=d),
                                                     matrix(NA,nrow=sum(index>y_len),ncol=d))
                                        }
                                        if(use1){
                                          d=ncol(y$y1) # dimension of curves
                                          y_len=nrow(y$y1)
                                          index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
                                          y$y1=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                                                     matrix(y$y1[index[index<=y_len],],ncol=d),
                                                     matrix(NA,nrow=sum(index>y_len),ncol=d))
                                        }
                                        y=.select_domain(y,v_dom,use0,use1)
                                        return(y)
                                      },Y,s_k,SIMPLIFY=FALSE)
                   d=unlist(lapply(Y_in_motifs,.diss_d0_d1_L2,.select_domain(v_clean,v_dom,use0,use1),w,alpha))
                   return(d)
                 },S_k,V_dom,V_clean,MoreArgs=list(Y))
  output=list(Y0=Y0,Y1=Y1,
              V0=lapply(V,function(v) v$v0),V1=lapply(V,function(v) v$v1),
              V0_clean=lapply(V_clean,function(v) v$v0),V1_clean=lapply(V_clean,function(v) v$v1),
              P=P,P_clean=P_clean,S=S,S_clean=S_clean,
              D=D,D_clean=D_clean,iter=iter,J_iter=J_iter,BC_dist_iter=BC_dist_iter)
  if(return_options){
    output=c(output,list(standardize=standardize,K=K,c=c,c_max=c_max,diss=diss,alpha=alpha,w=w,m=m,
                         iter_max=iter_max,stop_criterion=stop_criterion,quantile=quantile,tol=tol,
                         iter4elong=iter4elong,tol4elong=tol4elong,max_elong=max_elong,trials_elong=trials_elong,deltaJk_elong=deltaJk_elong,max_gap=max_gap,
                         iter4clean=iter4clean,tol4clean=tol4clean))
  }
  if(return_init){
    output=c(output,list(P0=P0,S0=S0))
  }
  end=proc.time()
  #message('output: ',round((end-start)[3],2))
  

  ### return output ####################################################################################
  return(output)
}


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

probKMA_silhouette <- function(probKMA_results,align=FALSE,plot=TRUE){
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
  curves_in_motifs=apply(probKMA_results$P_clean,2,function(P_clean_k) which(P_clean_k==1))# for each column(motif P_clean_k) the function which(P_clean_k == 1) returns a list where for each column we have indices of curves assigned to that motif.
  if(!is.null(ncol(curves_in_motifs))){ #If there are curves assigned to the motifs
    curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs)))#the curves_in_motifs variable is split into a list where each element of the list corresponds to a motif (column) and contains the indices of curves assigned to that motif.
  }
  curves_in_motifs_number=colSums(probKMA_results$P_clean) #The curves_in_motifs_number vector stores the count of curves assigned to each motif.
  S_clean_k=mapply(function(s_k,curves_in_motif) s_k[curves_in_motif],
                            split(probKMA_results$S_clean,rep(seq_len(K),each=N)),curves_in_motifs,SIMPLIFY=FALSE)# This function takes the shift warping matrix s_k and extracts only the rows corresponding to curves assigned to the motif using the indices from curves_in_motif
                                                                                                                  # The result is stored in the S_clean_k list, where each element corresponds to a motif and contains the shift warping matrices for the curves assigned to that motif
                                                                                                                  # The number of rows in the shift warping matrix will be equal to the number of curves assigned to that motif, and the number of columns will be determined by the length of the shift warping vectors
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
                            },curves_in_motifs,S_clean_k,V_dom,SIMPLIFY=FALSE),recursive=FALSE)#The result is a matrix-like object Y_in_motifs that represents the curves after applying the shift warping transformations and domain adjustments
  Y_motifs=rep.int(1:K,curves_in_motifs_number)#Y_motifs associates each curve with its corresponding motif index.
  YY=combn(Y_in_motifs,2,simplify=FALSE) #It generates a list of matrices where each column represents a pair of motifs
  YY=array(unlist(YY,recursive=FALSE),dim=c(2,length(YY)))#: This line converts the list of matrices YY into a two-dimensional array. Each column of the array corresponds to a pair of motifs.
  YY_lengths=combn(V_length[Y_motifs],2) #This line computes all possible combinations of pairs of elements from the V_length vector, considering only the motifs associated with each curve.
                                         #The first row of the matrix contains the lengths of the motifs for the first motif in each pair, and the second row contains the lengths of the motifs for the second motif in each pair.
  swap=YY_lengths[1,]<YY_lengths[2,]
  YY[,swap]=rbind(YY[2,swap],YY[1,swap])#This line swaps the columns of the YY array based on the condition in the swap vector.
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
    # probKMA_result$c = vector: Cluster assignments for each curve
    # Y_motifs = vector: Motif assignments for each curve
    # probKMA_results$c[Y_motifs]: each curve is associated with the cluster to which its motif belongs
    cc_motifs=apply(combn(probKMA_results$c[Y_motifs],2),2,min)# This gives the minimum index among the pairs of motif indices, and it is used as an indicator of the overlap required for alignment
    SD=mapply(.find_min_diss,YY[1,],YY[2,],cc_motifs,
              MoreArgs=list(alpha=alpha,w=w,d,use0,use1),SIMPLIFY=TRUE)
  }
  YY_D=matrix(0,nrow=length(Y_motifs),ncol=length(Y_motifs))
  YY_D[lower.tri(YY_D)]=SD[2,] #this step effectively fills in the lower triangular part of the matrix with the distances between the pairs of curve pieces.
  YY_D=YY_D+t(YY_D) #Here,the matrix is made symmetric by adding its transpose to itself. This step ensures that the matrix YY_D is symmetric.
  
  # compute intra-cluster distances
  intra=Reduce(cbind,lapply(Y_motifs,function(motif) Y_motifs==motif))# This matrix computes the intra-cluster distances. It is a binary matrix, where each row corresponds to a specific motif and each column corresponds to a different curve piece
  diag(intra)=FALSE
  a=colSums(intra*YY_D)/(curves_in_motifs_number[Y_motifs]-1)#This vector calculates the average intra-cluster distance for each motif
  
  # compute inter-cluster distances
  b_k=Reduce(rbind,lapply(1:(K-1),
                          function(k){
                            inter=Reduce(cbind,lapply(Y_motifs,
                                                      function(motif) Y_motifs==ifelse(motif+k>K,(motif+k)%%K,motif+k)))
                            b_k=colSums(inter*YY_D)/(curves_in_motifs_number[ifelse(Y_motifs+1>K,(Y_motifs+1)%%K,Y_motifs+1)])
                            return(b_k)})) # b_k is a matrix where each row corresponds to a different motif, and each column corresponds to a different curve piece
  if(is.matrix(b_k)){
    b=apply(b_k,2,min)
  }else{
    b=b_k
  } #b is a vector that calculates the minimum value along each column of b_k, resulting in the minimum inter-cluster distance for each curve piece.
  
  # compute silhouette
  silhouette=(b-a)/pmax(a,b)
  silhouette[is.nan(silhouette)]=0
  
  # compute average silhouette per cluster
  silhouette_average=rep(NA,K) #store the average silhouette value for each cluster
  for(k in seq_len(K)){
    silhouette_k=silhouette[Y_motifs==k] # Extract the silhouette values for the curve pieces that belong to cluster k
    curves_in_motifs[[k]]=curves_in_motifs[[k]][order(silhouette_k,decreasing=TRUE)] # Reorder the curve pieces within the cluster based on their silhouette values(FOR PLOT)
    silhouette[Y_motifs==k]=sort(silhouette_k,decreasing=TRUE)#Sort the silhouette values in descending order for the curve pieces in cluster k(FOR PLOT)
    silhouette_average[k]=mean(silhouette_k) # compute the mean  
  }

  ### plot silhouette ########################################################################################
  if(plot){
    n=length(silhouette)
    sil=rev(silhouette)
    y=barplot(sil,space=c(0,rev(diff(Y_motifs))),xlab='Silhouette index',names='',
              xlim=c(min(0,min(sil)-0.05),1.2),horiz=TRUE,las=1,mgp=c(2.5,1,0),col='gray')
    text(ifelse(sil>=0,sil+0.03,sil-0.03),ifelse(sil>0,y,y+0.2),labels=rev(unlist(curves_in_motifs)),cex=0.5)
    title(main='Silhouette plot')
    title(sub=paste("Average silhouette width:",round(mean(silhouette_average),2)),adj=0)
    mtext(paste("n =",n),adj=0)
    mtext(substitute(K~~"motifs",list(K=K)),adj=1)
    mtext(expression(paste(j," : ",n[j]," | avg ",s[i])),adj=1.04,line=-1.2)
    y=rev(y)
    for(k in seq_len(K)){
      y_k=mean(y[Y_motifs==k])
      text(1,y_k,
           paste0(k,": ",curves_in_motifs_number[k]," | ",format(silhouette_average[k],digits=1,nsmall=2)),xpd=NA,adj=0.1)
    }
  }
  
  return(list(silhouette=silhouette,motifs=Y_motifs,curves=unlist(curves_in_motifs),
              silhouette_average=silhouette_average))
}


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
  core_number <- detectCores()
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
      cl_find=makeCluster(worker_number,timeout=60*60*24*30)
      clusterExport(cl_find,c('name','names_var','Y0','Y1','probKMA_options',
                              'probKMA','probKMA_plot','probKMA_silhouette','.compute_motif', # delete these in the package
                              '.mapply_custom','.diss_d0_d1_L2','.domain','.select_domain','.find_min_diss','.compute_Jk'),envir=environment()) 
      clusterCall(cl_find,function()library(parallel,combinat))
      on.exit(stopCluster(cl_find))
    }else{
      cl_find=NULL
    }
  }
  
  ### run probKMA ##########################################################################################
  i_c_K=expand.grid(seq_len(n_init),c,K) #all possible combinations
  results=.mapply_custom(cl_find,function(K,c,i){
                                    dir.create(paste0(name,"_K",K,"_c",c),showWarnings=FALSE) #"name"_K"_c"
                                    files=list.files(paste0(name,"_K",K,"_c",c)) #list files in the directory that was just created 
                                    message("K",K,"_c",c,'_random',i) #i prints the message
                                    if(paste0('random',i,'.RData') %in% files){#checks if a specific RData file named 'random' followed by the value of i with extension .RData is present in the list of files in the directory. If it is present, it loads the data from that RData file.
                                      load(paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
                                      return(list(probKMA_results=probKMA_results,
                                                  time=time,silhouette=silhouette))
                                    }else{
                                      iter=iter_max=1
                                      while(iter==iter_max){
                                        start=proc.time()
                                        probKMA_results=do.call(probKMA,c(list(Y0=Y0,Y1=Y1,K=K,c=c),probKMA_options)) #execute ProbKMA with default arguments
                                        end=proc.time()
                                        time=end-start
                                        iter=probKMA_results$iter
                                        iter_max=probKMA_results$iter_max
                                        if(iter==iter_max)
                                          warning('Maximum number of iteration reached. Re-starting.')
                                      }
                                      pdf(paste0(name,"_K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10) #open pdf file
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
  results=split(results,list(factor(i_c_K[,2],c),factor(i_c_K[,3],K))) #results is being split into subgroups based on two factors: the values of c and the value of K.This operation creates subgroups of results based on unique combinations of c and K
  results=split(results,rep(K,each=length(c))) #each subgroup corresponds to a specific value of K.
  
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
    # calculate the motif length from a nested list structure of results
    motif_length=mapply(function(results){
                          motif_length=mapply(function(results){
                                                motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                                                           function(results){
                                                                                             unlist(lapply(results$probKMA_results$V0,nrow)) #The lapply function iterates over each element in the current level of the nested results.
                                                                                                                                             #For each element in the nested results, it accesses the probKMA_results$V0 attribute and calculates the number of rows in each matrix using nrow
                                                                                           })))# non ho capito il senso dei primi due mapply
                                                return(as.matrix(motif_length))
                                              },results,SIMPLIFY=FALSE)
                        },results,SIMPLIFY=FALSE)#The result should be a matrix where each row corresponds to a specific level of nesting in the results structure, and each column corresponds to a different motif's length.
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


.probKMA_silhouette_filter <- function(probKMA_results,silhouette,sil_threshold=0.5,size_threshold=2){
  # Filter the motifs found by probKMA on the basis of a threshold on the average silhouette index
  # probKMA_results: output of probKMA function (with return_options=TRUE)
  # silhouette: output of probKMA_silhouette function
  # sil_threshold: threshold on the average silhouette index
  # size_threshold: threshold on the size of the motif (number of curves in the cluster)
  
  index_sil=which(silhouette$silhouette_average>=sil_threshold)#indices of clusters (motifs) with silhouette average values greater than or equal to a specified sil_threshold.
  index_size=which(colSums(probKMA_results$P_clean)>=size_threshold) #index_size is a vector containing the indices of clusters (motifs) with number of curves in the cluster greater than or equal to a specified size_threshold.
  index=intersect(index_sil,index_size)
  if(length(index)==0)
    return(NULL)
  return(list(V0_clean=probKMA_results$V0_clean[index],V1_clean=probKMA_results$V1_clean[index],
              D=probKMA_results$D[,index],D_clean=probKMA_results$D_clean[,index],P=probKMA_results$P[,index],P_clean=probKMA_results$P_clean[,index],
              c=probKMA_results$c[index],K=rep(probKMA_results$K,length(index))))#If there are clusters that satisfy both criteria, the function returns a list containing specific elements of probKMA_results corresponding to the selected clusters
}


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
                                  })#This nested loop structure loads the filtered probKMA_results for each combination of K, c, and n_init, applies the filtering function .probKMA_silhouette_filter, and collects the resulting motifs.
                          })
                 })
  motifs=unlist(unlist(motifs,recursive=FALSE),recursive=FALSE)#This line flattens the nested list structure into a single list of motif results.
  if(is.null(unlist(motifs)))
    stop("No motif present after filtering. Please re-run the function with less stringent parameters.")
  V0_clean=unlist(lapply(motifs,function(motifs) motifs$V0_clean),recursive=FALSE)# Extract V0_clean from each motif
  V_clean_length=unlist(lapply(V0_clean,length)) # length of a matrix is the product of rows and columns 
  index=order(V_clean_length,decreasing=TRUE) # order from the longest to the shortest
  V0_clean=V0_clean[index] # reodered in descreasing order
  V1_clean=unlist(lapply(motifs,function(motifs) motifs$V1_clean),recursive=FALSE)[index] # reodered 
  D_clean=Reduce(cbind,lapply(motifs,function(motifs) motifs$D_clean))[,index]# horizontally concatenated D_clean
  P_clean=Reduce(cbind,lapply(motifs,function(motifs) motifs$P_clean))[,index]
  c=Reduce(c,lapply(motifs,function(motifs) motifs$c))[index]
  K=Reduce(c,lapply(motifs,function(motifs) motifs$K))[index]

  ### output ##################################################################################################
  load(paste0(name,"_K",K[1],"_c",c[1],'/random',1,'.RData'))
  return(list(V0_clean=V0_clean,V1_clean=V1_clean,D_clean=as.matrix(D_clean),P_clean=as.matrix(P_clean),c=c,K=K,
              Y0=probKMA_results$Y0,Y1=probKMA_results$Y1,
              diss=probKMA_results$diss,alpha=probKMA_results$alpha,w=probKMA_results$w,max_gap=probKMA_results$max_gap))
}


.find_occurrences <- function(v,Y,R,alpha,w,c_k,use0,use1){
  # Find occurrences of a motif in a set of curves (dimesion=d), with dissimilarity lower than R.
  # Return curve id, shift and dissimilarity.
  # v: list of 2 elements, v0, v1, matrices with d columns.
  # Y: list of N lists of two elements, Y0, Y1, matrices with d columns.
  # R: maximum dissimilarity allowed. 
  # alpha: if diss_fun=diss_d0_d1_L2, weight coefficient between d0_L2 and d1_L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  
  .diss_d0_d1_L2 <- function(y,v,w,alpha){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
    # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
    # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
    
    .diss_L2 <- function(y,v,w){
      # Dissimilarity index for multidimensional curves (dimension=d).
      # L2 distance with normalization on common support.
      # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
      # w: weights for the dissimilarity index in the different dimensions (w>0).
      
      sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
    }
    
    if(alpha==0){
      return(.diss_L2(y[[1]],v[[1]],w))
    }else if(alpha==1){
      return(.diss_L2(y[[2]],v[[2]],w))
    }else{
      return((1-alpha)*.diss_L2(y[[1]],v[[1]],w)+alpha*.diss_L2(y[[2]],v[[2]],w))
    }
  }
  # returns a logical vector where FALSE means that the corresponing rows is full of NA
  .domain <- function(v,use0){
    if(use0){
      rowSums(!is.na(v[[1]]))!=0
    }else{
      rowSums(!is.na(v[[2]]))!=0
    }
  }
  # Select only the domain where the curve is defined 
  .select_domain <- function(v,v_dom,use0,use1){
    if(use0)
      v[[1]]=as.matrix(v[[1]][v_dom,])
    if(use1)
      v[[2]]=as.matrix(v[[2]][v_dom,])
    return(v)
  }
  
  v_dom=.domain(v,use0)
  v_len=length(v_dom)
  SD_motif=lapply(Y,
                  function(y){
                    y_len=unlist(lapply(y,nrow))[1] # lenght of the domain(curve)
                    s_rep=seq_len(y_len-v_len+1) # possible starting points for the alignment between y and the motif v.
                    y_rep=lapply(s_rep,
                                 function(i){
                                   y_rep=list(y0=NULL,y1=NULL)
                                   if(use0)
                                     y_rep$y0=as.matrix(y$y0[i-1+seq_len(v_len),])
                                   if(use1)
                                     y_rep$y1=as.matrix(y$y1[i-1+seq_len(v_len),])
                                   y_rep=.select_domain(y_rep,v_dom,use0,use1)
                                 })#a sub-curve y_rep is extracted from the curve y, spanning the same length as the motif v
                    valid=unlist(lapply(lapply(y_rep,.domain,use0),sum))>=c_k #Valid sub-curves are selected based on the minimum length requirement c_k
                    s_rep=s_rep[valid]
                    y_rep=y_rep[valid]
                    d_rep=unlist(lapply(y_rep,.diss_d0_d1_L2,.select_domain(v,v_dom,use0,use1),w,alpha)) #The dissimilarity between the sub-curve y_rep and the motif v is calculated using the .diss_d0_d1_L2(returns a vector)
                    d_rep_R=c(FALSE,d_rep<=R,FALSE) #each element indicates whether the dissimilarity is less than or equal to a threshold R.
                    diff_d_rep_R=diff(d_rep_R) #The differences in d_rep_R are computed using the diff function to identify the transitions from FALSE to TRUE and vice versa.
                    start=which(diff_d_rep_R==1) # compute different allignment intervals
                    end=which(diff_d_rep_R==(-1))-1
                    SD_motif=mapply(function(start,end){
                                      index=(start:end)[which.min(d_rep[start:end])] #For each alignment interval, the index with the minimum dissimilarity is selected,
                                      return(c(s_rep[index],d_rep[index]))
                                    },start,end) #The resulting SD_motif list contains the starting point and corresponding dissimilarity values for each valid alignment interval.
                    return(SD_motif)
                  })
  if(!is.null(unlist(SD_motif))){
    v_occurrences=cbind(rep(seq_along(SD_motif),unlist(lapply(SD_motif,length))/2), # Calculates the length of each element in the list SD_motif and divides the lengths by 2, as each element in SD_motif represents both a shift and a dissimilarity value.
                        matrix(unlist(SD_motif),ncol=2,byrow=TRUE))                 # It assigns each shift-dissimilarity pair to its corresponding curve number using rep. Creates a matrix from the flattened SD_motif data, with two columns for shift and dissimilarity values.
    row.names(v_occurrences)=NULL                                                   # combines the two results into a single matrix where the first column is the curve number, and the second and third columns are the shift and dissimilarity values, respectively.
    colnames(v_occurrences)=c('curve','shift','diss')
  }else{
    v_occurrences=c()
  }
  return(v_occurrences)
}


cluster_candidate_motifs <- function(filter_candidate_motifs_results,motif_overlap=0.6,
                                     k_knn=3,votes_knn_Rall=0.5,votes_knn_Rm=0.5,worker_number=NULL){
  # Determine a global radius Rall (based on distances between all motifs and all curves), 
  # and cluster candidate motifs based on their distance, requiring groups to be more than 2*Rall apart. 
  # Then for each group m=1,...,M, determine a group-specific radius Rm (based on distances between motifs of the same group and all curves). 
  # filter_candidate_motifs_results: output of filter_candidate_motifs function.
  # motif_overlap: minimum overlap required between candidate motifs in their distance computation, in % of the shortest motif.
  # k_knn: number of neighbors to be used in k-nearest neighbors, in order to determine Rall and Rm when the two groups (curves with/without motif) are not separated.
  # votes_knn_Rall: threshold on the percentage of votes for class 1 (curve has the motif) in k-nearest neighbors, in order to determine Rall.
  # votes_knn_Rm: threshold on the percentage of votes for class 1 (curve has the motif) in k-nearest neighbors, in order to determine Rm.
  # plot: if TRUE, histohrams and dendrogram are drawn. 
  # ask: if TRUE the user is prompted before a new plot is drawn.
  # worker_number: number of CPU cores to be used for parallelization (default number of CPU cores). If worker_number=1, the function is run sequentially. 
  
  ### set parallel jobs ###################################################################################
  library(parallel)
  core_number <- detectCores()
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
    cl_search=makeCluster(worker_number,timeout=60*60*24*30)
    on.exit(stopCluster(cl_search))
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
    pred_knn=knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=k_knn,prob=TRUE)
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
      pred_knn=knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=k_knn,prob=TRUE)
      R_m[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=votes_knn_Rm)[1]]
    }
  }
  }
  
  ### output ###############################################################################################
  return(c(list(VV_D=VV_D,VV_S=VV_S,k_knn=k_knn,votes_knn_Rall=votes_knn_Rall,R_all=R_all,hclust_res=hclust_res,votes_knn_Rm=votes_knn_Rm,R_m=R_m),filter_candidate_motifs_results)) 
}


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


motifs_search <- function(cluster_candidate_motifs_results,
                          R_all=cluster_candidate_motifs_results$R_all,R_m=NULL,
                          use_real_occurrences=FALSE,length_diff=Inf,worker_number=NULL){
  # Find occurrences of the candidate motifs in the curves and sort them according to their frequencies and radius.
  # In each group (as defined by cutting the dendrogram at high 2*Rall), we choose the motif 
  # with highest frequency and lower mean dissimilarity (the one ranking best in both dimensions). 
  # Additional motifs can be chosen in a group, if their lengths differ enough from the length of the first motif chosen.
  # A candidate motif matches a piece of curve if their dissimilarity is less than the corresponding R_m.
  # cluster_candidate_motifs_results: output of cluster_candidate_motifs function.
  # R_all: global radius, used to cut the dendrogram (requiring groups to be more than 2*Rall apart).
  # R_m: vector with group-specific radii, used to find motif occurrences. The length of the vector must match the number of clusters
  #      obtained cutting the dendrogram at height 2*Rall. If NULL, Rm is determined in each group (based on distances between motifs of the same group and all curves).
  # use_real_occurrences: if TRUE, find occurrences for all candidate motifs and uses real frequency and mean dissimilarity to choose motifs 
  #                       in groups (more accurate, but time consuming). Otherwise, uses approximate frequency and mean dissimilarity (default).
  # length_diff: minimum difference in length among motifs of the same group, required in ordered to keep more than one motif, in % of the most frequent motif.
  # worker_number: number of CPU cores to be used for parallelization (default number of CPU cores). If worker_number=1, the function is run sequentially. 
  
  ### set parallel jobs ###################################################################################
  library(parallel)
  core_number <- detectCores()
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
    cl_search=makeCluster(worker_number,timeout=60*60*24*30)
    on.exit(stopCluster(cl_search))
  }else{
    cl_search=NULL
  }

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
  
  ### cut dendrogram and check group-specific radius Rm ##################################################
  if(is.null(cluster_candidate_motifs_results$hclust_res)){
    V_hclust=1
  }else{
    V_hclust=cutree(cluster_candidate_motifs_results$hclust_res,h=2*R_all) # cut at high 2*R_all
  }
  n_hclust=max(V_hclust)
  
  # re-compute or check group-specific radius Rm
  if(is.null(R_m)){
    R_m=rep(NA,n_hclust)
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_length_i=V_length[index_i]
      V_P_i=as.matrix(cluster_candidate_motifs_results$P_clean[,index_i])
      V_D_i=as.matrix(cluster_candidate_motifs_results$D_clean[,index_i])
      
      dataR=data.frame(P=as.vector(V_P_i),D=as.vector(V_D_i))
      if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
        R_m[i_hclust]=max(dataR$D[dataR$P==1])
      }else{
        D_new=seq(0,max(V_D_i),length.out=10000)
        pred_knn=knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=cluster_candidate_motifs_results$k_knn,prob=TRUE)
        R_m[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=cluster_candidate_motifs_results$votes_knn_Rm)[1]]
      }
    }
  }
  if(length(R_m)!=n_hclust)
    stop(paste0('The length of the vector R_m must match the number of clusters: ',n_hclust))
  
  ### select candidate motifs and find occurrences #######################################################
  if(use_real_occurrences){
    ### find candidate motifs in the curves ####################################################################
    c_k=floor(V_length*(1-max_gap))
    c_k[c_k<cluster_candidate_motifs_results$c]=cluster_candidate_motifs_results$c[c_k<cluster_candidate_motifs_results$c]
    V_R_m=R_m[V_hclust]
    # find occurrences
    V_occurrences=.mapply_custom(cl_search,.find_occurrences,
                                 V,V_R_m,c_k,MoreArgs=list(Y=Y,alpha=alpha,w=w,use0,use1),SIMPLIFY=FALSE)
    not_null=which(!unlist(lapply(V_occurrences,is.null)))
    V=V[not_null]
    V_length=V_length[not_null]
    V_R_m=V_R_m[not_null]
    V_hclust=V_hclust[not_null]
    
    ### select candidate motifs in each group ##############################################################
    select=c()
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      
      V_frequencies_i=unlist(lapply(V_occurrences[index_i],nrow)) # real frequency 
      V_mean_diss_i=unlist(lapply(V_occurrences[index_i],function(x) mean(x[,3]))) # real average distance (actual radius)
      # order based of frequency and average distance
      V_order_i=order(rank(-V_frequencies_i)+rank(V_mean_diss_i)) # sum of ranks in the two dimensions
      V_frequencies_i=V_frequencies_i[V_order_i]
      V_mean_diss_i=V_mean_diss_i[V_order_i]
      index_i_ordered=index_i[V_order_i]
      V_length_i=V_length[index_i_ordered]
      
      # select motifs to keep
      keep=rep_len(TRUE,length(index_i))
      for(i in seq_len(length(index_i))){
        if(keep[i]){
          select=c(select,index_i_ordered[i])
          keep[i]=FALSE
          # motifs with length different enough from length of selected motif
          length_diff_enough=keep&((V_length_i<(V_length_i[i]*(1-length_diff)))|(V_length_i>(V_length_i[i]*(1+length_diff))))
          keep[!length_diff_enough]=FALSE
        }
      }
    }
    V=V[select]
    V_length=V_length[select]
    V_R_m=V_R_m[select]
    V_occurrences=V_occurrences[select]
    
    V_frequencies=unlist(lapply(V_occurrences,nrow)) # real frequency
    V_mean_diss=unlist(lapply(V_occurrences,function(x) mean(x[,3]))) # real average distance (actual radius)
    # sort final motifs based of frequency and average distance
    V_order=order(rank(-V_frequencies)+rank(V_mean_diss)) # sum of ranks in the two dimensions
    V=V[V_order]
    V_occurrences=V_occurrences[V_order]
    V_length=V_length[V_order]
    V_R_m=V_R_m[V_order]
    V_frequencies=V_frequencies[V_order]
    V_mean_diss=V_mean_diss[V_order]
    
    index_final=not_null[select][V_order]
  }else{
    ### select candidate motifs in each group ##############################################################
    select=c()
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_D_i=as.matrix(filter_candidate_motifs_results$D_clean[,index_i])
      
      V_frequencies_approx_i=colSums(V_D_i<=R_m[i_hclust]) # approximate frequency 
      V_mean_diss_approx_i=apply(V_D_i,2,function(x) mean(x[x<=R_m[i_hclust]])) # approximate average distance (actual radius)
      # order based of frequency and average distance
      V_order_i=order(rank(-V_frequencies_approx_i)+rank(V_mean_diss_approx_i)) # sum of ranks in the two dimensions
      V_frequencies_approx_i=V_frequencies_approx_i[V_order_i]
      V_mean_diss_approx_i=V_mean_diss_approx_i[V_order_i]
      index_i_ordered=index_i[V_order_i]
      V_length_i=V_length[index_i_ordered]
      
      # select motifs to keep
      keep=rep_len(TRUE,length(index_i))
      for(i in seq_len(length(index_i))){
        if(keep[i]){
          select=c(select,index_i_ordered[i])
          keep[i]=FALSE
          # motifs with length different enough from length of selected motif
          length_diff_enough=keep&((V_length_i<(V_length_i[i]*(1-length_diff)))|(V_length_i>(V_length_i[i]*(1+length_diff))))
          keep[!length_diff_enough]=FALSE
        }
      }
    }
    V=V[select]
    V_length=V_length[select]
    V_R_m=R_m[V_hclust[select]]
    
    ### find candidate motifs in the curves ####################################################################
    c_k=floor(V_length*(1-max_gap))
    c_k[c_k<cluster_candidate_motifs_results$c[select]]=cluster_candidate_motifs_results$c[select][c_k<cluster_candidate_motifs_results$c[select]]
    # find occurrences
    V_occurrences=.mapply_custom(cl_search,.find_occurrences,
                                 V,V_R_m,c_k,MoreArgs=list(Y=Y,alpha=alpha,w=w,use0,use1),SIMPLIFY=FALSE)
    not_null=which(!unlist(lapply(V_occurrences,is.null)))
    V=V[not_null]
    V_occurrences=V_occurrences[not_null]
    V_length=V_length[not_null]
    V_R_m=V_R_m[not_null]
    V_frequencies=unlist(lapply(V_occurrences,nrow)) # real frequency
    V_mean_diss=unlist(lapply(V_occurrences,function(x) mean(x[,3]))) # real average distance (actual radius)
    # sort final motifs based of frequency and average distance
    V_order=order(rank(-V_frequencies)+rank(V_mean_diss)) # sum of ranks in the two dimensions
    V=V[V_order]
    V_occurrences=V_occurrences[V_order]
    V_length=V_length[V_order]
    V_R_m=V_R_m[V_order]
    V_frequencies=V_frequencies[V_order]
    V_mean_diss=V_mean_diss[V_order]
    
    index_final=select[not_null][V_order]
  }
  
  ### output ##################################################################################################
  if(cluster_candidate_motifs_results$diss=='d0_L2'){
    V0=lapply(V,function(v) v$v0)
    V1=cluster_candidate_motifs_results$V1_clean[index_final]
  }else if(cluster_candidate_motifs_results$diss=='d1_L2'){
    V0=cluster_candidate_motifs_results$V0_clean[index_final]
    V1=lapply(V,function(v) v$v1)
  }else{
    V0=lapply(V,function(v) v$v0)
    V1=lapply(V,function(v) v$v1)
  }
  return(list(V0=V0,V1=V1,
              V_length=V_length,V_occurrences=V_occurrences,V_frequencies=V_frequencies,V_mean_diss=V_mean_diss,
              Y0=cluster_candidate_motifs_results$Y0,Y1=cluster_candidate_motifs_results$Y1,R_motifs=V_R_m))
}


motifs_search_plot <- function(motifs_search_results,ylab='',freq_threshold=5,top_n='all',plot_curves=TRUE){
  # Plot the results of motifs_search.
  # motifs_search_results: output of motifs_search function.
  # ylab: a vector of length d, with the titles for the y axis for each dimension.
  # freq_threshold: plot only motifs with frequency at least equal to freq_threshold.
  # top_n: if 'all', plot all motifs found. If top_n is an integer, then all top top_n motifs are plotted.
  # plot_curves: if TRUE, plot all the curves with coloured motifs.
  
  ### check input ############################################################################################
  # check freq_threshold
  if(max(motifs_search_results$V_frequencies)<freq_threshold)
    stop('There are no motifs with frequency at least equal to freq_threshold.')
  # check top_n
  if(top_n!='all'){
    if(length(top_n)!=1)
      stop('top_n not valid.')
    if(top_n%%1!=0)
      stop('top_n should be an integer.')
    if(top_n<1)
      stop('top_n should be at least 1.')
  }

  ### select motifs to plot ##################################################################################
  d=ncol(motifs_search_results$Y0[[1]])
  N=length(motifs_search_results$Y0)
  index_plot=which(motifs_search_results$V_frequencies>=freq_threshold)
  if(top_n!='all'){
    if(length(index_plot)>top_n)
      index_plot=index_plot[seq_len(top_n)]
  }
  K=length(index_plot)
  V0=motifs_search_results$V0[index_plot]
  V1=motifs_search_results$V1[index_plot]
  V_dom=lapply(V0,function(v) rowSums(!is.na(v))!=0)
  V_length=motifs_search_results$V_length[index_plot]
  V_occurrences=motifs_search_results$V_occurrences[index_plot]
  V_frequencies=motifs_search_results$V_frequencies[index_plot]
  V_mean_diss=motifs_search_results$V_mean_diss[index_plot]
  R_motifs=motifs_search_results$R_motifs[index_plot]

  ### plot motifs ############################################################################################
  layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
  lapply(seq_len(d),
         function(j){
           par(mar=c(3,4,4,2)+0.1)
           plot(V0[[1]][,j],type='l',col=rainbow(K),lwd=5,lty=1,main=ylab[j],xlim=c(1,max(V_length)),
                ylab=ylab[j],ylim=c(min(unlist(V0)),max(unlist(V0))))
           mapply(function(v,k) points(v[,j],type='l',col=rainbow(K)[k+1],lwd=5,lty=1,ylab=ylab),
                  V0[-1],seq_len(K-1))
           par(mar=c(0,0,0,0))
           return()})
  plot.new()
  legend('left',paste('motif',seq_len(K)),col=rainbow(K),lwd=7,lty=1,bty="n",xpd=TRUE)
  if(!is.null(V1[[1]])){
    layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             plot(V1[[1]][,j],type='l',col=rainbow(K),lwd=5,lty=1,main=paste(ylab[j],'derivative'),xlim=c(1,max(V_length)),
                  ylab=ylab[j],ylim=c(min(unlist(V1)),max(unlist(V1))))
             mapply(function(v,k) points(v[,j],type='l',col=rainbow(K)[k+1],lwd=5,lty=1,ylab=ylab),
                    V1[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             return()})
    plot.new()
    legend('left',paste('motif',seq_len(K)),col=rainbow(K),lwd=7,lty=1,bty="n",xpd=TRUE)
  }

  ### plot motifs with matched curves ########################################################################
  if(is.null(V1[[1]])){
    mapply(function(v,v_dom,v_occurrences,v_frequencies,k,R_motif){
             Y_inters_k=mapply(function(y,s_k_i,v_dom){
                          v_len=length(v_dom)
                          Y_inters_k=as.matrix(as.matrix(y[s_k_i-1+seq_len(v_len),])[v_dom,])
                          return(Y_inters_k)},
                          motifs_search_results$Y0[v_occurrences[,'curve']],v_occurrences[,'shift'],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
             layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
             lapply(seq_len(d),
                    function(j){
                      par(mar=c(3,4,4,2)+0.1)
                      y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
                      y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                      matplot(y_plot,type='l',col=v_occurrences[,'curve']+1,lwd=round(-4/R_motif*v_occurrences[,'diss']+5,2),
                              lty=1,ylab=ylab[j],main=paste0('Motif ',k,' (',v_frequencies,' occurrences) - ',ylab[j]))
                      points(v[,j],type='l',col='black',lwd=7,lty=1)
                      par(mar=c(0,0,0,0))
                      plot.new()
                      legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
                      })
             return()},V0,V_dom,V_occurrences,V_frequencies,seq_len(K),R_motifs)
  }else{
    mapply(function(v0,v1,v_dom,v_occurrences,v_frequencies,k,R_motif){
             Y0_inters_k=mapply(
                           function(y,s_k_i,v_dom){
                             v_len=length(v_dom)
                             Y_inters_k=as.matrix(as.matrix(y[s_k_i-1+seq_len(v_len),])[v_dom,])
                             return(Y_inters_k)},
                           motifs_search_results$Y0[v_occurrences[,'curve']],v_occurrences[,'shift'],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
             Y1_inters_k=mapply(
                           function(y,s_k_i,v_dom){
                             v_len=length(v_dom)
                             Y_inters_k=as.matrix(as.matrix(y[s_k_i-1+seq_len(v_len),])[v_dom,])
                             return(Y_inters_k)},
                           motifs_search_results$Y1[v_occurrences[,'curve']],v_occurrences[,'shift'],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
             layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
             lapply(seq_len(d),
                    function(j){
                      par(mar=c(3,4,4,2)+0.1)
                      y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y0_inters_k))
                      y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                      matplot(y_plot,type='l',col=v_occurrences[,'curve']+1,lwd=round(-4/R_motif*v_occurrences[,'diss']+5,2),
                              lty=1,ylab=ylab[j],main=paste0('Motif ',k,' (',v_frequencies,' occurrences) - ',ylab[j]))
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
                      matplot(y_plot,type='l',col=v_occurrences[,'curve']+1,lwd=round(-4/R_motif*v_occurrences[,'diss']+5,2),
                              lty=1,ylab=ylab[j],main=paste0('Motif ',k,' (',v_frequencies,' occurrences) - ',ylab[j],' derivative'))
                      points(v1[,j],type='l',col='black',lwd=7,lty=1)
                      par(mar=c(0,0,0,0))
                      plot.new()
                      legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
                    })
            return()},V0,V1,V_dom,V_occurrences,V_frequencies,seq_len(K),R_motifs)
    }

  ### plot curves with motifs ################################################################################
  if(plot_curves){
    if(is.null(motifs_search_results$Y1[[1]])){
      mapply(function(y0,i){
               s_i=lapply(V_occurrences,function(occurrences) occurrences[occurrences[,'curve']==i,'shift'])
               motifs_in_curve=rep(seq_len(K),unlist(lapply(s_i,length)))
               s_i=unlist(s_i)
               Y_inters_k=mapply(function(v_dom,s_i_k,y){
                                   v_len=length(v_dom)
                                   Y_inters_k=matrix(NA,nrow=length(v_dom),ncol=d)
                                   Y_inters_k[v_dom,]=as.matrix(as.matrix(y[s_i_k-1+seq_len(v_len),])[v_dom,])
                                   return(Y_inters_k)},
                                 V_dom[motifs_in_curve],s_i,MoreArgs=list(y0),SIMPLIFY=FALSE)
               layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
               lapply(seq_len(d),
                      function(j){
                        par(mar=c(3,4,4,2)+0.1)
                        plot(y0[,j],type='l',main=paste('Region',i,'-',ylab[j]),ylab=ylab[j])
                        for(k in seq_along(motifs_in_curve)){
                          lines(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]),Y_inters_k[[k]][,j],col=rainbow(K)[motifs_in_curve[k]],lwd=5)
                        }
                      })
               plot.new()
               if(length(motifs_in_curve)==0){
                 legend_text=''
               }else{
                 legend_text=paste('motif',unique(motifs_in_curve))
               }
               legend('left',legend_text,col=rainbow(K)[unique(motifs_in_curve)],lwd=7,lty=1,bty="n",xpd=TRUE,title='Motifs')
               return()},motifs_search_results$Y0,seq_len(N))
    }else{
      mapply(function(y0,y1,i){
                s_i=lapply(V_occurrences,function(occurrences) occurrences[occurrences[,'curve']==i,'shift'])
                motifs_in_curve=rep(seq_len(K),unlist(lapply(s_i,length)))
                s_i=unlist(s_i)
                Y0_inters_k=mapply(function(v_dom,s_i_k,y){
                                     v_len=length(v_dom)
                                     Y_inters_k=matrix(NA,nrow=length(v_dom),ncol=d)
                                     Y_inters_k[v_dom,]=as.matrix(as.matrix(y[s_i_k-1+seq_len(v_len),])[v_dom,])
                                     return(Y_inters_k)},
                                   V_dom[motifs_in_curve],s_i,MoreArgs=list(y0),SIMPLIFY=FALSE)
                Y1_inters_k=mapply(function(v_dom,s_i_k,y){
                                     v_len=length(v_dom)
                                     Y_inters_k=matrix(NA,nrow=length(v_dom),ncol=d)
                                     Y_inters_k[v_dom,]=as.matrix(as.matrix(y[s_i_k-1+seq_len(v_len),])[v_dom,])
                                     return(Y_inters_k)},
                                   V_dom[motifs_in_curve],s_i,MoreArgs=list(y1),SIMPLIFY=FALSE)
                layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
                lapply(seq_len(d),
                       function(j){
                         par(mar=c(3,4,4,2)+0.1)
                         plot(y0[,j],type='l',main=paste('Region',i,'-',ylab[j]),ylab=ylab[j],xlab='')
                         for(k in seq_along(motifs_in_curve)){
                           lines(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]),Y0_inters_k[[k]][,j],col=rainbow(K)[motifs_in_curve[k]],lwd=5)
                         }
                       })
                plot.new()
                if(length(motifs_in_curve)==0){
                  legend_text=''
                }else{
                  legend_text=paste('motif',unique(motifs_in_curve))
                }
                legend('left',legend_text,col=rainbow(K)[unique(motifs_in_curve)],lwd=7,lty=1,bty="n",xpd=TRUE,title='Motifs')
                lapply(seq_len(d),
                       function(j){
                         par(mar=c(3,4,4,2)+0.1)
                         plot(y1[,j],type='l',main=paste('Region',i,'-',paste(ylab[j],'derivative')),ylab=ylab[j])
                         for(k in seq_along(motifs_in_curve)){
                           lines(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]),Y1_inters_k[[k]][,j],col=rainbow(K)[motifs_in_curve[k]],lwd=5)
                         }
                       })
                plot.new()
                if(length(motifs_in_curve)==0){
                  legend_text=''
                }else{
                  legend_text=paste('motif',unique(motifs_in_curve))
                }
                legend('left',legend_text,col=rainbow(K)[unique(motifs_in_curve)],lwd=7,lty=1,bty="n",xpd=TRUE,title='Motifs')
                return()},motifs_search_results$Y0,motifs_search_results$Y1,seq_len(N))
    }
  }

  return()
}

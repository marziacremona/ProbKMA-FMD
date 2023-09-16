library(microbenchmark)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # mi metto nel path corrente 

source("../ProbKMA-FMD_functions.r") # "../ProbKMA-FMD_functions_modified.r" se vuoi quello modificato 

# comparison d0_d1_L2
use0 = TRUE
use1 = TRUE
alpha = 0.5
max_gap = 0
# we need data.RData and resLeftRight_arguments.RData
# resLeftRight_arguments.RData are arguments with d0_d1_L2
Y = mapply(function(y0,y1) list(y0=y0,y1=y1),output_prof$Y0,output_prof$Y1,SIMPLIFY=FALSE) # better to be takes using debug mode
m = 2 
deltaJk_elong = 0.05
w = 1
elong_prof <- function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
  if(length(len_elong_k)==0){
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
  s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k) s_k-len_elong_k),(length(len_elong_k)+1):1)[-1]
  v_dom_elong_left_right=unlist(lapply(c(0,len_elong_k),
                                       function(len_elong_k_left)
                                         lapply(c(0,len_elong_k[len_elong_k<=(max(len_elong_k)-len_elong_k_left)]),
                                                function(len_elong_k_right) 
                                                  c(rep_len(TRUE,len_elong_k_left),v_dom_k,rep_len(TRUE,len_elong_k_right)))),
                                recursive=FALSE)[-1]
  v_elong_left_right=mapply(.compute_motif,v_dom_elong_left_right,s_k_elong_left_right,
                            MoreArgs=list(p_k,Y,m,use0,use1),SIMPLIFY=FALSE)
  start_with_NA=unlist(lapply(v_elong_left_right,length))>2
  v_elong_left_right=v_elong_left_right[!start_with_NA]
  s_k_elong_left_right=s_k_elong_left_right[!start_with_NA]
  Jk_before=.compute_Jk(v_new_k,s_k,p_k,Y,alpha,w,m,use0=use0,use1=use1)
  c_k_after=floor(unlist(lapply(lapply(v_elong_left_right,.domain,use0),length))*(1-max_gap))
  c_k_after[c_k_after<c]=c
  Jk_after=unlist(mapply(.compute_Jk,v_elong_left_right,s_k_elong_left_right,c_k_after,
                         MoreArgs=list(p_k=p_k,Y=Y,alpha=alpha,w=w,m=m,keep_k=keep_k,use0=use0,use1=use1)))
  best_elong=which.min((Jk_after-Jk_before)/Jk_before)
  if(length(best_elong)>0){
    elongate=((Jk_after-Jk_before)/Jk_before)[best_elong]<deltaJk_elong
  }else{
    elongate=FALSE
  }
  if(elongate){
    return(list(v_new=v_elong_left_right[[best_elong]],
                v_dom=v_dom_elong_left_right[[best_elong]],
                s_k=s_k_elong_left_right[[best_elong]]))
  }else{
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
}
result <- elong_prof(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c)
my_result = elongation_rcpp(v_new_k, v_dom_k, s_k, p_k, len_elong_k,keep_k, c, .compute_Jk, .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong)

setequal(result$v_new,my_result$v_new)
setequal(result$v_dom,my_result$v_dom)
setequal(result$s_k,my_result$s_k)

# comparison d0_L2
use0 = TRUE
use1 = FALSE
alpha = 0.0
max_gap = 0
# we need data.RData and resLeftRight_arguments_d0L2.RData
Y=lapply(output_prof$Y0,function(y0) list(y0=y0,y1=NULL))
m = 2 
deltaJk_elong = 0.05
w = 1
elong_prof <- function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
  if(length(len_elong_k)==0){
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
  s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k) s_k-len_elong_k),(length(len_elong_k)+1):1)[-1]
  v_dom_elong_left_right=unlist(lapply(c(0,len_elong_k),
                                       function(len_elong_k_left)
                                         lapply(c(0,len_elong_k[len_elong_k<=(max(len_elong_k)-len_elong_k_left)]),
                                                function(len_elong_k_right) 
                                                  c(rep_len(TRUE,len_elong_k_left),v_dom_k,rep_len(TRUE,len_elong_k_right)))),
                                recursive=FALSE)[-1]
  v_elong_left_right=mapply(.compute_motif,v_dom_elong_left_right,s_k_elong_left_right,
                            MoreArgs=list(p_k,Y,m,use0,use1),SIMPLIFY=FALSE)
  start_with_NA=unlist(lapply(v_elong_left_right,length))>2
  v_elong_left_right=v_elong_left_right[!start_with_NA]
  s_k_elong_left_right=s_k_elong_left_right[!start_with_NA]
  Jk_before=.compute_Jk(v_new_k,s_k,p_k,Y,alpha,w,m,use0=use0,use1=use1)
  c_k_after=floor(unlist(lapply(lapply(v_elong_left_right,.domain,use0),length))*(1-max_gap))
  c_k_after[c_k_after<c]=c
  Jk_after=unlist(mapply(.compute_Jk,v_elong_left_right,s_k_elong_left_right,c_k_after,
                         MoreArgs=list(p_k=p_k,Y=Y,alpha=alpha,w=w,m=m,keep_k=keep_k,use0=use0,use1=use1)))
  best_elong=which.min((Jk_after-Jk_before)/Jk_before)
  if(length(best_elong)>0){
    elongate=((Jk_after-Jk_before)/Jk_before)[best_elong]<deltaJk_elong
  }else{
    elongate=FALSE
  }
  if(elongate){
    return(list(v_new=v_elong_left_right[[best_elong]],
                v_dom=v_dom_elong_left_right[[best_elong]],
                s_k=s_k_elong_left_right[[best_elong]]))
  }else{
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
}
result <- elong_prof(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c)
my_result = elongation_rcpp(v_new_k, v_dom_k, s_k, p_k, len_elong_k,keep_k, c, .compute_Jk, .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong)

setequal(result$v_new,my_result$v_new)
setequal(result$v_dom,my_result$v_dom)
setequal(result$s_k,my_result$s_k)

# comparison d1_L2
use0 = FALSE
use1 = TRUE
alpha = 1.0
max_gap = 0
# we need data.RData and resLeftRight_arguments_d1L2.RData
Y = lapply(output_prof$Y1,function(y1) list(y0=NULL,y1=y1))
m = 2 
deltaJk_elong = 0.05
w = 1
elong_prof <- function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
  if(length(len_elong_k)==0){
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
  s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k) s_k-len_elong_k),(length(len_elong_k)+1):1)[-1]
  v_dom_elong_left_right=unlist(lapply(c(0,len_elong_k),
                                       function(len_elong_k_left)
                                         lapply(c(0,len_elong_k[len_elong_k<=(max(len_elong_k)-len_elong_k_left)]),
                                                function(len_elong_k_right) 
                                                  c(rep_len(TRUE,len_elong_k_left),v_dom_k,rep_len(TRUE,len_elong_k_right)))),
                                recursive=FALSE)[-1]
  v_elong_left_right=mapply(.compute_motif,v_dom_elong_left_right,s_k_elong_left_right,
                            MoreArgs=list(p_k,Y,m,use0,use1),SIMPLIFY=FALSE)
  start_with_NA=unlist(lapply(v_elong_left_right,length))>2
  v_elong_left_right=v_elong_left_right[!start_with_NA]
  s_k_elong_left_right=s_k_elong_left_right[!start_with_NA]
  Jk_before=.compute_Jk(v_new_k,s_k,p_k,Y,alpha,w,m,use0=use0,use1=use1)
  c_k_after=floor(unlist(lapply(lapply(v_elong_left_right,.domain,use0),length))*(1-max_gap))
  c_k_after[c_k_after<c]=c
  Jk_after=unlist(mapply(.compute_Jk,v_elong_left_right,s_k_elong_left_right,c_k_after,
                         MoreArgs=list(p_k=p_k,Y=Y,alpha=alpha,w=w,m=m,keep_k=keep_k,use0=use0,use1=use1)))
  best_elong=which.min((Jk_after-Jk_before)/Jk_before)
  if(length(best_elong)>0){
    elongate=((Jk_after-Jk_before)/Jk_before)[best_elong]<deltaJk_elong
  }else{
    elongate=FALSE
  }
  if(elongate){
    return(list(v_new=v_elong_left_right[[best_elong]],
                v_dom=v_dom_elong_left_right[[best_elong]],
                s_k=s_k_elong_left_right[[best_elong]]))
  }else{
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
}
result <- elong_prof(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c)
my_result = elongation_rcpp(v_new_k, v_dom_k, s_k, p_k, len_elong_k,keep_k, c, .compute_Jk, .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong)

setequal(result$v_new,my_result$v_new)
setequal(result$v_dom,my_result$v_dom)
setequal(result$s_k,my_result$s_k)

# time comparison 
times_comparison <- microbenchmark(prof = elong_prof(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c),mia = elongation_rcpp(v_new_k, v_dom_k, s_k, p_k, len_elong_k,keep_k, c, .compute_Jk, .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong))
print(times_comparison)



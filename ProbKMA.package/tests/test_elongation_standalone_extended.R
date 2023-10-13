source("../ProbKMA-FMD_functions.r") # "../ProbKMA-FMD_functions_modified.r" se vuoi quello modificato 
#library(microbenchmark)


# comparison d0_d1_L2
use0 = TRUE
use1 = TRUE
alpha = 0.5
max_gap = 0
set.seed(123)
# we need data_before_elongation.RData
load(paste0('len200_sd0.1.RData')) 
m = 2 
deltaJk_elong = 0.05
w = 1
elong_prof <- function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
  if(length(len_elong_k)==0){
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
  s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k_i) s_k-len_elong_k_i),(length(len_elong_k)+1):1)[-1]
  save(s_k_elong_left_right, file = "s_k_elong.RData")
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

result <- elong_prof(V_new[[1]],V_dom[[1]],S_k[[1]],P_k[[1]],len_elong[[1]],Keep_k[[1]],c[[1]])
my_result = elongation_rcpp(V_new[[1]],V_dom[[1]],S_k[[1]],P_k[[1]],len_elong[[1]],Keep_k[[1]],c[[1]], .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong)

setequal(result$v_new,my_result$v_new)

result <- elong_prof(V_new[[2]],V_dom[[2]],S_k[[2]],P_k[[2]],len_elong[[2]],Keep_k[[2]],c[[2]])
my_result = elongation_rcpp(V_new[[2]],V_dom[[2]],S_k[[2]],P_k[[2]],len_elong[[2]],Keep_k[[2]],c[[2]], .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong)

setequal(result$v_new,my_result$v_new)

#library(microbenchmark)

#risultati <- microbenchmark(prima = elongation_rcpp_prima(V_new[[2]],V_dom[[2]],S_k[[2]],P_k[[2]],len_elong[[2]],Keep_k[[2]],c[[2]],.domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong),
                            dopo = elongation_rcpp_dopo(V_new[[2]],V_dom[[2]],S_k[[2]],P_k[[2]],len_elong[[2]],Keep_k[[2]],c[[2]],.domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong))
#print(risultati)


diss = 'd1_L2' 
alpha = 1
max_gap = 0 
use0 = FALSE
use1 = TRUE
m = 2 
deltaJk_elong = 0.05
w = 1
elong_prof <- function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
  if(length(len_elong_k)==0){
    return(list(v_new=v_new_k,
                v_dom=v_dom_k,
                s_k=s_k))
  }
  s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k_i) s_k-len_elong_k_i),(length(len_elong_k)+1):1)[-1]
  save(s_k_elong_left_right, file = "s_k_elong.RData")
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

result <- elong_prof(V_new[[2]],V_dom[[2]],S_k[[2]],P_k[[2]],len_elong[[2]],Keep_k[[2]],c[[2]])
my_result = elongation_rcpp(V_new[[2]],V_dom[[2]],S_k[[2]],P_k[[2]],len_elong[[2]],Keep_k[[2]],c[[2]], .domain, .compute_motif, use0, use1,w, alpha,max_gap, Y,m, deltaJk_elong)
setequal(my_result, result)

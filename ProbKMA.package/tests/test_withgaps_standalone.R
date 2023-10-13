source("../ProbKMA-FMD_functions.r")
# env_before_gaps_part.RData to be used

V_dom_false <- V_dom
for (i in 1:length(V_dom_false[[1]])) {
  if (i %% 5 == 0) {
    V_dom_false[[1]][i] <- FALSE
    V_dom_false[[2]][i] <- FALSE
  }
}
save(V_dom_false,file = "dom_with_false.RData")

with_gaps=which(unlist(lapply(V_dom_false,function(V_dom_false) sum(!V_dom_false)!=0)))
if(length(with_gaps)>0){
  V_dom_false_filled=lapply(V_dom_false[with_gaps],function(V_dom_false) rep_len(TRUE,length(V_dom_false)))
  V_filled=mapply(.compute_motif,V_dom_false_filled,S_k[with_gaps],P_k[with_gaps],MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
  Jk_before=mapply(.compute_Jk,
                   V_new[with_gaps],S_k[with_gaps],P_k[with_gaps],
                   MoreArgs=list(Y=Y,alpha=alpha,w=w,m=m,use0=use0,use1=use1))
  Jk_after=mapply(.compute_Jk,
                  V_filled,S_k[with_gaps],P_k[with_gaps],
                  MoreArgs=list(Y=Y,alpha=alpha,w=w,m=m,use0=use0,use1=use1))
  fill=(Jk_after-Jk_before)/Jk_before<deltaJk_elong # output is true, true
  V_dom_false[with_gaps[fill]]=V_dom_false_filled[fill]
  V_new[with_gaps[fill]]=V_filled[fill]
}

V_dom_false_prof <- V_dom_false
V_new_prof <- V_new
rm(V_dom_false)
rm(V_new)
# load dom_with_false.RData
# load env_before_gaps_part.RData


# chiamo la mia funzione corrotta con il return e confronto V_dom_false e V_new
list_output <- elongate_motifs(V_new,
                               V_dom_false,
                               S_k,
                               P_k,
                               Y,
                               w, 
                               m, 
                               use0,
                               use1,
                               alpha,
                               c,
                               c_max, 
                               max_elong, 
                               deltaJk_elong,
                               trials_elong,
                               D,
                               K,
                               max_gap) 
# return List::create(V_dom, V_new,with_gaps, fill_); to be add in elongate_motifs
setequal(V_new_prof, list_output[[2]])
setequal(V_dom_false_prof, list_output[[1]])
setequal(fill, list_output[[4]])
setequal(with_gaps, list_output[[3]] + 1)

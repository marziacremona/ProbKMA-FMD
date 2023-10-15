
compute_motif <- function(v_dom,s_k,p_k,Y,m,use0,use1){
  # Compute the new motif v_new.
  # v_dom: TRUE for x in supp(v).
  # s_k: shift vector for motif k.
  # p_k: membership vector for motif k.
  # Y: list of N lists of two elements, Y0=y_i(x), Y1=y'_i(x), matrices with d columns, for d-dimensional curves.
  
  .compute_v_new <- function(Y_inters_k,Y_inters_supp,v_dom,v_len,p_k,d,m){
    v_new=matrix(NA,nrow=v_len,ncol=d)
    if(length(Y_inters_k)==1){
      v_new[v_dom,]=Y_inters_k[[1]]
      return(v_new)
    }
    Y_inters_supp=Reduce(rbind,Y_inters_supp)
    coeff_k=p_k[p_k>0]^m/(rowSums(Y_inters_supp)) # NB: divide for the length of the interval, not for the squared length!
    coeff_x=colSums(Y_inters_supp*coeff_k)
    coeff_x[colSums(Y_inters_supp)==0]=NA
    v_new[v_dom,]=Reduce('+',mapply('*',Y_inters_k,coeff_k,SIMPLIFY=FALSE))/coeff_x
    return(v_new)
  }
  
  if(sum(p_k)==0){
    stop('Motif with no members! Degenerate cluster!')
  }
  v_len=length(v_dom)
  d=unlist(lapply(Y[[1]],ncol))[1] # dimension of curves
  Y_inters_k=mapply(function(y,s_k_i,d,use0,use1){
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
    return(select_domain(y_inters_k,v_dom,use0,use1))
  },Y[p_k>0],s_k[p_k>0],MoreArgs=list(d,use0,use1),SIMPLIFY=FALSE)
  Y_inters_supp=lapply(Y_inters_k,domain,use0)
  Y_inters_k=mapply(function(y_inters_k,y_inters_supp,use0,use1){
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
  range_v_new=range(which(domain(v_new,use0)))
  v_dom_new=rep(FALSE,v_len)
  v_dom_new[(range_v_new[1]):(range_v_new[2])]=TRUE
  v_new=select_domain(v_new,v_dom_new,use0,use1)
  if(range_v_new[1]>1){
    return(c(v_new,list(shift=range_v_new[1]-1)))
  }
  return(v_new)
}

find_diss_aligned <- function(y,v,alpha,w,aligned,d,use0,use1)
{
  out <- .find_diss_aligned_rcpp(y, v, w, alpha, aligned, d, 
                                 use0, use1, domain, select_domain, 
                                 diss_d0_d1_L2)
}
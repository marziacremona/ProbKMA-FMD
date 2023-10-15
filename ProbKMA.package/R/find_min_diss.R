
#' @export
find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1){
  # Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
  # Return shift and dissimilarity.
  # y: list of two elements y0=y(x), y1=y'(x), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x), matrices with d columns.
  # alpha: weight coefficient between d0.L2 and d1.L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  
  out<-.find_diss(y,v,w,alpha,c_k,d,use0,use1,domain,select_domain,diss_d0_d1_L2)
}


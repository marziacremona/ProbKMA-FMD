
.find_occurrences <- function(v,Y,R,alpha,w,c_k,use0,use1){
  # Find occurrences of a motif in a set of curves (dimesion=d), with dissimilarity lower than R.
  # Return curve id, shift and dissimilarity.
  # v: list of 2 elements, v0, v1, matrices with d columns.
  # Y: list of N lists of two elements, Y0, Y1, matrices with d columns.
  # R: maximum dissimilarity allowed. 
  # alpha: if diss_fun=diss_d0_d1_L2, weight coefficient between d0_L2 and d1_L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.

  v_occurrences <- .find_occurrences_cpp(v,Y,R,alpha,w,c_k,use0,use1,
                                         diss_d0_d1_L2,
                                         domain,
                                         select_domain)
  if(!is.null(result))
  {
    row.names(v_occurrences)=NULL
    colnames(v_occurrences)=c('curve','shift','diss')
  }
  
}
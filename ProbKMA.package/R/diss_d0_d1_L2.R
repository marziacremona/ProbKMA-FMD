#' @title .diss_d0_d1_L2
#'
#' @description Dissimilarity index for multidimensional curves (dimension=d).
#' Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
#'
#' @param y list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param alpha weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
#' @export
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
diss_d0_d1_L2 <- function(y,v,w,alpha){
  # Dissimilarity index for multidimensional curves (dimension=d).
  # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
  # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
  
  if(alpha==0){
    return(diss_L2(y[[1]],v[[1]],w))
  }else if(alpha==1){
    return(diss_L2(y[[2]],v[[2]],w))
  }else{
    return((1-alpha)*diss_L2(y[[1]],v[[1]],w)+alpha*diss_L2(y[[2]],v[[2]],w))
  }
 
}
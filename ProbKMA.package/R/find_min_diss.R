#' @title .find_min_diss
#'
#' @description Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
#'
#' @param y list of two elements y0=y(x), y1=y'(x), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x), matrices with d columns.
#' @param alpha weight coefficient between d0.L2 and d1.L2.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param c_k minimum length of supp(y_shifted) and supp(v) intersection.
#' @return Shift warping and dissimilarity
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1){

  out<-.find_diss(y,v,w,alpha,c_k,d,use0,use1,domain,select_domain,diss_d0_d1_L2)
}


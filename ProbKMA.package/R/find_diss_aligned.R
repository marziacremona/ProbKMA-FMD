#' @title .find_diss
#'
#' @description Find dissimilarity between multidimensional curves (dimension=d), without alignment unless
#' their lengths are different. To be used by probKMA_silhouette fucntion.
#'
#' @param y list of two elements y0=y(x), y1=y'(x), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x), matrices with d columns.
#' @param alpha weight coefficient between d0.L2 and d1.L2.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param aligned if TRUE, curves are already aligned. Else, the shortest curve is aligned inside the longest.
#' @return shift and dissimilarity
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
find_diss_aligned <- function(y,v,alpha,w,aligned,d,use0,use1)
{
  out <- .find_diss_aligned_rcpp(y, v, w, alpha, aligned, d, 
                                 use0, use1, domain, select_domain, 
                                 diss_d0_d1_L2)
}
#' @title diss_L2
#'
#' @description Dissimilarity index for multidimensional curves (dimension=d)
#' L2 distance with normalization on common support.
#'
#' @param y list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
diss_L2 <- function(y,v,w){
 
  sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
}
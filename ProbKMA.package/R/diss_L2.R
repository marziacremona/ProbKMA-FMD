
diss_L2 <- function(y,v,w){
  # Dissimilarity index for multidimensional curves (dimension=d).
  # L2 distance with normalization on common support.
  # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  
  sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
}
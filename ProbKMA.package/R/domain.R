#' @title .domain
#'
#' @description Length of the domain of the curve without NA
#'
#' @param v list of two elements: v0=v(x), v1=v'(x), matrices with d columns.
#' @param use0 boolean: use0=TRUE, v0 is considered else v1 is considered.
#' @return vector of boolean whose elements are TRUE iff at least one element in the row of v is not NA
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
domain <- function(v,use0){
  if(use0){
    rowSums(!is.na(v[[1]]))!=0
  }else{
    rowSums(!is.na(v[[2]]))!=0
  }
}
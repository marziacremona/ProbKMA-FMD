#' @title .select_domain
#'
#' @description Select the portion of a motif without NA
#'
#' @param v list of two elements: v0=v(x), v1=v'(x), matrices with d columns.
#' @param v_dom boolean vector: domain of v0 and v1.
#' @param use0 boolean: if use0=TRUE, v0 is selected.
#' @param use1 boolean: if use1=TRUE, v1 is selected. 
#' 
#' @return v list: containing the selected portion of v and v'.
#' 
#' @author Marzia Angela Cremona  & Francesca Chiaromonte
#' @export
select_domain <- function(v,v_dom,use0,use1){
  if(use0)
    v[[1]]=as.matrix(v[[1]][v_dom,])
  if(use1)
    v[[2]]=as.matrix(v[[2]][v_dom,])
  return(v)
}
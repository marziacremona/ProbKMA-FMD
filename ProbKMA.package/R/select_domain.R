
#' @export
select_domain <- function(v,v_dom,use0,use1){
  if(use0)
    v[[1]]=as.matrix(v[[1]][v_dom,])
  if(use1)
    v[[2]]=as.matrix(v[[2]][v_dom,])
  return(v)
}

#' Sum of vector elements
#' 
#' @description
#' `sum` returns the sum of all the values present in its arguments.
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
#' @export
compute_jk <- function(v,s_k, p_k,Y,
                       alpha,w,m,use0, use1,
                       domain,select_domain,diss_d0_d1_L2,
                       c_k = NULL,
                       keep_k = NULL)
{
  out<-.Call('_ProbKMA_package_compute_Jk_rcpp', PACKAGE = 'probKMA.package', 
             v,s_k,p_k,Y,alpha,w,m,
             use0,use1,domain,select_domain,
             diss_d0_d1_L2,c_k = NULL,keep_k = NULL)

}

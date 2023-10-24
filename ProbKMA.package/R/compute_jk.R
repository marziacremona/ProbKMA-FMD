#' @title .compute_Jk
#'
#' @description Compute the objective function J for the motif k.
#'
#' @param v list of two elements: v0=v(x), v1=v'(x), matrices with d columns.
#' @param s_k integer vector: shift vector for motif k.
#' @param p_k numeric vector: membership vector for motif k.
#' @param Y list of N lists of two elements: Y0=y_i(x), Y1=y'_i(x), matrices with d columns, for d-dimensional curves.
#' @param alpha double: weight coefficient between d0_L2 and d1_L2.
#' @param m>1 integer vector: weighting exponent in least-squares functional.
#' @param w numeric vector: weights for the dissimilarity index in the different dimensions (w>0).
#' @param c_k integer vector: minimum length of supp(y_shifted) and supp(v) intersection.
#' @param keep_k boolean vector: check c_k only when keep=TRUE for y_shifted.
#' @return The value of the objective function J
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
compute_jk <- function(v,s_k, p_k,Y,
                       alpha,w,m,use0, use1,
                       domain,select_domain,diss_d0_d1_L2,
                       c_k = NULL,
                       keep_k = NULL)
{
  out<-.compute_Jk_rcpp(v,s_k,p_k,Y,alpha,w,m,
                        use0,use1,domain,select_domain,
                        diss_d0_d1_L2,c_k)
}

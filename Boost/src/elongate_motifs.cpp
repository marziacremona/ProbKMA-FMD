// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// void function which modifies V_new, V_dom, S_k. 
// return the matrix S to be initialized in the R code 
// [[Rcpp::export]]
mat elongate_motifs(List & V_new,
                    List & V_dom,
                    List & S_k,
                    const List & P_k,
                    const List & Y,
                    const vec & w,
                    int m,
                    bool use0,
                    bool use1,
                    double alpha,
                    const ivec & c,
                    int c_max, 
                    double max_elong, 
                    double deltaJk_elong,
                    int trials_elong,
                    const mat & D,
                    unsigned int K,
                    double max_gap) 
{
  
  // with_gaps is a vector of indexes and contains all the indexes of the curves with gaps in the domain
  const unsigned int V_dom_size = V_dom.size();
  IntegerVector with_gaps;
  IntegerVector len_dom(V_dom_size);
  for (unsigned int i = 0; i < V_dom_size; ++i){
    const LogicalVector & v_dom = V_dom[i];
    
    // compute the lengths of the domains
    len_dom[i] = v_dom.size();
    
    if (sum(!v_dom)!= 0) {
      with_gaps.push_back(i);
    }
  }
  const unsigned int with_gaps_size = with_gaps.size(); 
  // if the are motifs with gaps try to fill the gaps and compute the perf_idnex before and after
  if (with_gaps_size > 0){
    
    Function compute_motif(".compute_motif");
    List V_dom_filled(with_gaps_size);
    
    // fill the domains of the motifs with gaps
    for (unsigned int i = 0; i < with_gaps_size; ++i){
      const LogicalVector & v_dom = V_dom[with_gaps[i]];
      V_dom_filled[i] = rep(true,v_dom.size());
    }
    
    // recompute the motifs with the filled domains
    List V_filled(with_gaps_size);
    for (unsigned int i = 0; i < with_gaps_size; ++i){
      V_filled[i] = as<List>(compute_motif(V_dom_filled[i],S_k[with_gaps[i]],P_k[with_gaps[i]],Y,m,use0,use1));
    }
    
    // compute the perf.indexes before and after the filling
    vec Jk_before(with_gaps_size);
    vec Jk_after(with_gaps_size);
    Function compute_Jk(".compute_Jk");
    for (unsigned int i = 0; i < with_gaps_size; ++i){
      Jk_before(i) = as<double>(compute_Jk(V_new[with_gaps[i]],
                                S_k[with_gaps[i]],
                                P_k[with_gaps[i]],
                                Y,
                                alpha,
                                w,
                                m,
                                use0,
                                use1));
      Jk_after(i) = as<double>(compute_Jk(V_filled[i],
                               S_k[with_gaps[i]],
                               P_k[with_gaps[i]],
                               Y,
                               alpha,
                               w,
                               m,
                               use0,
                               use1));
    }
    
    // if filling the domain improves the perf. index over a certain threshold replace the domain and the motifs with the filled one
    const uvec & fill = (Jk_after-Jk_before)/Jk_before < deltaJk_elong;
    for (unsigned int i = 0; i < with_gaps_size; ++i){
      if (fill[i]){
        V_dom[with_gaps[i]] = V_dom_filled[i];
        V_new[with_gaps[i]] = V_filled[i];
      }
    }
   }
    IntegerVector len_max_elong(V_dom_size);
    for (unsigned int i = 0; i < V_dom_size; ++i){
      int temp_elongation = std::floor(len_dom[i]*max_elong);
      len_max_elong[i] = std::min(temp_elongation,c_max - len_dom[i]);
    }
    List len_elong(V_dom_size);
    for (unsigned int i = 0; i < V_dom_size; ++i){
      len_elong[i] =  (len_max_elong[i] <= trials_elong) ? 
                      regspace<ivec>(1, len_max_elong[i]): 
                      regspace<ivec>(1, (len_max_elong[i] - 1) / (double)(trials_elong - 1), len_max_elong[i]); // to be checked
    }
    // vector of probabilities for the quantile function
    vec prob(1,fill::value(0.25));
    // compute the quantile of the distance matrix
    vec quantile = arma::quantile(vectorise(D), prob);
    // keep will me a matrix whose value i,j will be D(i,j) < quantile(0)
    umat keep = D < quantile(0);
    // col-wise sum of the matrix keep
    const uvec & col_sum_keep = sum(keep, 0);
    // vector of bool = true iff col_sum_keep[i]==0
    const uvec & col_sum_keep_zero = (col_sum_keep==0);
    // empty_k stores the indexes of the col that have col_sum_keep = 0 
    const uvec & empty_k = find(col_sum_keep_zero);
    
    for (auto k : empty_k){
      const unsigned int min_col_k_D = index_min(D.col(k));
      keep(min_col_k_D, k) = true;
    } 
    // create a list with the cols of keep
    List Keep_k(keep.n_cols);
    unsigned int i = 0;
    keep.each_col([&Keep_k,&i](const uvec& v){Keep_k[i++] = v;});
    
    List V_new_temp(V_dom_size);
    List V_dom_temp(V_dom_size);
    List S_k_temp(V_dom_size);
    Function elongation_rcpp("elongation_rcpp"); // to be passed differently
    Function domain(".domain");
    Function compute_motif(".compute_motif");
    for (unsigned int i = 0; i < V_dom_size; ++i){
      const List & result = elongation_rcpp(V_new[i],
                                            V_dom[i],
                                            S_k[i],
                                            P_k[i],
                                            len_elong[i],
                                            Keep_k[i],
                                            c[i],
                                            domain,
                                            compute_motif,
                                            use0,
                                            use1,
                                            w, 
                                            alpha,
                                            max_gap,
                                            Y,
                                            m,
                                            deltaJk_elong); // discrepancy in c-type, to be checked
      V_new_temp[i] = result[0];
      V_dom_temp[i] = result[1];
      S_k_temp[i]   = result[2];
    }
    V_new = V_new_temp;
    V_dom = V_dom_temp;
    S_k = S_k_temp;
    mat S = mat(as<mat>(S_k));
    return S.reshape(S.n_elem/K,K);
}
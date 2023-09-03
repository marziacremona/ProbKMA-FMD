#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


// [[Rcpp::export]]
double compute_Jk_rcpp(const List & v,
                     const ivec & s_k,
                     const vec & p_k,
                     const List & Y,
                     double alpha,
                     const vec & w,
                     int m,
                     bool use0,
                     bool use1,
                     const Function & domain,
                     const Function & select_domain,
                     const Function& diss_d0_d1_L2,
                     Nullable<int> c_k = R_NilValue,
                     Nullable<LogicalVector> keep_k = R_NilValue){
  
  // domain of the centroid
  // uvec non gli piacciono a select_domain 
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  
  // length of the domain
  const unsigned int v_len = v_dom.size();
  
  // select the part of the domain of the centroid
  List v_new = select_domain(v, v_dom, use0, use1);

  const List & first_y = Y[0]; // Y is list of list, first_y is the first list, first_y[0] is the first list of first_y 
  const mat & first_y0 = as<mat>(first_y[0]);
  
  // dimensionality of the curves
  const unsigned int d = first_y0.n_cols;
  
  // length of the curves 
  const int y_len = first_y0.n_rows;
  
  const unsigned int Y_size = Y.size();
  
  // curves shifted as s_k says 
  List Y_inters_k(Y_size);
  
  for (unsigned int i = 0; i < Y_size; ++i){
    List y_inters_k = List::create(Named("y0") = R_NilValue,
                                   Named("y1") = R_NilValue);
    int s_k_i = s_k[i];
    
    ivec index = regspace<ivec>(1, v_len - std::max(0, 1-s_k_i))
                 +std::max(1,s_k_i)-1;

    unsigned int index_size = index.size();
    
    List y_i = Y[i];
    
    mat new_y0(index_size + std::max(0, 1-s_k_i), d);
    
    if (use0){
      mat temp_y0_i = as<mat>(y_i[0]);
      new_y0.fill(datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j]  <= y_len){
          uword index_row = std::max(0, 1-s_k_i) + j;
          new_y0.row(index_row) =  temp_y0_i.row(index[j] - 1);
          }
      }
      y_inters_k["y0"] = new_y0;
    }
    
    mat new_y1(index_size + std::max(0, 1-s_k_i), d);
    
    if (use1){
      mat temp_y1_i = as<mat>(y_i[1]);
      new_y1.fill(datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j] <= y_len){
          uword index_row = std::max(0, 1-s_k_i) + j;
          new_y1.row(index_row) =  temp_y1_i.row(index[j] - 1);
          }
      }
      y_inters_k["y1"] = new_y1;
    }
    Y_inters_k[i] = as<List>(select_domain(y_inters_k,v_dom,use0,use1));
  }
  
  // @TODO: check other implementations for lines 83 to 96
  // nullable objects of Rcpp have to be converted to usual objects (necessary as in the prof function)
  // @TODO: check output of this part of the code because in the test keep_k.isNotNull() && c_k.isNotNull() == FALSE
  if(keep_k.isNotNull() && c_k.isNotNull()){
    
    NumericVector supp_inters_length;
    
    LogicalVector keep_k_notnull = as<LogicalVector>(keep_k);
    
    for (uword i = 0; i < Y_size; ++i){
      if (keep_k_notnull[i]){
        LogicalVector domain_y_inters_k  = as<LogicalVector>(domain(Y_inters_k[i],use0));
        supp_inters_length.push_back(sum(domain_y_inters_k));
      }
    }
    
    int c_k_notnull = as<int>(c_k);
    
    LogicalVector check_lengths = supp_inters_length < c_k_notnull;
    
    if (is_true(any(check_lengths))){
      return NA_REAL;
    }
  }
  
  vec dist(Y_size);
  for (uword i = 0; i < Y_size; ++i){
    dist(i) = as<double>(diss_d0_d1_L2(Y_inters_k[i], v_new, w, alpha));
  }
  
  vec result = dist % pow(p_k,m);
  
  return accu(result.elem(find_finite(result))); //@TODO: improve this because così però escludo non solo i Nan ma anche gli infiniti
}




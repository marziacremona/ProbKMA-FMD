// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <numeric>
#include <ranges>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


IntegerVector myseq(int first, int last);

// [[Rcpp::export]]
NumericVector find_diss_rcpp(const List &y,const List &v,  
                            unsigned int c_k,
                            const NumericVector & w, 
                            double alpha, 
                            unsigned int d,bool use0,bool use1,
                            const Function & domain, 
                            const Function & select_domain,
                            const Function & diss_d0_d1_L2)
{
  // Rcout << "ciao sono dentro " << std::endl;
  
  // Convert domain and select_domain
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
  int v_len = v_dom.size();
  int y_len = as<mat>(y[0]).n_rows;
  IntegerVector s_rep = myseq(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
  std::size_t s_rep_size = s_rep.size();
  List y_rep(s_rep_size);
  
  // Convert y[0] and y[1] to mat objects
  const int index_size = v_len;
  mat temp_y0, temp_y1;
  if (use0) {
    temp_y0 = as<mat>(y[0]);
  }
  if (use1) {
    temp_y1 = as<mat>(y[1]);
  }
  auto index_range = std::views::iota(0,index_size);
  
  for (unsigned int i = 0; i < s_rep_size; ++i) {
    IntegerVector index = s_rep[i] - 1 + seq_len(v_len);
    List y_rep_i = List::create(Named("y0") = R_NilValue, Named("y1") = R_NilValue);
    auto j_true = index_range
      | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));});
    if (use0) {
      mat new_y0(index_size, d);
      new_y0.fill(datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y0,&temp_y0,&index](int j){new_y0.row(j) = temp_y0.row(index[j] - 1);});
      y_rep_i["y0"] = new_y0;
    }
    if (use1) {
      mat new_y1(index_size, d);
      new_y1.fill(datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y1,&temp_y1,&index](int j){new_y1.row(j) = temp_y1.row(index[j] - 1);});
      y_rep_i["y1"] = new_y1;
    }
    y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
    y_rep[i] = y_rep_i;
  }
  const std::size_t y_rep_size = y_rep.size();
  IntegerVector length_inter(y_rep_size);
  
  mat temp_y;
  int i = 0;
  for (const Rcpp::List& y_rep_i : y_rep) {
    if (use0) 
      temp_y = as<mat>(y_rep_i["y0"]);
    else
      temp_y = as<mat>(y_rep_i["y1"]);
    uvec non_na_indices = find_nan(temp_y.col(0)); 
    length_inter[i] = temp_y.col(0).n_elem - non_na_indices.n_elem;
    ++i;
  }
  
  LogicalVector valid = length_inter >= c_k;
  if (sum(valid) == 0) {        
    valid[length_inter == max(length_inter)] = true;  // valid[which_max(length_inter)] = true;
  }
  
  s_rep = s_rep[valid];
  y_rep = y_rep[valid];
  
  const unsigned int y_rep_size_new = y_rep.size();
    
  //NumericVector d_rep(y_rep_size_new);
  
  double min_d = std::numeric_limits<double>::max();
  int min_s = 0;
  
  for (int i = 0; i < y_rep_size_new; i++) {
    double dist = as<double>(diss_d0_d1_L2(y_rep[i], v_new, w, alpha));
    //d_rep[i] = dist;
    if (dist < min_d){
      min_d = dist;
      min_s = s_rep[i];
    }
  }
  return NumericVector::create(min_s, min_d); 
}
  
IntegerVector myseq(int first, int last) {
    IntegerVector y(abs(last - first) + 1);
    if (first < last) 
      std::iota(y.begin(), y.end(), first);
    else {
      std::iota(y.begin(), y.end(), last);
      std::reverse(y.begin(), y.end());
    }
    return y;
  }


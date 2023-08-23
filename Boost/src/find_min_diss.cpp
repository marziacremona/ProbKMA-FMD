// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <iostream>
#include <fstream>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

IntegerVector myseq(int first, int last);

// [[Rcpp::export]]
NumericVector find_diss(const List &y,const List &v,  
                        const NumericVector & w, double alpha, unsigned int c_k,
                        unsigned int d, bool use0, 
                        bool use1, const Function & domain, 
                        const Function & select_domain,
                        const Function & diss_d0_d1_L2)
{
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  Rcpp::List v_new = select_domain(as<List>(v), v_dom, use0, use1);
  int v_len = v_dom.size();
  int y_len = as<mat>(y[0]).n_rows; 
  IntegerVector s_rep = myseq(1 - (v_len - c_k),y_len - v_len + 1 + (v_len - c_k));

  List y_rep(s_rep.size());
  
  for (unsigned int i = 0; i < s_rep.size(); i++) {
    Rcpp::checkUserInterrupt();
    IntegerVector index = s_rep[i] - 1 + myseq(1,v_len);
    List y_rep_i = List::create(Named("y0") = R_NilValue, Named("y1") = R_NilValue);
    if (use0) {
      mat temp_y0 = as<mat>(y[0]); 
      mat new_y0(index.size(), d); //  sum(index <= 0) + sum((index > 0) & (index <= y_len)) + sum(index > y_len) == index.size()??
      new_y0.fill(datum::nan); 
      for (unsigned int j = 0; j < index.size(); j++) {
        if ((index[j] > 0) && (index[j] <= y_len)) {
          new_y0.row(j) = temp_y0.row(index[j] - 1); // try also with new_y0.row(sum(index <= 0) + j - sum(index <= index[j])) 
        }
      }
      y_rep_i["y0"] = new_y0;
    }
    if (use1) {
      mat temp_y1 = as<mat>(y[1]);
      mat new_y1(index.size(), d); // sum(index <= 0) + sum((index > 0) & (index <= y_len)) + sum(index > y_len)
      new_y1.fill(datum::nan);
      for (int j = 0; j < index.size(); j++) {
        Rcpp::checkUserInterrupt();
        if ((index[j] > 0) && (index[j] <= y_len)) {
          new_y1.row(j) = temp_y1.row(index[j] - 1);  //try also with new_y0.row(sum(index <= 0) + j - sum(index <= index[j])) 
        }
      }
      y_rep_i["y1"] = new_y1;
    }
    y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
    y_rep[i] = y_rep_i;
    }
  IntegerVector length_inter(y_rep.size());
  for (int i = 0; i < y_rep.size(); i++) {
    Rcpp::checkUserInterrupt();
    List temp_y_rep_i = as<List>(y_rep[i]);
    mat temp_y;
    if (use0) {
      temp_y = as<mat>(temp_y_rep_i["y0"]);
    } else {
      temp_y = as<mat>(temp_y_rep_i["y1"]);
    }
    uvec non_na_indices = find_nan(temp_y.col(0)); 
    length_inter[i] = temp_y.col(0).n_elem - non_na_indices.n_elem;
  }
  LogicalVector valid = length_inter >= c_k;
  if (sum(valid) == 0) {        
    valid[length_inter == max(length_inter)] = true;  // valid[which_max(length_inter)] = true;
  }
  
  s_rep = s_rep[valid]; 
  y_rep = y_rep[valid]; 
  
  NumericVector d_rep(y_rep.size());
  
  for (int i = 0; i < y_rep.size(); i++) {
    Rcpp::checkUserInterrupt();
    double dist = as<double>(diss_d0_d1_L2(as<List>(y_rep[i]), v_new, w, alpha)); 
    d_rep[i] = dist;
  }
  return NumericVector::create(s_rep[which_min(d_rep)], min(d_rep)); 
}

// more efficient version 
// [[Rcpp::export]]
NumericVector find_diss_efficient(const List &y,const List &v,
                        const NumericVector & w, double alpha, unsigned int c_k,
                        unsigned int d, bool use0, 
                        bool use1, const Function & domain, 
                        const Function & select_domain,
                        const Function & diss_d0_d1_L2)
{
  // Convert domain and select_domain
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
  int v_len = v_dom.size();
  int y_len = as<mat>(y[0]).n_rows;
  IntegerVector s_rep = myseq(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
  
  List y_rep(s_rep.size());
  
  // Convert y[0] and y[1] to mat objects
  mat temp_y0, temp_y1;
  if (use0) {
    temp_y0 = as<mat>(y[0]);
  }
  if (use1) {
    temp_y1 = as<mat>(y[1]);
  }
  
  for (unsigned int i = 0; i < s_rep.size(); i++) {
    IntegerVector index = s_rep[i] - 1 + seq_len(v_len);
    List y_rep_i = List::create(Named("y0") = R_NilValue, Named("y1") = R_NilValue);
    if (use0) {
      mat new_y0(index.size(), d);
      new_y0.fill(datum::nan);
      for (unsigned int j = 0; j < index.size(); j++) {
        if ((index[j] > 0) && (index[j] <= y_len)) {
          new_y0.row(j) = temp_y0.row(index[j] - 1); //TO BE CHECKED 
        }
      }
      y_rep_i["y0"] = new_y0;
    }
    if (use1) {
      mat new_y1(index.size(), d);
      new_y1.fill(datum::nan);
      for (int j = 0; j < index.size(); j++) {
        if ((index[j] > 0) && (index[j] <= y_len)) {
          new_y1.row(j) = temp_y1.row(index[j] - 1); // TO BE CHECKED 
        }
      }
      y_rep_i["y1"] = new_y1;
    }
    y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
    y_rep[i] = y_rep_i;
  }
  
  IntegerVector length_inter(y_rep.size());
  
  for (int i = 0; i < y_rep.size(); i++) {
    List temp_y_rep_i = as<List>(y_rep[i]);
    mat temp_y;
    if (use0) {
      temp_y = as<mat>(temp_y_rep_i["y0"]);
    } else {
      temp_y = as<mat>(temp_y_rep_i["y1"]);
    }
    uvec non_na_indices = find_finite(temp_y.col(0));
    length_inter[i] = non_na_indices.n_elem;
  }
  
  LogicalVector valid = length_inter >= c_k;
  if (sum(valid) == 0) {        
    valid[length_inter == max(length_inter)] = true;  // valid[which_max(length_inter)] = true;
  }
  
  s_rep = s_rep[valid];
  y_rep = y_rep[valid];
  
  NumericVector d_rep(y_rep.size());
  
  for (int i = 0; i < y_rep.size(); i++) {
    double dist = as<double>(diss_d0_d1_L2(as<List>(y_rep[i]), v_new, w, alpha));
    d_rep[i] = dist;
  }
  
  return NumericVector::create(s_rep[which_min(d_rep)], min(d_rep)); 
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


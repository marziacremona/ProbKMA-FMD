#ifndef COMPUTE_JK_HPP
#define COMPUTE_JK_HPP
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>
double compute_Jk_rcpp(const Rcpp::List & v,
                       const arma::ivec & s_k,
                       const arma::vec & p_k,
                       const Rcpp::List & Y,
                       double alpha,
                       const arma::vec & w,
                       int m,
                       bool use0,
                       bool use1,
                       const Rcpp::Function & domain,
                       const Rcpp::Function & select_domain,
                       const Rcpp::Function& diss_d0_d1_L2,
                       Rcpp::Nullable<int> c_k = R_NilValue,
                       Rcpp::Nullable<Rcpp::LogicalVector> keep_k = R_NilValue);




#endif //COMPUTE_JK_HPP
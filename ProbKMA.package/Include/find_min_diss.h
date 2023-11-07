#ifndef FIND_MIN_DISS_H
#define FIND_MIN_DISS_H
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>
#include <limits>
#include "utilities.h"


arma::vec find_diss(const Rcpp::List &y,const Rcpp::List &v,  
                    const arma::vec & w, 
                    double alpha, unsigned int c_k,
                    unsigned int d,bool use0,bool use1,
                    const Rcpp::Function & domain, 
                    const Rcpp::Function & select_domain,
                    const Rcpp::Function & diss_d0_d1_L2);

arma::vec find_diss_aligned_rcpp(const Rcpp::List &y,
                                 const Rcpp::List &v,  
                                 const arma::vec & w, 
                                 double alpha,
                                 bool aligned,
                                 unsigned int d,
                                 bool use0,
                                 bool use1,
                                 const Rcpp::Function & domain,
                                 const Rcpp::Function & select_domain,
                                 const Rcpp::Function & diss_d0_d1_L2);


Rcpp::List find_shift_warp_min(const Rcpp::List & Y, 
                               const Rcpp::List & V_new,
                               const arma::vec & w,
                               const arma::ivec & c_k,
                               unsigned int K,
                               unsigned int d,
                               double max_gap,
                               double alpha,
                               bool use0,
                               bool use1,
                               const Rcpp::Function & domain,
                               const Rcpp::Function & select_domain,
                               const Rcpp::Function & diss_d0_d1_L2);


#endif // FIND_MIN_DISS_H
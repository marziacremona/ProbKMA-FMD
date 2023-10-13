#ifndef FIND_MIN_DISS_H
#define FIND_MIN_DISS_H
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>
#include <limits>
#include "utilities.h"

Rcpp::NumericVector find_diss(const Rcpp::List &y,const Rcpp::List &v,  
                              const arma::vec & w, 
                              double alpha, unsigned int c_k,
                              unsigned int d,bool use0,bool use1,
                              const Rcpp::Function & domain, 
                              const Rcpp::Function & select_domain,
                              const Rcpp::Function & diss_d0_d1_L2);

#endif // FIND_MIN_DISS_H
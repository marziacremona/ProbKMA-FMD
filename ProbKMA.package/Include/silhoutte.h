#ifndef SILHOUTTE_H
#define SILHOUTTE_H

#include "RcppArmadillo.h"
#include "find_min_diss.h"
#include <limits>
#include <ranges>
#include <algorithm>
#include "utilities.h"

Rcpp::NumericVector find_diss_aligned_rcpp(const Rcpp::List &y,
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


Rcpp::List probKMA_silhouette_rcpp(const Rcpp::List & probKMA_results,
                                   const Rcpp::Function & domain,
                                   const Rcpp::Function & select_domain,
                                   const Rcpp::Function & diss_d0_d1_L2,
                                   bool align = false);


#endif // SILHOUTTE_H



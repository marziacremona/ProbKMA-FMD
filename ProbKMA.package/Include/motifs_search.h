#ifndef MOTIFS_SEARCH_H
#define MOTIFS_SEARCH_H

#include "RcppArmadillo.h"
#include "utilities.h"
#include "find_occurrences.h"
#include <vector>
#include <ranges>
#include <algorithm>
#include <list>

Rcpp::List motifs_search_cpp(const Rcpp::List& Y, // list of list of matrices
                             const Rcpp::List& V,  // list of list of matrices 
                             const Rcpp::List& V0_clean,
                             const Rcpp::List& V1_clean,
                             const Rcpp::List& V_dom, // list of LogicalVector
                             const arma::vec& V_length, // vector of V_dom_i's length 
                             const arma::mat& P_clean,
                             const arma::mat& D_clean,
                             arma::uvec V_hclust, //vector of indices related to clusters 
                             const double alpha, 
                             const bool use0, const bool use1,
                             const arma::vec& w, 
                             const arma::vec& c,
                             const double max_gap,  
                             const double d, // n_col of Y[[0]]
                             const double N, // n_row of D
                             const double K, // n_col of D
                             const double R_all,
                             const arma::vec& R_m,
                             bool use_real_occurrences,
                             double length_diff,
                             Rcpp::Function diss_d0_d1_L2, 
                             Rcpp::Function domain,
                             Rcpp::Function select_domain);


#endif  // MOTIFS_SEARCH_H
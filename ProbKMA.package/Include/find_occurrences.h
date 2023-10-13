#ifndef FIND_OCCURRENCES_H
#define FIND_OCCURRENCES_H
#include "RcppArmadillo.h"
#include <vector>
#include <ranges>
#include <algorithm>
#include <list>

arma::mat find_occurrences_cpp(const Rcpp::List& v,
                               const Rcpp::List& Y,
                               const double R,
                               const double alpha,
                               const arma::vec& w,
                               const double c_k,
                               const bool use0,
                               const bool use1,
                               Rcpp::Function diss_d0_d1_L2,
                               Rcpp::Function domain,
                               Rcpp::Function select_domain);

# endif // FIND_OCCURRENCES_H
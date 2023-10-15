#ifndef SILHOUTTE_H
#define SILHOUTTE_H

#include "RcppArmadillo.h"
#include "find_min_diss.h"
#include <limits>
#include <ranges>
#include <algorithm>
#include "utilities.h"



Rcpp::List probKMA_silhouette_rcpp(const Rcpp::List & probKMA_results,
                                   const Rcpp::Function & domain,
                                   const Rcpp::Function & select_domain,
                                   const Rcpp::Function & diss_d0_d1_L2,
                                   bool align = false);


#endif // SILHOUTTE_H



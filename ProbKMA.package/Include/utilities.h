#ifndef UTILITIES_HPP
#define UTILITIES_HPP
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>
namespace util
{
  Rcpp::List repeat_elements(const arma::imat& A,const arma::ivec & times);
}



#endif //UTILITIES_HPP
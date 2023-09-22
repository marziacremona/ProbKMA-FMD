#include "RcppArmadillo.h"
using namespace Rcpp;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

arma::uvec domain(const std::pair<arma::mat,arma::mat> &v,
                  bool use0){
  if (use0){
    arma::uvec result(v.first.n_rows,arma::fill::zeros);
    for (arma::uword i=0; i < v.first.n_rows; ++i){
      const arma::uvec & finite_row = find_finite(v.first.row(i));
      if(finite_row.n_elem)
        result(i) = 1;
    }
    return find(result == 1);
  } else {
    arma::uvec result(v.second.n_rows,arma::fill::zeros);
    for (arma::uword i=0; i < v.second.n_rows; ++i){
      const arma::uvec & finite_row = find_finite(v.second.row(i));
      if(finite_row.n_elem) 
        result(i) = 1;
    }
    return find(result == 1);
  }
}

void select_domain(std::pair<arma::mat,arma::mat> &v,
                   const arma::uvec &v_dom,
                   bool use0,
                   bool use1){
  if(use0)
    v.first = v.first.rows(v_dom);
  if(use1)
    v.second = v.second.rows(v_dom);
}

// [[Rcpp::export]]
List test_domain_select_domain(List & v,
                               bool use0,
                               bool use1){
  arma::mat v0 = as<arma::mat>(v[0]);
  arma::mat v1 = as<arma::mat>(v[1]);
  std::pair<arma::mat,arma::mat> v_pair(std::make_pair<arma::mat,arma::mat>(std::move(v0),std::move(v1)));
  arma::uvec my_domain = domain(v_pair,use0);
  select_domain(v_pair,my_domain,use0,use1);
  return List::create(my_domain);
}


#ifndef ELONGATE_MOTIFS_HPP
#define ELONGATE_MOTIFS_HPP
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>
#include "utilities.h"
#include "compute_jk.h"

Rcpp::List elongation_rcpp(const Rcpp::List & v_new_k, 
                           const arma::uvec& v_dom_k,  
                           const arma::ivec& s_k, 
                           const arma::vec & p_k, 
                           const arma::ivec& len_elong_k, 
                           const arma::uvec& keep_k,  
                           double c, 
                           const Rcpp::Function& domain,
                           const Rcpp::Function& compute_motif,
                           bool use0, bool use1,
                           const arma::vec& w, 
                           double alpha, double max_gap,  
                           const Rcpp::List& Y, int m, double deltaJk_elong);

void elongate_motifs(Rcpp::List & V_new,
                      Rcpp::List & V_dom,
                      Rcpp::List & S_k,
                      const Rcpp::List & P_k,
                      const Rcpp::List & Y,
                      const arma::vec & w, //in this case is a number
                      int m, //in this case is a number
                      bool use0,
                      bool use1,
                      double alpha,
                      const arma::ivec & c,
                      const arma::ivec & c_max, // is a vector to be understood
                      double max_elong, 
                      double deltaJk_elong,
                      int trials_elong,
                      const arma::mat & D,
                      unsigned int K,
                      double max_gap);


#endif //ELONGATE_MOTIFS_HPP
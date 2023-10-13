#include "find_min_diss.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

Rcpp::NumericVector find_diss(const Rcpp::List &y,const Rcpp::List &v,  
                        const arma::vec & w, 
                        double alpha, unsigned int c_k,
                        unsigned int d,bool use0,bool use1,
                        const Rcpp::Function & domain, 
                        const Rcpp::Function & select_domain,
                        const Rcpp::Function & diss_d0_d1_L2)
{
  // Convert domain and select_domain
  Rcpp::LogicalVector v_dom = Rcpp::as<Rcpp::LogicalVector>(domain(v,use0));
  Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
  int v_len = v_dom.size();
  int y_len = Rcpp::as<arma::mat>(y[0]).n_rows;
  Rcpp::IntegerVector s_rep = util::myseq(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
  std::size_t s_rep_size = s_rep.size();
  Rcpp::List y_rep(s_rep_size);
  
  // Convert y[0] and y[1] to mat objects
  const int index_size = v_len;
  arma::mat temp_y0, temp_y1;
  if (use0) {
    temp_y0 = Rcpp::as<arma::mat>(y[0]);
  }
  if (use1) {
    temp_y1 = Rcpp::as<arma::mat>(y[1]);
  }
  auto index_range = std::views::iota(0,index_size);
  
  for (unsigned int i = 0; i < s_rep_size; ++i) {
    Rcpp::IntegerVector index = s_rep[i] - 1 + Rcpp::seq_len(v_len);
    Rcpp::List y_rep_i = Rcpp::List::create(Rcpp::Named("y0") = R_NilValue, Rcpp::Named("y1") = R_NilValue);
    auto j_true = index_range
    | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));});
    if (use0) {
      arma::mat new_y0(index_size, d);
      new_y0.fill(arma::datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y0,&temp_y0,&index](int j){new_y0.row(j) = temp_y0.row(index[j] - 1);});
      y_rep_i["y0"] = new_y0;
    }
    if (use1) {
      arma::mat new_y1(index_size, d);
      new_y1.fill(arma::datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y1,&temp_y1,&index](int j){new_y1.row(j) = temp_y1.row(index[j] - 1);});
      y_rep_i["y1"] = new_y1;
    }
    y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
    y_rep[i] = y_rep_i;
  }
  const unsigned int y_rep_size = y_rep.size();
  Rcpp::IntegerVector length_inter(y_rep_size);
  
  arma::mat temp_y;
  int i = 0;
  for (const Rcpp::List& y_rep_i : y_rep) { // @TODO: risolvere warning
    if (use0) 
      temp_y = Rcpp::as<arma::mat>(y_rep_i["y0"]);
    else
      temp_y = Rcpp::as<arma::mat>(y_rep_i["y1"]);
    arma::uvec non_na_indices = find_nan(temp_y.col(0)); 
    length_inter[i] = temp_y.col(0).n_elem - non_na_indices.n_elem;
    ++i;
  }
  
  Rcpp::LogicalVector valid = length_inter >= c_k;
  if (sum(valid) == 0) {        
    valid[length_inter == max(length_inter)] = true; 
  }
  
  s_rep = s_rep[valid];
  y_rep = y_rep[valid];
  
  const unsigned int y_rep_size_valid = y_rep.size(); 
  
  double min_d = std::numeric_limits<double>::max();
  int min_s = 0;
  
  for (unsigned int i = 0; i < y_rep_size_valid; i++) {
    double dist = Rcpp::as<double>(diss_d0_d1_L2(y_rep[i], v_new, w, alpha));
    if (dist < min_d){
      min_d = dist;
      min_s = s_rep[i];
    }
  }
  return Rcpp::NumericVector::create(min_s, min_d); 
}
  



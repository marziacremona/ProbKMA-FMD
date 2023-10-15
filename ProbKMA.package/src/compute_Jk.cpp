#include "compute_jk.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

//[[Rcpp::export(.compute_Jk_rcpp)]]
double compute_Jk_rcpp(const Rcpp::List & v,
                       const arma::ivec & s_k,
                       const arma::vec & p_k,
                       const Rcpp::List & Y,
                       double alpha,
                       const arma::vec & w,
                       int m,
                       bool use0,
                       bool use1,
                       const Rcpp::Function & domain,
                       const Rcpp::Function & select_domain,
                       const Rcpp::Function& diss_d0_d1_L2,
                       Rcpp::Nullable<int> c_k,
                       Rcpp::Nullable<Rcpp::LogicalVector> keep_k){
  
  // domain of the centroid
  // uvec non gli piacciono a select_domain 
  Rcpp::LogicalVector v_dom = Rcpp::as<Rcpp::LogicalVector>(domain(v,use0));
  
  // length of the domain
  const unsigned int v_len = v_dom.size();
  
  // select the part of the domain of the centroid
  Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
  
  const Rcpp::List & first_y = Y[0]; // Y is list of list, first_y is the first list, first_y[0] is the first list of first_y 
  const arma::mat & first_y0 = use0 ? as<arma::mat>(first_y[0]) : as<arma::mat>(first_y[1]);
  
  // dimensionality of the curves
  const unsigned int d = first_y0.n_cols;
  
  const unsigned int Y_size = Y.size();
  
  // curves shifted as s_k says 
  Rcpp::List Y_inters_k(Y_size);
  
  for (unsigned int i = 0; i < Y_size; ++i){
   
    Rcpp::List y_inters_k = Rcpp::List::create(Rcpp::Named("y0") = R_NilValue,
                                               Rcpp::Named("y1") = R_NilValue);
    int s_k_i = s_k[i];
    
    arma::ivec index = arma::regspace<arma::ivec>(1, v_len - std::max(0, 1-s_k_i))
      +std::max(1,s_k_i)-1;
    
    unsigned int index_size = index.size();
    
    Rcpp::List y_i = Y[i];
    
    // length of the curve 
    const int y_len = as<arma::mat>(y_i[0]).n_rows;
    
    arma::mat new_y0(index_size + std::max(0, 1-s_k_i), d);
    
    if (use0){
      arma::mat temp_y0_i = Rcpp::as<arma::mat>(y_i[0]);
      new_y0.fill(arma::datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j]  <= y_len){
          arma::uword index_row = std::max(0, 1-s_k_i) + j;
          new_y0.row(index_row) =  temp_y0_i.row(index[j] - 1);
        }
      }
      y_inters_k["y0"] = new_y0;
    }
    
    arma::mat new_y1(index_size + std::max(0, 1-s_k_i), d);
    
    if (use1){
      arma::mat temp_y1_i = Rcpp::as<arma::mat>(y_i[1]);
      new_y1.fill(arma::datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j] <= y_len){
          arma::uword index_row = std::max(0, 1-s_k_i) + j;
          new_y1.row(index_row) =  temp_y1_i.row(index[j] - 1);
        }
      }
      y_inters_k["y1"] = new_y1;
    }
    Y_inters_k[i] = Rcpp::as<Rcpp::List>(select_domain(y_inters_k,v_dom,use0,use1));
  }
  
  // @TODO: check other implementations for lines 83 to 96
  // nullable objects of Rcpp have to be converted to usual objects (necessary as in the prof function)
  // @TODO: check output of this part of the code because in the test keep_k.isNotNull() && c_k.isNotNull() == FALSE
  if(keep_k.isNotNull() && c_k.isNotNull()){
    
    Rcpp::NumericVector supp_inters_length;
    
    Rcpp::LogicalVector keep_k_notnull = Rcpp::as<Rcpp::LogicalVector>(keep_k);
    
    for (arma::uword i = 0; i < Y_size; ++i){
      if (keep_k_notnull[i]){
        Rcpp::LogicalVector domain_y_inters_k  = Rcpp::as<Rcpp::LogicalVector>(domain(Y_inters_k[i],use0));
        supp_inters_length.push_back(sum(domain_y_inters_k));
      }
    }
    
    int c_k_notnull = as<int>(c_k);
    
    Rcpp::LogicalVector check_lengths = supp_inters_length < c_k_notnull;
    
    if (is_true(any(check_lengths))){
      return NA_REAL;
    }
  }
  
  arma::vec dist(Y_size);
  for (arma::uword i = 0; i < Y_size; ++i){
    dist(i) = Rcpp::as<double>(diss_d0_d1_L2(Y_inters_k[i], v_new, w, alpha));
  }
  
  arma::vec result = dist % pow(p_k,m);
  
  return arma::accu(result.elem(arma::find_finite(result))); //@TODO: improve this because cos? per? escludo non solo i Nan ma anche gli infiniti
}



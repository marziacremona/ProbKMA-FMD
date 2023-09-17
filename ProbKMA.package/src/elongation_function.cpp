// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <gperftools/profiler.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

List repeat_elements(const arma::imat& A,const arma::ivec & times) {
  arma::uword times_size = times.n_elem;
  List result((times_size*(times_size+1))/2 - 1);
  std::size_t i = 0;
  for(arma::uword j = 0;j < times_size;++j)
  {
    const arma::imat& B = repmat(A.col(j),1,times[j]);
    B.each_col([&result,&i](const arma::ivec& v){result[i++] = v;});
  }
  return result;
}



double compute_Jk_rcpp(const List & v,
                       const arma::ivec & s_k,
                       const arma::vec & p_k,
                       const List & Y,
                       double alpha,
                       const arma::vec & w,
                       int m,
                       bool use0,
                       bool use1,
                       Nullable<int> c_k = R_NilValue,
                       Nullable<LogicalVector> keep_k = R_NilValue){
  
  Function domain(".domain");
  Function select_domain(".select_domain");
  Function diss_d0_d1_L2(".diss_d0_d1_L2");
  
  // domain of the centroid
  // arma::uvec non gli piacciono a select_domain 
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  
  // length of the domain
  const unsigned int v_len = v_dom.size();
  
  // select the part of the domain of the centroid
  List v_new = select_domain(v, v_dom, use0, use1);
  
  const List& first_y = Y[0]; // Y is list of list, first_y is the first list, first_y[0] is the first list of first_y 
  const arma::mat& first_y0 = use0 ? as<arma::mat>(first_y[0]):as<arma::mat>(first_y[1]);
  
  // dimensionality of the curves
  const unsigned int d = first_y0.n_cols;
  
  // length of the curves 
  const int y_len = first_y0.n_rows;
  
  const unsigned int Y_size = Y.size();
  
  // curves shifted as s_k says 
  List Y_inters_k(Y_size);
  for (unsigned int i = 0; i < Y_size; ++i){
    List y_inters_k = List::create(Named("y0") = R_NilValue,
                                   Named("y1") = R_NilValue);
    const int& s_k_i = s_k[i];
    const arma::ivec& index = regspace<arma::ivec>(1, v_len - std::max(0, 1-s_k_i)) + std::max(1,s_k_i) - 1;
    int index_size = index.size();
    const List& y_i = Y[i];
    
    arma::mat new_y01(index_size + std::max(0, 1-s_k_i), d);
    if (use0){
      const arma::mat& temp_y0_i = as<arma::mat>(y_i[0]);
      new_y01.fill(datum::nan);
      
      auto filtered_j = std::views::iota(0,index_size)
        | std::views::filter([&y_len,index](int j){return index[j] <= y_len;});
      for(int j : filtered_j)
        new_y01.row(std::max(0, 1-s_k_i) + j) =  temp_y0_i.row(index[j] - 1);
      
      y_inters_k["y0"] = new_y01;
    }

    if (use1){
      const arma::mat& temp_y1_i = as<arma::mat>(y_i[1]);
      new_y01.fill(datum::nan);
      
      auto filtered_j = std::views::iota(0,index_size)
        | std::views::filter([&y_len,index](int j){return index[j] <= y_len;});
      for(int j : filtered_j)
        new_y01.row(std::max(0, 1-s_k_i) + j) =  temp_y1_i.row(index[j] - 1);
      
      y_inters_k["y1"] = new_y01;
    }
    Y_inters_k[i] = as<List>(select_domain(y_inters_k,v_dom,use0,use1));
  }
  
  // nullable objects of Rcpp have to be converted to usual objects (necessary as in the prof function)
  // @TODO: check output of this part of the code because in the test keep_k.isNotNull() && c_k.isNotNull() == FALSE
  if(keep_k.isNotNull() && c_k.isNotNull()){
    
    LogicalVector keep_k_notnull = as<LogicalVector>(keep_k);
    auto filtered_Y_inters = std::views::iota(0,Y.size()) 
        | std::views::filter([&keep_k_notnull](int j){return keep_k_notnull[j];})
        | std::views::transform([&Y_inters_k](int j){return Y_inters_k[j];});
    ProfilerStart("myprof.log");
    NumericVector supp_inters_length(std::ranges::distance(filtered_Y_inters));
    for (unsigned int i = 0; auto y_filtered:filtered_Y_inters){
      arma::uvec domain_y_inters_k  = as<arma::uvec>(domain(y_filtered,use0));
      supp_inters_length[i++] = sum(domain_y_inters_k);
    }
    ProfilerStop();
    int c_k_notnull = as<int>(c_k);
    LogicalVector check_lengths = supp_inters_length < c_k_notnull;
    
    if (is_true(any(check_lengths))) return NA_REAL;
  }
  
  vec dist(Y_size);
  for (arma::uword i = 0; i < Y_size; ++i)
  {
    dist[i] = as<double>(diss_d0_d1_L2(Y_inters_k[i], v_new, w, alpha));
  }
  vec result = dist % pow(p_k,m);
  
  return sum(result.elem(find_finite(result))); //@TODO: improve this because cos� per� escludo non solo i Nan ma anche gli infiniti
}


// [[Rcpp::export]]
List elongation_rcpp(const List & v_new_k, 
                     const arma::uvec& v_dom_k,  
                     const arma::ivec& s_k, 
                     const vec & p_k, 
                     const arma::ivec& len_elong_k, 
                     const arma::uvec& keep_k,  
                     double c, 
                     const Function& domain,
                     const Function& compute_motif,
                     bool use0, bool use1,
                     const vec& w, 
                     double alpha, double max_gap,  
                     const List& Y, int m, double deltaJk_elong) 
{
  if(len_elong_k.empty()) {
    return List::create(Named("v_new") = v_new_k,
                        Named("v_dom") = as<LogicalVector>(wrap(v_dom_k)),
                        Named("s_k") = s_k);
  }
  
  // new vec with zero at the top
  arma::ivec len_elong_k_zero(len_elong_k.size() + 1, fill::zeros);
  std::copy(len_elong_k.begin(), len_elong_k.end(), len_elong_k_zero.begin() + 1);
  
  // create a matrix whose column_i contains the vector s_k - len_elong_k_zero[i]
  arma::uword len_elong_k_zero_size = len_elong_k_zero.size();
  arma::imat s_k_elong_left_right_temp(s_k.n_elem, len_elong_k_zero_size);
  
  for (arma::uword i=0; i < len_elong_k_zero_size;++i) {
    s_k_elong_left_right_temp.col(i) = s_k - len_elong_k_zero(i);
  }
  
  // create a sequence of integer from len_elong_k_zero.size() to 1
  arma::ivec reversedSequence = arma::regspace<arma::ivec>(len_elong_k_zero_size,-1,1);
  reversedSequence(0) -= 1;
  
  // repeat each col of s_k_elong_left_right a number of times specified by reversedSequence and return a list 
  List s_k_elong_left_right =  repeat_elements(s_k_elong_left_right_temp, reversedSequence);
  
  //  v_dom_elong_left_right will be a vector of LogicalVector containing all the elongated version of v_dom_k
  const std::size_t v_dom_k_len = v_dom_k.n_elem;
  std::size_t max_len_elong_k = len_elong_k.back();
  List v_dom_elong_left_right((max_len_elong_k+1)*(max_len_elong_k+2)/2); 
  
  for(std::size_t i = 0;const int len_elong_k_left:len_elong_k_zero)
  {
    arma::uvec temp(len_elong_k_left + v_dom_k_len + len_elong_k_zero[max_len_elong_k],fill::ones);
    std::copy(v_dom_k.begin(), v_dom_k.end(), temp.begin() + len_elong_k_left);
    for(unsigned int j = 0; j < max_len_elong_k+1;++j)
    {
      const arma::uvec& temp_subwiew = temp(regspace<arma::uvec>(0,1,len_elong_k_left + v_dom_k_len + len_elong_k_zero[j]-1)); 
      v_dom_elong_left_right[i++] = temp_subwiew;
    }
    --max_len_elong_k;
  }
  
  // remove the last element (if any) of v_dom_elong_left_right (reason why we have choose std::vector in place of List)
  if (v_dom_elong_left_right.size()) 
    v_dom_elong_left_right.erase(0);
  
  // create the list containing all the possible v_dom_k elongated using compute_motif
  const int v_elong_left_right_size = s_k_elong_left_right.size();
  List v_elong_left_right(v_elong_left_right_size); 
  for (int i = 0; i < v_elong_left_right_size; i++) {
    v_elong_left_right[i] = as<List>(compute_motif(as<LogicalVector>(wrap(v_dom_elong_left_right[i])), s_k_elong_left_right[i], p_k, Y, m, use0, use1));
  }
  
  // create a LogicalVector start_with_NA whose elements are true iff the correspodent elemet of v_elong_left_right has lenght >2
  LogicalVector start_with_NA(v_elong_left_right_size);
  auto check_length = [](const List & v_elong_left_right_i){return(v_elong_left_right_i.size()>2);};
  std::transform(v_elong_left_right.begin(),v_elong_left_right.end(),start_with_NA.begin(),check_length);
      
  // filter the domain, centroid and shifts that are not in NA positions
  auto not_NA_index = std::views::iota(0,v_elong_left_right_size) 
    | std::views::filter([&start_with_NA](int index_j){return(!start_with_NA[index_j]);});
  auto filtered_elong = not_NA_index 
                        | std::views::transform([&v_elong_left_right](int j){return v_elong_left_right[j];});
  auto filtered_s_k = not_NA_index 
                      | std::views::transform([&s_k_elong_left_right](int j){return s_k_elong_left_right[j];});
  
  v_elong_left_right = List(filtered_elong.begin(),filtered_elong.end());
  s_k_elong_left_right = List(filtered_s_k.begin(),filtered_s_k.end());
  
  // compute performance index before elongation
  double Jk_before = compute_Jk_rcpp(v_new_k, s_k, p_k, Y, alpha, w, m, use0, use1);
  
  // compute performance indexes for all possible elongations
  vec c_k_after(v_elong_left_right_size);
  vec Jk_after(v_elong_left_right_size);
  
  for (arma::uword i = 0; i < v_elong_left_right_size; i++) {
    const arma::uvec& domain_elong = as<arma::uvec>(domain(v_elong_left_right[i], use0));
    int c_i =  std::max(floor(domain_elong.n_elem*(1 - max_gap)),c); 
    c_k_after[i] = c_i; // avoid if, next line
    const List& v_i = v_elong_left_right[i];
    const arma::ivec& s_i = s_k_elong_left_right[i];
    Jk_after[i] = compute_Jk_rcpp(v_i, s_i, p_k, Y, alpha, w, m,use0 , use1,wrap(c_i),as<LogicalVector>(wrap(keep_k)));
  }

  // find the best elongation in terms of perf. index
  vec diff_perc = ((Jk_after-Jk_before)/Jk_before);
  arma::uword best_elong = index_min(diff_perc);

  // check that the min really exists
  bool elongate = false;
  if (best_elong < v_elong_left_right_size)
    elongate = diff_perc(best_elong) < deltaJk_elong;
  
  // evaluate if elongate or not
  if(elongate) {
    return List::create(Named("v_new") = v_elong_left_right[best_elong],
                        Named("v_dom") = as<LogicalVector>(wrap(v_dom_elong_left_right[best_elong])),
                        Named("s_k")   = s_k_elong_left_right[best_elong]);
  } else {
    return List::create(Named("v_new") = v_new_k,
                        Named("v_dom") = as<LogicalVector>(wrap(v_dom_k)),
                        Named("s_k")   = s_k);
  }
}


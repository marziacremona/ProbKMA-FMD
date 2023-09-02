// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

List repeat_elements(const imat& A,const ivec & times) {
  uword times_size = times.n_elem;
  List result((times_size*(times_size+1))/2 - 1);
  std::size_t i = 0;
  for(uword j = 0;j < times_size;++j)
  {
    const imat& B = repmat(A.col(j),1,times[j]);
    B.each_col([&result,&i](const ivec& v){result[i++] = v;});
  }
  return result;
}

double compute_Jk_rcpp(const List & v,
                       const ivec & s_k,
                       const vec & p_k,
                       const List & Y,
                       double alpha,
                       const vec & w,
                       int m,
                       bool use0,
                       bool use1,
                       Nullable<int> c_k = R_NilValue,
                       Nullable<LogicalVector> keep_k = R_NilValue){
  
  Function domain(".domain");
  Function select_domain(".select_domain");
  Function diss_d0_d1_L2(".diss_d0_d1_L2");
  
  // domain of the centroid
  // uvec non gli piacciono a select_domain 
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  
  // length of the domain
  const unsigned int v_len = v_dom.size();
  
  // select the part of the domain of the centroid
  List v_new = select_domain(v, v_dom, use0, use1);
  
  const List & first_y = Y[0]; // Y is list of list, first_y is the first list, first_y[0] is the first list of first_y 
  const mat & first_y0 = as<mat>(first_y[0]);
  
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
    int s_k_i = s_k[i];
    
    const ivec & index = regspace<ivec>(1, v_len - std::max(0, 1-s_k_i)) + std::max(1,s_k_i) - 1;
    
    unsigned int index_size = index.size();
    
    const List & y_i = Y[i];
    
    if (use0){
      const mat & temp_y0_i = as<mat>(y_i[0]);
      mat new_y0(index_size + std::max(0, 1-s_k_i), d);
      new_y0.fill(datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j]  <= y_len){
          new_y0.row(std::max(0, 1-s_k_i) + j) =  temp_y0_i.row(index[j] - 1);
        }
      }
      y_inters_k["y0"] = new_y0;
    }

    if (use1){
      const mat & temp_y1_i = as<mat>(y_i[1]);
      mat new_y1(index_size + std::max(0, 1-s_k_i), d);
      new_y1.fill(datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j] <= y_len){
          new_y1.row(std::max(0, 1-s_k_i) + j) =  temp_y1_i.row(index[j] - 1);
        }
      }
      y_inters_k["y1"] = new_y1;
    }
    Y_inters_k[i] = as<List>(select_domain(y_inters_k,v_dom,use0,use1));
  }
  
  // @TODO: check other implementations for lines 83 to 96
  // nullable objects of Rcpp have to be converted to usual objects (necessary as in the prof function)
  // @TODO: check output of this part of the code because in the test keep_k.isNotNull() && c_k.isNotNull() == FALSE
  if(keep_k.isNotNull() && c_k.isNotNull()){
    
    NumericVector supp_inters_length;
    
    LogicalVector keep_k_notnull = as<LogicalVector>(keep_k);
    
    for (uword i = 0; i < Y_size; ++i){
      if (keep_k_notnull[i]){
        LogicalVector domain_y_inters_k  = as<LogicalVector>(domain(Y_inters_k[i],use0));
        supp_inters_length.push_back(sum(domain_y_inters_k));
      }
    }
    
    int c_k_notnull = as<int>(c_k);
    
    LogicalVector check_lengths = supp_inters_length < c_k_notnull;
    
    if (is_true(any(check_lengths))){
      return NA_REAL;
    }
  }
  
  vec dist(Y_size);
  for (uword i = 0; i < Y_size; ++i){
    dist(i) = as<double>(diss_d0_d1_L2(Y_inters_k[i], v_new, w, alpha));
  }
  
  vec result = dist % pow(p_k,m);
  
  return accu(result.elem(find_finite(result))); //@TODO: improve this because cos� per� escludo non solo i Nan ma anche gli infiniti
}




// [[Rcpp::export]]
List elongation_rcpp(const List & v_new_k, 
                     const uvec& v_dom_k,  
                     const ivec& s_k, 
                     const vec & p_k, 
                     const ivec& len_elong_k, 
                     const uvec& keep_k,  
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
  ivec len_elong_k_zero(len_elong_k.size() + 1, fill::zeros);
  std::copy(len_elong_k.begin(), len_elong_k.end(), len_elong_k_zero.begin() + 1);
  
  // create a matrix whose column_i contains the vector s_k - len_elong_k_zero[i]
  uword len_elong_k_zero_size = len_elong_k_zero.size();
  imat s_k_elong_left_right_temp(s_k.n_elem, len_elong_k_zero_size);
  
  for (uword i=0; i < len_elong_k_zero_size;++i) {
    s_k_elong_left_right_temp.col(i) = s_k - len_elong_k_zero(i);
  }
  
  // create a sequence of integer from len_elong_k_zero.size() to 1
  ivec reversedSequence = regspace<ivec>(len_elong_k_zero_size,-1,1);
  reversedSequence(0) -= 1;
  
  // repeat each col of s_k_elong_left_right a number of times specified by reversedSequence and return a list 
  List s_k_elong_left_right =  repeat_elements(s_k_elong_left_right_temp, reversedSequence);
  
  //  v_dom_elong_left_right will be a vector of LogicalVector containing all the elongated version of v_dom_k
  const std::size_t v_dom_k_len = v_dom_k.n_elem;
  std::size_t max_len_elong_k = len_elong_k.back();
  List v_dom_elong_left_right((max_len_elong_k+1)*(max_len_elong_k+2)/2); 
  std::size_t i = 0;
  for(const int len_elong_k_left:len_elong_k_zero)
  {
    uvec temp(len_elong_k_left + v_dom_k_len + len_elong_k_zero[max_len_elong_k],fill::ones);
    std::copy(v_dom_k.begin(), v_dom_k.end(), temp.begin() + len_elong_k_left);
    for(unsigned int j = 0; j < max_len_elong_k+1;++j)
    {
      const uvec& temp_subwiew = temp(regspace<uvec>(0,1,len_elong_k_left + v_dom_k_len + len_elong_k_zero[j]-1)); //v_dom_k_len = 0 invalido
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
  List new_v(std::distance(not_NA_index.begin(),not_NA_index.end()));
  List new_s(std::distance(not_NA_index.begin(),not_NA_index.end()));
  i = 0;
  std::for_each(not_NA_index.begin(),not_NA_index.end(),
                [&new_v,&new_s,&v_elong_left_right,
                 &s_k_elong_left_right,&i]
                 (int j){new_v[i] = as<List>(v_elong_left_right[j]);
                         new_s[i++] = s_k_elong_left_right[j];});
  v_elong_left_right = new_v;
  s_k_elong_left_right = new_s;
  
  // compute performance index before elongation
  double Jk_before = compute_Jk_rcpp(v_new_k, s_k, p_k, Y, alpha, w, m, use0, use1);
  
  // compute performance indexes for all possible elongations
  vec c_k_after(v_elong_left_right_size);
  vec Jk_after(v_elong_left_right_size);
  
  for (uword i = 0; i < v_elong_left_right_size; i++) {
    const uvec& domain_elong = as<uvec>(domain(v_elong_left_right[i], use0));
    int c_i =  std::max(floor(domain_elong.n_elem*(1 - max_gap)),c); 
    c_k_after(i) = c_i; // avoid if, next line
    const List& v_i = v_elong_left_right[i];
    const ivec& s_i = s_k_elong_left_right[i];
    Jk_after(i) = compute_Jk_rcpp(v_i, s_i, p_k, Y, alpha, w, m,use0 , use1,wrap(c_i),as<LogicalVector>(wrap(keep_k)));
  }

  // find the best elongation in terms of perf. index
  vec diff_perc = ((Jk_after-Jk_before)/Jk_before);
  uword best_elong = index_min(diff_perc);

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


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// [[Rcpp::export]]
List repeat_elements(const mat& A,const ivec & times) {
  List result;
  uword times_size = times.n_elem;
  for(uword j = 0;j< times_size;++j)
  {
    const mat& B = repmat(A.col(j),1,times[j]);
    B.each_col([&result](const vec& v){result.push_back(v);});
  }
  return result;
}


// [[Rcpp::export]]
List elongation_rcpp(const List & v_new_k, // ok 
                     const uvec& v_dom_k, // ok 
                     const vec& s_k, // ok 
                     const vec & p_k, // ok 
                     const vec& len_elong_k, // ok 
                     const uvec& keep_k, // ok 
                     double c, // ok 
                     const Function& compute_Jk,
                     const Function& domain,
                     const Function& compute_motif,
                     bool use0, bool use1,
                     const vec& w, 
                     double alpha, double max_gap, //ok 
                     const List& Y, int m, double deltaJk_elong) // ok 
{
  if(len_elong_k.empty()) {
    return List::create(Named("v_new") = v_new_k,
                        Named("v_dom") = v_dom_k,
                        Named("s_k") = s_k);
  }
  
  // new list with zero at the top
  List len_elong_k_zero(len_elong_k.begin(),len_elong_k.end()); // list because has push_front 
  len_elong_k_zero.push_front(0);

  // create a matrix whose column i contains the vector s_k - len_elong_k_zero[i]
  uword len_elong_k_zero_size = len_elong_k_zero.size();
  mat s_k_elong_left_right_temp(s_k.n_elem, len_elong_k_zero_size);
  //!! double len_elong_k_zero_i;
  for (uword i=0; i < len_elong_k_zero_size;++i) {
   //!! len_elong_k_zero_i = len_elong_k_zero[i]; // it is a way to force the conversion to double
    s_k_elong_left_right_temp.col(i) = s_k - static_cast<double>(len_elong_k_zero[i]);
  }
  // create a sequence of integer from len_elong_k_zero.size() to 1
  //!! ivec sequence = conv_to<Col<int>>::from(regspace(len_elong_k_zero_size, 1)); // @TODO : conv_to<Col<int>>::from to be understood
  //!! ivec reversedSequence = reverse(sequence);
  ivec reversedSequence = regspace<ivec>(len_elong_k_zero_size,-1,1);
  
  // repeat each col of s_k_elong_left_right a number of times specified by reversedSequence and return a list 
  List  s_k_elong_left_right =  repeat_elements(s_k_elong_left_right_temp, reversedSequence);
  
  // remove the last element (if any) of _k_elong_left_right
  if (s_k_elong_left_right.size() > 0) {
    s_k_elong_left_right.erase(s_k_elong_left_right.end()-1); // @TODO : maybe a different data-struct is more efficient
  }
  
  //  v_dom_elong_left_right will be a vector of LogicalVector containing all the elongated version of v_dom_k
  List v_dom_elong_left_right; // we should use uvec for consistency?!
  const int v_dom_k_len = v_dom_k.n_elem;
  const int max_len_elong_k = len_elong_k.max(); 
  for (const int len_elong_k_left: len_elong_k_zero) {  
    for (const int len_elong_k_right: len_elong_k_zero) {
      if (len_elong_k_right <= max_len_elong_k - len_elong_k_left) {
        LogicalVector temp(len_elong_k_left + v_dom_k_len + len_elong_k_right);
        std::fill_n(temp.begin(), len_elong_k_left, true);
        std::copy(v_dom_k.begin(), v_dom_k.end(), temp.begin() + len_elong_k_left);
        std::fill_n(temp.end() -len_elong_k_right, len_elong_k_right, true);
        v_dom_elong_left_right.push_back(temp);
      }
    }
  }
  
  // remove the last element (if any) of v_dom_elong_left_right (reason why we have choose std::vector in place of List)
  if (v_dom_elong_left_right.size()) 
    v_dom_elong_left_right.erase(v_dom_elong_left_right.end()-1);
 
  // create the list containing all the possible v_dom_k elongated using compute_motif
  const int v_elong_left_right_size = s_k_elong_left_right.size();
  List v_elong_left_right(v_elong_left_right_size); 
  for (int i = 0; i < v_elong_left_right_size; i++) {
    v_elong_left_right[i] = as<List>(compute_motif(v_dom_elong_left_right[i], s_k_elong_left_right[i], p_k, Y, m, use0, use1));
  }
  // v_elong_left_right is a List of List of two vectors v0 and v1.
  
  // create a LogicalVector start_with_NA whose elements are true iff the correspodent elemet of v_elong_left_right has lenght >2
  LogicalVector start_with_NA(v_elong_left_right_size);
  auto check_length = [](const List & v_elong_left_right_i){return(v_elong_left_right_i.size()>2);};
  std::transform(v_elong_left_right.begin(),v_elong_left_right.end(),start_with_NA.begin(),check_length);
  // filter the domain, centroid and shifts that are not in NA positions
  //!!std::vector<LogicalVector> new_v_dom; // @TODO: avoid unuseless copies
  List new_v_dom;
  List new_v;
  List new_s;
  auto not_NA_index = std::views::iota(0,v_elong_left_right_size) 
    | std::views::filter([&start_with_NA](int index_j){return(!start_with_NA[index_j]);});
  std::for_each(not_NA_index.begin(),not_NA_index.end(),
                [&new_v_dom,&new_v,
                 &new_s,&v_elong_left_right,
                 &v_dom_elong_left_right,
                 &s_k_elong_left_right](int j){new_v.push_back(as<List>(v_elong_left_right[j]));
                                new_v_dom.push_back(v_dom_elong_left_right[j]);
                                new_s.push_back(s_k_elong_left_right[j]);});
// new_v is a list of list of two vectors
// new_s is a list of list of two vectors
// new_v_dom is a list of list of two vectors

  v_dom_elong_left_right = new_v_dom;
  v_elong_left_right = new_v;
  s_k_elong_left_right = new_s;
  // compute performance index before elongation
  double Jk_before = as<double>(compute_Jk(v_new_k, s_k, p_k, Y, alpha, w, m, use0, use1));
  // compute performance indexes for all possible elongations
  vec c_k_after(v_elong_left_right_size);
  vec Jk_after(v_elong_left_right_size);
  for (uword i = 0; i < v_elong_left_right_size; i++) {
    const uvec& domain_elong = as<uvec>(domain(v_elong_left_right[i], use0));
    c_k_after(i) = floor(domain_elong.n_elem*(1 - max_gap));
    if (c_k_after(i) < c) {
      c_k_after(i) = c;
    }
    
    //!!const LogicalVector & v_i = v_dom_elong_left_right[i];
    const List& v_i = v_elong_left_right[i];
    const vec & s_i = s_k_elong_left_right[i];
    double c_i = c_k_after(i);
    Jk_after(i) = as<double>(compute_Jk(v_i, s_i, p_k, Y, alpha, w, m,use0 , use1,c_i,keep_k));
  }
  // find the best elongation in terms of perf. index
  vec diff_perc = ((Jk_after-Jk_before)/Jk_before);
  uword best_elong = index_min(diff_perc); //@TODO : ci sarebbe da fare un controllo su questo index_min se non c'? un min come fa la prof
  bool elongate = (diff_perc)(best_elong) < deltaJk_elong;
  // evaluate if elongate or not
  if(elongate) {
    return List::create(Named("v_new") = v_elong_left_right[best_elong],
                        Named("v_dom") = v_dom_elong_left_right[best_elong],
                        Named("s_k")   = s_k_elong_left_right[best_elong]);
  } else {
    return List::create(Named("v_new") = v_new_k,
                        Named("v_dom") = v_dom_k,
                        Named("s_k")   = s_k);
  }
  
}


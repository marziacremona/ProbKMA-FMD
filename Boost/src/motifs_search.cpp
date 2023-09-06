#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

class Comparator {
private:
  const Rcpp::NumericVector& ref;
  
  bool is_na(double x) const 
  {
    return Rcpp::traits::is_na<REALSXP>(x);    
  }
  
public:
  Comparator(const Rcpp::NumericVector& ref_)
    : ref(ref_)
  {}
  
  bool operator()(const int ilhs, const int irhs) const
  {
    double lhs = ref[ilhs], rhs = ref[irhs]; 
    if (is_na(lhs)) return false;
    if (is_na(rhs)) return true;
    return lhs < rhs;
  }
};

Rcpp::NumericVector avg_rank(Rcpp::NumericVector x)
{
  R_xlen_t sz = x.size();
  Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
  std::sort(w.begin(), w.end(), Comparator(x));
  
  Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      r[w[i + k]] = i + (n + 1) / 2.;
    }
  }
  
  return r;
}


template <int RTYPE>
IntegerVector order_impl(const Vector<RTYPE>& x, bool desc) {
  auto n = x.size();
  IntegerVector idx = no_init(n);
  std::iota(idx.begin(), idx.end(), static_cast<size_t>(1));
  if (desc) {
    auto comparator = [&x](size_t a, size_t b){ return x[a - 1] > x[b - 1]; };
    std::stable_sort(idx.begin(), idx.end(), comparator);
  } else {
    auto comparator = [&x](size_t a, size_t b){ return x[a - 1] < x[b - 1]; };
    std::stable_sort(idx.begin(), idx.end(), comparator);
    // simulate na.last
    size_t nas = 0;
    for (size_t i = 0; i < n; ++i, ++nas)
      if (!Vector<RTYPE>::is_na(x[idx[i] - 1])) break;
      std::rotate(idx.begin(), idx.begin() + nas, idx.end());
  }
  return idx;
}

IntegerVector order2(SEXP x, bool desc = false) {
  switch(TYPEOF(x)) {
  case INTSXP: return order_impl<INTSXP>(x, desc);
  case REALSXP: return order_impl<REALSXP>(x, desc);
  case STRSXP: return order_impl<STRSXP>(x, desc);
  default: stop("Unsupported type.");
  }
}



// [[Rcpp::export]]
List motifs_search(const List & Y, // list of matrices
                   const uvec& V,  // list of matrices 
                   const List& V_dom, // list of LogicalVector
                   const uvec V_length, // vector of V_dom_i's length 
                   const mat& P_clean,
                   const mat& D_clean,
                   const uvec& V_hclust, //vector of indices related to clusters 
                   const double alpha, 
                   const bool use0, const bool use1,
                   const vec& w, 
                   const vec& c,
                   const double max_gap,  
                   const double d, // n_col of Y[[0]]
                   const double N, // n_row of D
                   const double K, // n_col of D
                   const double R_all,
                   const vec R_m,
                   bool use_real_occurrences,
                   double length_diff) 
{
  unsigned int n_hclust = V_hclust.max();
  if(use_real_occurrences)
  {
    uvec c_k = floor(V_length*(1-max_gap));
    uvec index = find(c_k < c);
    for(auto i:index)
      c_k[i] = c[i];
    vec V_R_m = R_m.elem(V_hclust);
  // find occurrences
    Function find_occurrences(".find_occurrences");
    List V_occurrences(std::max(V.size(),V_R_m.size()));
    NumericVector not_null;
    for(int i = 0; i < V_occurrences.size();++i) // i am assuming that V_R_m has not necessary the same length of V and c_k
    {
      V_occurrences[i] = find_occurrences(V[i%V.size()],Y,
                                          V_R_m[i%V_R_m.size()],
                                          alpha,w,c_k[i%c_k.size()],
                                          use0,use1); // return a list of matrix or empty vector c()
      if(Nullable<mat>(V_occurrences[i]).isNotNull())
        not_null.push_back(i);
    }
    List V_new(not_null.size());
    uvec V_lenght_new(not_null.size());
    vec V_R_m_new(not_null.size());
    uvec V_hclust_new(not_null.size());
    for(auto i:not_null)
    {
      V_new[i] = V[i];
      V_lenght_new[i] = V_length[i];
      V_R_m_new[i] = V_R_m[i];
      V_hclust_new[i] = V_hclust[i];
    }
// select candidate motifs in each group ##############################################################
    uvec index_i;
    IntegerVector V_frequencies_i;
    IntegerVector V_mean_diss_i;
    for(int i = 1;i <= n_hclust;++i)
    {
      index_i = find(V_hclust_new==i);
      for(auto j:index_i)
      {
        V_frequencies_i.push_back(as<mat>(V_occurrences[j]).n_rows); 
        V_mean_diss_i.push_back(mean(as<mat>(V_occurrences[j]).col(2))); 
      }
//  order based of frequency and average distance
       //Rcpp::rank(V_frequencies_i);
    }
  }
  return List::create();
}

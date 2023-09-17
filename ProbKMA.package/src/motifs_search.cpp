#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <string>
#include <list>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

template<typename T>
class Comparator {
private:
  const T& ref;
  
  bool is_na(double x) const 
  {
    return Rcpp::traits::is_na<REALSXP>(x);    
  }
  
public:
  Comparator(const T& ref_)
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

template<typename T>
vec avg_rank(const T& x)
{
  R_xlen_t sz = x.size();
  Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
  std::sort(w.begin(), w.end(), Comparator(x));
  
  vec r;
  r.set_size(sz);
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      r[w[i + k]] = i + (n + 1) / 2.;
    }
  }
  
  return r;
}

template<typename V,typename T>
T order2(const V& x, bool desc = false) {
  std::size_t n = x.size();
  T idx;
  if constexpr(std::is_same<V, arma::uvec>::value)
  {
    idx.set_size(n); 
  }
  else
  {
    idx = T(n);
  }
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
      if (!Vector<REALSXP>::is_na(x[idx[i] - 1])) break;
      std::rotate(idx.begin(), idx.begin() + nas, idx.end());
  }
  return idx;
}


// [[Rcpp::export]]
List motifs_search_cpp(const List& Y, // list of list of matrices
                   const List& V,  // list of list of matrices 
                   const List& V0_clean,
                   const List& V1_clean,
                   const List& V_dom, // list of LogicalVector
                   const vec& V_length, // vector of V_dom_i's length 
                   const mat& P_clean,
                   const mat& D_clean,
                   uvec V_hclust, //vector of indices related to clusters 
                   const double alpha, 
                   const bool use0, const bool use1,
                   const vec& w, 
                   const vec& c,
                   const double max_gap,  
                   const double d, // n_col of Y[[0]]
                   const double N, // n_row of D
                   const double K, // n_col of D
                   const double R_all,
                   const vec& R_m,
                   bool use_real_occurrences,
                   double length_diff) 
{
  // transform index from R to c++
  unsigned int n_hclust = V_hclust.max();
  V_hclust -= 1;
  
  // prepare output data_structure
  uvec select(n_hclust*V_hclust.size()); // TODO:C'è da controllare che select non sia nullo altrimenti è un vettore vuoto e quindi restituisco risultati nulli 
  std::size_t V_size = V.size();
  List V_final;
  List V_occurrences(V_size);
  vec V_length_final;
  vec V_R_m;
  NumericVector V_frequencies_final;
  NumericVector V_mean_diss_final;
  IntegerVector index_final;
  Function find_occurrences_cpp("find_occurrences_cpp");
  int count = 0;

  if (use_real_occurrences)
  {
    uvec c_k = conv_to<uvec>::from(floor(V_length*(1-max_gap)));
    uvec index = find(c_k < c);
    std::for_each(index.begin(),index.end(),[&c_k,&c](uword i){c_k[i] = c[i];}); // provare togliere return
    V_R_m = arma::conv_to<vec>::from(R_m.elem(V_hclust));
  
  // find occurrences
    uvec not_null(V_size);
    for(int i = 0; i < V_size;++i) // i am assuming that V_R_m has not necessary the same length of V and c_k
    {
      V_occurrences[i] = find_occurrences_cpp(V[i],Y,
                                          V_R_m[i],
                                          alpha,w,c_k[i], 
                                          use0,use1); // return a list of matrix or empty vector c()
      if(Nullable<mat>(V_occurrences[i]).isNotNull())
        not_null[count++] = i;
    }
    not_null.resize(count); // resize the vector 
    V_final = V[as<IntegerVector>(wrap(not_null))];
    V_length_final = arma::conv_to<vec>::from(V_length.elem(not_null));
    V_R_m = arma::conv_to<vec>::from(V_R_m.elem(not_null));
    V_hclust = arma::conv_to<uvec>::from(V_hclust.elem(not_null));
  
    // select candidate motifs in each group
    count = 0;
    for(int i = 0;i < n_hclust;++i)
    {
      uvec index_i = find(V_hclust==i);
      std::size_t index_i_size = index_i.size();
      auto range_rows = index_i | std::views::transform([&V_occurrences](uword j)
        {return as<mat>(V_occurrences[j]).n_rows;});
      auto range_mean = index_i | std::views::transform([&V_occurrences](uword j)
        {return mean(as<mat>(V_occurrences[j]).col(2));});
      NumericVector V_frequencies_i(range_rows.begin(),range_rows.end());
      NumericVector V_mean_diss_i(range_mean.begin(),range_mean.end());

      // order based of frequency and average distanceà
      const vec& avgR = avg_rank(-1 * V_frequencies_i)+avg_rank(V_mean_diss_i);
      IntegerVector V_order_i = order2<vec,IntegerVector>(avgR) - 1;
      
      V_frequencies_i=V_frequencies_i[V_order_i];
      V_mean_diss_i=V_mean_diss_i[V_order_i];
      uvec index_i_ordered=index_i.elem(as<uvec>(V_order_i));
      vec V_length_i=V_length.elem(index_i_ordered);
      
      // select motifs to keep
      uvec keep(index_i_size,fill::ones);
      for(int i = 0;i<index_i_size;++i)
      {
        if(keep[i])
        {
         select[count++] = index_i_ordered[i];
         keep[i] = 0;
         // motifs with length different enough from length of selected motif
         const uvec& lhs = V_length_i < (static_cast<double>(V_length_i[i])*(1-length_diff));
         const uvec& rhs = V_length_i > (static_cast<double>(V_length_i[i])*(1+length_diff));
         keep = keep && (keep && (lhs || rhs));
        }
      }
    }
    select.resize(count);
    V_final = V_final[as<IntegerVector>(wrap((select)))];
    V_length_final = V_length_final.elem(select);
    V_R_m = V_R_m.elem(select);
    V_occurrences = V_occurrences[as<IntegerVector>(wrap((select)))];
    int V_occurrences_size = count; // set the new size 

    auto range_rows = std::views::iota(0,V_occurrences_size) | std::views::transform([&V_occurrences](int j)
    {return as<mat>(V_occurrences[j]).n_rows;});
    auto range_mean = std::views::iota(0,V_occurrences_size) | std::views::transform([&V_occurrences](int j)
    {return mean(as<mat>(V_occurrences[j]).col(2));});
    
    NumericVector V_frequencies(range_rows.begin(),range_rows.end());
    NumericVector V_mean_diss(range_mean.begin(),range_mean.end());
    
    const vec& avgR = avg_rank(-1 * V_frequencies)+avg_rank(V_mean_diss);
    IntegerVector V_order = order2<vec,IntegerVector>(avgR) - 1;
    V_final = V_final[V_order];
    V_occurrences = V_occurrences[V_order];
    V_length_final = arma::conv_to<vec>::from(V_length_final.elem(as<uvec>(V_order)));
    V_R_m = arma::conv_to<vec>::from(V_R_m.elem(as<uvec>(V_order)));
    V_frequencies_final = V_frequencies[V_order];
    V_mean_diss_final = V_mean_diss[V_order];
    index_final = as<IntegerVector>(wrap(arma::conv_to<uvec>::from(not_null.elem(select.elem(as<uvec>(V_order))))));//not_null[select][V_order]
  }
  else
  {
    for(int i_hclust = 0;i_hclust < n_hclust;++i_hclust)
    {
      const uvec& index_i = find(V_hclust==i_hclust);
      std::size_t index_i_size = index_i.size();
      const auto& V_D_i = D_clean.cols(index_i);
      const auto& Logic = V_D_i <= R_m[i_hclust];
      vec V_frequencies_approx_i = arma::conv_to<vec>::from(sum(Logic, 0).t());
      vec V_mean_diss_approx_i = conv_to<vec>::from(sum(V_D_i % Logic, 0).t() / V_frequencies_approx_i);
      
      vec avgR = avg_rank(arma::conv_to<vec>::from(-1 * V_frequencies_approx_i))+avg_rank(V_mean_diss_approx_i);
      uvec V_order_i = order2<vec,uvec>(avgR) - 1; // order starts to count from 1 not zero
      
      V_frequencies_approx_i=V_frequencies_approx_i.elem(V_order_i); 
      V_mean_diss_approx_i=V_mean_diss_approx_i.elem(V_order_i);
      uvec index_i_ordered=index_i.elem(V_order_i);
      vec V_length_i=V_length.elem(index_i_ordered); 
      
      // select motifs to keep
      uvec keep(index_i_size,fill::ones);
      for(int i = 0;i<index_i_size;++i)
      {
        if(keep[i])
        {
          select[count++] = index_i_ordered[i];
          keep[i] = 0;
          // motifs with length different enough from length of selected motif
          const uvec& lhs = V_length_i < (static_cast<double>(V_length_i[i])*(1-length_diff));
          const uvec& rhs = V_length_i > (static_cast<double>(V_length_i[i])*(1+length_diff));
          keep = keep && (keep && (lhs || rhs));
        }
      }
    }
      select.resize(count);
      V_final = V[as<IntegerVector>(wrap((select)))];
      V_size = count;
      V_length_final = V_length.elem(select); 
      V_R_m = R_m.elem(V_hclust.elem(select));
      uvec c_k = conv_to<uvec>::from(floor(V_length_final*(1-max_gap)));
      const uvec& index_logic = find(c_k < c.elem(select)); 
      c_k.elem(index_logic) = arma::conv_to<arma::uvec>::from(c.elem(select.elem(index_logic)));
  
      // find occurrences
      count = 0;
      uvec not_null(V_size);      
      for(int i = 0; i < V_size;++i) 
      {
        V_occurrences[i] = find_occurrences_cpp(V_final[i],Y,
                                            V_R_m[i],
                                            alpha,w,c_k[i], // non so se c_k[i%c_k.size()]
                                            use0,use1); // return a list of matrix or empty vector c()
        if(Nullable<mat>(V_occurrences[i]).isNotNull())
          not_null[count++] = i;
      }
      
      not_null.resize(count); // resize the vector 
      V_final = V_final[as<IntegerVector>(wrap((not_null)))];
      V_occurrences = V_occurrences[as<IntegerVector>(wrap((not_null)))];
      V_length_final = V_length_final.elem(not_null);
      V_R_m = arma::conv_to<arma::vec>::from(V_R_m.elem(not_null));
      
      auto range_rows = std::views::iota(0,count) | std::views::transform([&V_occurrences](int j)
      {return as<mat>(V_occurrences[j]).n_rows;});
      auto range_mean = std::views::iota(0,count) | std::views::transform([&V_occurrences](int j)
      {return mean(as<mat>(V_occurrences[j]).col(2));});
      NumericVector V_frequencies(range_rows.begin(),range_rows.end());
      NumericVector V_mean_diss(range_mean.begin(),range_mean.end());

      const vec& avgR = avg_rank(-1 * V_frequencies)+avg_rank(V_mean_diss);
      IntegerVector V_order = (order2<vec,IntegerVector>(avgR) - 1);
      V_final = V_final[V_order];
      V_occurrences = V_occurrences[V_order];
      V_length_final = V_length_final.elem(as<uvec>(V_order));
      V_R_m = arma::conv_to<arma::vec>::from(V_R_m.elem(as<uvec>(V_order)));
      V_frequencies_final = V_frequencies[V_order];
      V_mean_diss_final = V_mean_diss[V_order];
      index_final = as<IntegerVector>(wrap(arma::conv_to<uvec>::from(select.elem(not_null.elem(as<uvec>(V_order)))))); //select[not_null][V_order]
  }
  return List::create(Named("V") = V_final,
                      Named("V_length")=IntegerVector(V_length_final.begin(),V_length_final.end()),
                      Named("V_occurrences")=V_occurrences,
                      Named("V_frequencies")=V_frequencies_final,
                      Named("V_mean_diss")=V_mean_diss_final,
                      Named("R_motifs")=NumericVector(V_R_m.begin(),V_R_m.end()),
                      Named("index_final")=index_final);
}




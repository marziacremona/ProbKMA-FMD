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
  const NumericVector& ref;
  
  bool is_na(double x) const 
  {
    return Rcpp::traits::is_na<REALSXP>(x);    
  }
  
public:
  Comparator(const NumericVector& ref_)
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

Rcpp::NumericVector avg_rank(NumericVector x)
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
List motifs_search_custom_cpp(const List& Y, // list of list of matrices
                   const List& Y0,
                   const List& Y1,
                   const List& V,  // list of list of matrices 
                   const List& V0_clean,
                   const List& V1_clean,
                   const List& V_dom, // list of LogicalVector
                   const uvec V_length, // vector of V_dom_i's length 
                   const mat& P_clean,
                   const mat& D_clean,
                   uvec V_hclust, //vector of indices related to clusters 
                   const std::string& diss,
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
  V_hclust = V_hclust - 1;
  List V_new_final;
  uvec V_length_final;
  List V_occurrences_final;
  NumericVector V_frequencies_final;
  NumericVector V_mean_diss_final;
  vec V_R_m_final;
  IntegerVector index_final;
  if(use_real_occurrences)
  {
    uvec c_k = floor(V_length*(1-max_gap));
    uvec index = find(c_k < c);
    for(auto i:index)
      c_k[i] = c[i];
    vec V_R_m = R_m.elem(V_hclust);
  // find occurrences
    Function find_occurrences(".find_occurrences");
    std::size_t V_R_m_size = V_R_m.size();
    std::size_t V_size = V.size();
    List V_occurrences(std::max(V_size,V_R_m_size));
    IntegerVector not_null;
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
    uvec V_length_new(not_null.size());
    vec V_R_m_new(not_null.size());
    uvec V_hclust_new(not_null.size());
    for(auto i:not_null)
    {
      V_new[i] = V[i];
      V_length_new[i] = V_length[i];
      V_R_m_new[i] = V_R_m[i];
      V_hclust_new[i] = V_hclust[i];
    }
// select candidate motifs in each group ##############################################################
    IntegerVector select; // TODO:C'è da controllare che select non sia nullo altrimenti è un vettore vuoto e quindi restituisco risultati nulli 
    
    for(int i = 0;i < n_hclust;++i)
    {
      uvec index_i = find(V_hclust_new==i);
      NumericVector V_frequencies_i(index_i.size());
      NumericVector V_mean_diss_i(index_i.size());
      for(int i = 0;auto j:index_i)
      {
        V_frequencies_i[i] = as<mat>(V_occurrences[j]).n_rows; 
        V_mean_diss_i[i++] = mean(as<mat>(V_occurrences[j]).col(2)); 
      }
//  order based of frequency and average distance

      NumericVector avgR = avg_rank(-1 * V_frequencies_i)+avg_rank(V_mean_diss_i);
      IntegerVector V_order_i = order2(avgR) - 1;
      
      V_frequencies_i=V_frequencies_i[V_order_i];
      V_mean_diss_i=V_mean_diss_i[V_order_i];
      uvec index_i_ordered=index_i.elem(as<uvec>(V_order_i));
      uvec V_length_i=V_length.elem(index_i_ordered);
      
// select motifs to keep

      LogicalVector keep=rep(true,index_i.size());
      
      for(int i = 0;i<index_i.size();++i)
      {
        if(keep[i])
        {
          select.push_back(index_i_ordered[i]);
          keep[i] = false;
        // motifs with length different enough from length of selected motif
        LogicalVector lhs = as<LogicalVector>(wrap(arma::conv_to<vec>::from(V_length_i) < (static_cast<double>(V_length_i[i])*(1-length_diff))));
        LogicalVector rhs = as<LogicalVector>(wrap(arma::conv_to<vec>::from(V_length_i) > (static_cast<double>(V_length_i[i])*(1+length_diff))));
        LogicalVector length_diff_enough = keep & (lhs | rhs);
        
        for(int i = 0;i<length_diff_enough.size();++i)
          if(!length_diff_enough[i]) keep[i] = false;
        }
      }
    }
    V_new = V_new[select];
    V_length_new = V_length_new.elem(as<uvec>(select));
    V_R_m_new = V_R_m_new.elem(as<uvec>(select));
    V_occurrences = V_occurrences[select];
    NumericVector V_frequencies(V_occurrences.size());
    NumericVector V_mean_diss(V_occurrences.size());
    for(int j = 0;j<V_occurrences.size();++j)
    {
      V_frequencies[j] = as<mat>(V_occurrences[j]).n_rows; 
      V_mean_diss[j] = mean(as<mat>(V_occurrences[j]).col(2)); 
    }
    const NumericVector& avgR = avg_rank(-1 * V_frequencies)+avg_rank(V_mean_diss);
    IntegerVector V_order = order2(avgR) - 1;
    V_new = V_new[V_order];
    V_occurrences = V_occurrences[V_order];
    V_length_new = V_length_new.elem(as<uvec>(V_order));
    V_R_m_new = V_R_m_new.elem(as<uvec>(V_order));
    V_frequencies = V_frequencies[V_order];
    V_mean_diss = V_mean_diss[V_order];
    const IntegerVector not_null_temp = not_null[select];;
    index_final = not_null_temp[V_order];
    
    V_new_final = V_new;
    V_length_final = V_length_new;
    V_occurrences_final = V_occurrences;
    V_frequencies_final = V_frequencies;
    V_mean_diss_final = V_mean_diss;
    V_R_m_final = V_R_m_new;

  }
  else
  {
    IntegerVector select; // TODO:C'è da controllare che select non sia nullo altrimenti è un vettore vuoto e quindi restituisco risultati nulli 
    
    for(int i_hclust = 0;i_hclust < n_hclust;++i_hclust)
    {
      uvec index_i = find(V_hclust==i_hclust);
      mat V_D_i = arma::conv_to<arma::mat>::from(D_clean.cols(index_i));
      const umat& Logic = V_D_i <= R_m[i_hclust];
      uvec V_frequencies_approx_i_ = sum(Logic).t(); // TODO:da convertire direttamente in NumericVector
 
      NumericVector V_frequencies_approx_i(V_frequencies_approx_i_.begin(),V_frequencies_approx_i_.end());
    
      NumericVector V_mean_diss_approx_i(V_D_i.n_cols);

      for(int j = 0;j < V_D_i.n_cols;++j)
      {
        const vec& x = V_D_i.col(j);    ///
        uvec index = find(x <= R_m(i_hclust));
        V_mean_diss_approx_i[j] = mean(x.elem(index)); ///
      }
      NumericVector avgR = avg_rank(-1 * V_frequencies_approx_i)+avg_rank(V_mean_diss_approx_i);
      IntegerVector V_order_i = order2(avgR) - 1; // order starts to count from 1 not zero
     
      
      V_frequencies_approx_i=V_frequencies_approx_i[V_order_i]; //
      V_mean_diss_approx_i=V_mean_diss_approx_i[V_order_i];
      uvec index_i_ordered=index_i.elem(as<uvec>(V_order_i));
      uvec V_length_i=V_length.elem(index_i_ordered);
      // select motifs to keep
      LogicalVector keep=rep(true,index_i.size());
      for(int i = 0;i<index_i.size();++i)
      {
        if(keep[i])
        {
          select.push_back(index_i_ordered[i]);
          keep[i] = false;
          // motifs with length different enough from length of selected motif
          LogicalVector lhs = as<LogicalVector>(wrap(arma::conv_to<vec>::from(V_length_i) < (static_cast<double>(V_length_i[i])*(1-length_diff))));
          LogicalVector rhs = as<LogicalVector>(wrap(arma::conv_to<vec>::from(V_length_i) > (static_cast<double>(V_length_i[i])*(1+length_diff))));
          LogicalVector length_diff_enough = keep & (lhs | rhs);
          
          for(int i = 0;i<length_diff_enough.size();++i)
            if(!length_diff_enough[i]) keep[i] = false;
          
        }
      }
    }
      List V_new = V[select];
      uvec V_length_new = V_length.elem(as<uvec>(select));
      vec V_R_m_new = R_m.elem(V_hclust.elem(as<uvec>(select)));
      /// find candidate motifs in the curves ####################################################################
      uvec c_k = floor(V_length_new*(1-max_gap));
      uvec index_logic = find(c_k < c.elem(as<uvec>(select))); 
      const uvec& temp = arma::conv_to<arma::uvec>::from(c.elem(as<uvec>(select)));
      const uvec& temp2 = arma::conv_to<arma::uvec>::from(temp.elem(index_logic));
      c_k.elem(index_logic) = arma::conv_to<arma::uvec>::from(temp.elem(index_logic));
      // find occurrences
      Function find_occurrences(".find_occurrences");
      std::size_t V_R_m_size = V_R_m_new.size();
      std::size_t V_size = V_new.size();
      List V_occurrences(std::max(V_size,V_R_m_size));
      IntegerVector not_null;
      for(int i = 0; i < V_occurrences.size();++i) // i am assuming that V_R_m has not necessary the same length of V and c_k
      {
        V_occurrences[i] = find_occurrences(V_new[i%V_size],Y,
                                            V_R_m_new[i%V_R_m_size],
                                                 alpha,w,c_k[i%c_k.size()],
                                                            use0,use1); // return a list of matrix or empty vector c()
        if(Nullable<mat>(V_occurrences[i]).isNotNull())
          not_null.push_back(i);
      }
      V_new = V_new[not_null];
      V_occurrences = V_occurrences[not_null];
      V_R_m_new = arma::conv_to<arma::vec>::from(V_R_m_new.elem(as<uvec>(not_null)));
      NumericVector V_frequencies(V_occurrences.size());
      NumericVector V_mean_diss(V_occurrences.size());
      for(int i = 0;i < V_occurrences.size();++i)
      {
        V_frequencies[i] = as<mat>(V_occurrences[i]).n_rows; 
        V_mean_diss[i] = mean(as<mat>(V_occurrences[i]).col(2)); 
      }
      const NumericVector& avgR = avg_rank(-1 * V_frequencies)+avg_rank(V_mean_diss);
      IntegerVector V_order = order2(avgR) - 1;
      V_new = V_new[V_order];
      V_occurrences = V_occurrences[V_order];
      V_length_new = V_length_new.elem(as<uvec>(V_order));
      V_R_m_new = arma::conv_to<arma::vec>::from(V_R_m_new.elem(as<uvec>(V_order)));
      V_frequencies = V_frequencies[V_order];
      V_mean_diss = V_mean_diss[V_order];
      const IntegerVector select_temp = select[not_null];
      index_final = select_temp[V_order];
      V_new_final = V_new;
      V_length_final = V_length_new;
      V_occurrences_final = V_occurrences;
      V_frequencies_final = V_frequencies;
      V_mean_diss_final = V_mean_diss;
      V_R_m_final = V_R_m_new;
  }
  List V0(V_new_final.size());
  List V1(V_new_final.size());
  if(diss=="d0_L2")
  {
    for(int i = 0; const List& v:V_new_final)
      V0[i++] = as<mat>(v[0]);
    
    V1 = V1_clean[index_final];
  }
  else if(diss=="d1_L2")
  {
    for(int i = 0; const List& v:V_new_final)
      V1[i++] = as<mat>(v[1]);
    
    V0 = V0_clean[index_final];
  }
  else
  {
    for(int i = 0; const List& v:V_new_final)
    {
      V0[i] = as<mat>(v[0]);
      V1[i++] = as<mat>(v[1]); 
    }
  }
  
  return List::create(Named("V0")=V0,Named("V1")=V1,
                      Named("V_length")=NumericVector(V_length_final.begin(),V_length_final.end()),Named("V_occurrences")=V_occurrences_final,
                      Named("V_frequencies")=V_frequencies_final,Named("V_mean_diss")=V_mean_diss_final,
                      Named("Y0")=Y0,Named("Y1")=Y1,Named("R_motifs")=V_R_m_final);
}

// V_frequencies IntegerVector
// V_R_m_final NumericVector
// V_length NumericVector
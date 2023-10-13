#include "motifs_search.h"

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

Rcpp::List motifs_search_cpp(const Rcpp::List& Y, // list of list of matrices
                             const Rcpp::List& V,  // list of list of matrices 
                             const Rcpp::List& V0_clean,
                             const Rcpp::List& V1_clean,
                             const Rcpp::List& V_dom, // list of LogicalVector
                             const arma::vec& V_length, // vector of V_dom_i's length 
                             const arma::mat& P_clean,
                             const arma::mat& D_clean,
                             arma::uvec V_hclust, //vector of indices related to clusters 
                             const double alpha, 
                             const bool use0, const bool use1,
                             const arma::vec& w, 
                             const arma::vec& c,
                             const double max_gap,  
                             const double d, // n_col of Y[[0]]
                             const double N, // n_row of D
                             const double K, // n_col of D
                             const double R_all,
                             const arma::vec& R_m,
                             bool use_real_occurrences,
                             double length_diff,
                             Rcpp::Function diss_d0_d1_L2, 
                             Rcpp::Function domain,
                             Rcpp::Function select_domain) 
{
  // transform index from R to C++
  unsigned int n_hclust = V_hclust.max();
  V_hclust -= 1;
  
  // prepare output data_structure
  arma::uvec select(n_hclust*V_hclust.size()); // TODO:C'è da controllare che select non sia nullo altrimenti è un vettore vuoto e quindi restituisco risultati nulli 
  std::size_t V_size = V.size();
  Rcpp::List V_final;
  Rcpp::List V_occurrences(V_size);
  arma::vec V_length_final;
  arma::vec V_R_m;
  Rcpp::NumericVector V_frequencies_final;
  Rcpp::NumericVector V_mean_diss_final;
  Rcpp::IntegerVector index_final;
  int count = 0;

  if (use_real_occurrences)
  {
    arma::uvec c_k = arma::conv_to<arma::uvec>::from(arma::floor(V_length*(1-max_gap)));
    arma::uvec index = find(c_k < c);
    std::for_each(index.begin(),index.end(),[&c_k,&c](arma::uword i){c_k[i] = c[i];}); // provare togliere return
    V_R_m = arma::conv_to<arma::vec>::from(R_m.elem(V_hclust));
  
  // find occurrences
    arma::uvec not_null(V_size);
    for(int i = 0; i < V_size;++i) // i am assuming that V_R_m has not necessary the same length of V and c_k
    {
      V_occurrences[i] = find_occurrences_cpp(V[i],Y,
                                          V_R_m[i],
                                          alpha,w,c_k[i], 
                                          use0,use1,
                                          diss_d0_d1_L2,
                                          domain,
                                          select_domain); // return a list of matrix or empty vector c()
      
      if(Rcpp::Nullable<arma::mat>(V_occurrences[i]).isNotNull())
        not_null[count++] = i;
    }
    not_null.resize(count); // resize the vector 
    V_final = V[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(not_null))];
    V_length_final = arma::conv_to<arma::vec>::from(V_length.elem(not_null));
    V_R_m = arma::conv_to<arma::vec>::from(V_R_m.elem(not_null));
    V_hclust = arma::conv_to<arma::uvec>::from(V_hclust.elem(not_null));
  
    // select candidate motifs in each group
    count = 0;
    for(int i = 0;i < n_hclust;++i)
    {
      arma::uvec index_i = arma::find(V_hclust==i);
      std::size_t index_i_size = index_i.size();
      auto range_rows = index_i | std::views::transform([&V_occurrences](arma::uword j)
      {return as<arma::mat>(V_occurrences[j]).n_rows;});
      auto range_mean = index_i | std::views::transform([&V_occurrences](arma::uword j)
      {return mean(as<arma::mat>(V_occurrences[j]).col(2));});
      Rcpp::NumericVector V_frequencies_i(range_rows.begin(),range_rows.end());
      Rcpp::NumericVector V_mean_diss_i(range_mean.begin(),range_mean.end());

      // order based of frequency and average distanceà
      const arma::vec& avgR = util::avg_rank<Rcpp::NumericVector,Comparator<Rcpp::NumericVector>>(-1 * V_frequencies_i)+util::avg_rank<Rcpp::NumericVector,Comparator<Rcpp::NumericVector>>(V_mean_diss_i);
      Rcpp::IntegerVector V_order_i = util::order2<arma::vec,Rcpp::IntegerVector>(avgR) - 1;
      
      V_frequencies_i=V_frequencies_i[V_order_i];
      V_mean_diss_i=V_mean_diss_i[V_order_i];
      arma::uvec index_i_ordered=index_i.elem(Rcpp::as<arma::uvec>(V_order_i));
      arma::vec V_length_i=V_length.elem(index_i_ordered);
      
      // select motifs to keep
      arma::uvec keep(index_i_size,arma::fill::ones);
      for(int i = 0;i<index_i_size;++i)
      {
        if(keep[i])
        {
         select[count++] = index_i_ordered[i];
         keep[i] = 0;
         // motifs with length different enough from length of selected motif
         const arma::uvec& lhs = V_length_i < (static_cast<double>(V_length_i[i])*(1-length_diff));
         const arma::uvec& rhs = V_length_i > (static_cast<double>(V_length_i[i])*(1+length_diff));
         keep = keep && (keep && (lhs || rhs));
        }
      }
    }
    select.resize(count);
    V_final = V_final[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap((select)))];
    V_length_final = V_length_final.elem(select);
    V_R_m = V_R_m.elem(select);
    V_occurrences = V_occurrences[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap((select)))];
    int V_occurrences_size = count; // set the new size 

    auto range_rows = std::views::iota(0,V_occurrences_size) | std::views::transform([&V_occurrences](int j)
    {return Rcpp::as<arma::mat>(V_occurrences[j]).n_rows;});
    auto range_mean = std::views::iota(0,V_occurrences_size) | std::views::transform([&V_occurrences](int j)
    {return mean(Rcpp::as<arma::mat>(V_occurrences[j]).col(2));});
    
    Rcpp::NumericVector V_frequencies(range_rows.begin(),range_rows.end());
    Rcpp::NumericVector V_mean_diss(range_mean.begin(),range_mean.end());
    
    const arma::vec& avgR = util::avg_rank<Rcpp::NumericVector,Comparator<Rcpp::NumericVector>>(-1 * V_frequencies)+util::avg_rank<Rcpp::NumericVector,Comparator<Rcpp::NumericVector>>(V_mean_diss);
    Rcpp::IntegerVector V_order = util::order2<arma::vec,Rcpp::IntegerVector>(avgR) - 1;
    V_final = V_final[V_order];
    V_occurrences = V_occurrences[V_order];
    V_length_final = arma::conv_to<arma::vec>::from(V_length_final.elem(Rcpp::as<arma::uvec>(V_order)));
    V_R_m = arma::conv_to<arma::vec>::from(V_R_m.elem(Rcpp::as<arma::uvec>(V_order)));
    V_frequencies_final = V_frequencies[V_order];
    V_mean_diss_final = V_mean_diss[V_order];
    index_final = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(arma::conv_to<arma::uvec>::from(not_null.elem(select.elem(Rcpp::as<arma::uvec>(V_order))))));//not_null[select][V_order]
  }
  else
  {
    for(int i_hclust = 0;i_hclust < n_hclust;++i_hclust)
    {
      const arma::uvec& index_i = find(V_hclust==i_hclust);
      std::size_t index_i_size = index_i.size();
      const auto& V_D_i = D_clean.cols(index_i);
      const auto& Logic = V_D_i <= R_m[i_hclust];
      arma::vec V_frequencies_approx_i = arma::conv_to<arma::vec>::from(sum(Logic, 0).t());
      arma::vec V_mean_diss_approx_i = arma::conv_to<arma::vec>::from(sum(V_D_i % Logic, 0).t() / V_frequencies_approx_i);
      
      arma::vec avgR = util::avg_rank<arma::vec,Comparator<arma::vec>>(arma::conv_to<arma::vec>::from(-1 * V_frequencies_approx_i))
                     + util::avg_rank<arma::vec,Comparator<arma::vec>>(V_mean_diss_approx_i);
      arma::uvec V_order_i = util::order2<arma::vec,arma::uvec>(avgR) - 1; // order starts to count from 1 not zero
      
      V_frequencies_approx_i=V_frequencies_approx_i.elem(V_order_i); 
      V_mean_diss_approx_i=V_mean_diss_approx_i.elem(V_order_i);
      arma::uvec index_i_ordered=index_i.elem(V_order_i);
      arma::vec V_length_i=V_length.elem(index_i_ordered); 
      
      // select motifs to keep
      arma::uvec keep(index_i_size,arma::fill::ones);
      for(int i = 0;i<index_i_size;++i)
      {
        if(keep[i])
        {
          select[count++] = index_i_ordered[i];
          keep[i] = 0;
          // motifs with length different enough from length of selected motif
          const arma::uvec& lhs = V_length_i < (static_cast<double>(V_length_i[i])*(1-length_diff));
          const arma::uvec& rhs = V_length_i > (static_cast<double>(V_length_i[i])*(1+length_diff));
          keep = keep && (keep && (lhs || rhs));
        }
      }
    }
  
      select.resize(count);
      V_final = V[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap((select)))];
      V_size = count;
      V_length_final = V_length.elem(select); 
      V_R_m = R_m.elem(V_hclust.elem(select));
      arma::uvec c_k = arma::conv_to<arma::uvec>::from(arma::floor(V_length_final*(1-max_gap)));
      const arma::uvec& index_logic = arma::find(c_k < c.elem(select)); 
      c_k.elem(index_logic) = arma::conv_to<arma::uvec>::from(c.elem(select.elem(index_logic)));
      
      // find occurrences
      count = 0;
      arma::uvec not_null(V_size); 
      for(int i = 0; i < V_size;++i) 
      {
        V_occurrences[i] = find_occurrences_cpp(V_final[i],Y,
                                                V_R_m[i],
                                                alpha,w,c_k[i], // non so se c_k[i%c_k.size()]
                                                use0,use1,
                                                diss_d0_d1_L2, 
                                                domain,
                                                select_domain); // return a list of matrix or empty vector c()
        if(Rcpp::Nullable<arma::mat>(V_occurrences[i]).isNotNull())
          not_null[count++] = i;
      }
      not_null.resize(count); // resize the vector 
      V_final = V_final[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap((not_null)))];
      V_occurrences = V_occurrences[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap((not_null)))];
      V_length_final = V_length_final.elem(not_null);
      V_R_m = arma::conv_to<arma::vec>::from(V_R_m.elem(not_null));
    
      auto range_rows = std::views::iota(0,count) | std::views::transform([&V_occurrences](int j)
      {return Rcpp::as<arma::mat>(V_occurrences[j]).n_rows;});
      auto range_mean = std::views::iota(0,count) | std::views::transform([&V_occurrences](int j)
      {return mean(Rcpp::as<arma::mat>(V_occurrences[j]).col(2));});
      Rcpp::NumericVector V_frequencies(range_rows.begin(),range_rows.end());
      Rcpp::NumericVector V_mean_diss(range_mean.begin(),range_mean.end());
      
      const arma::vec& avgR = util::avg_rank<Rcpp::NumericVector,Comparator<Rcpp::NumericVector>>(-1 * V_frequencies)
                            +util::avg_rank<Rcpp::NumericVector,Comparator<Rcpp::NumericVector>>(V_mean_diss);
      Rcpp::IntegerVector V_order = (util::order2<arma::vec,Rcpp::IntegerVector>(avgR) - 1);
      V_final = V_final[V_order];
      V_occurrences = V_occurrences[V_order];
      V_length_final = V_length_final.elem(Rcpp::as<arma::uvec>(V_order));
      V_R_m = arma::conv_to<arma::vec>::from(V_R_m.elem(Rcpp::as<arma::uvec>(V_order)));
      V_frequencies_final = V_frequencies[V_order];
      V_mean_diss_final = V_mean_diss[V_order];
      index_final = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(arma::conv_to<arma::uvec>::from(select.elem(not_null.elem(Rcpp::as<arma::uvec>(V_order)))))); //select[not_null][V_order]
  }
  V_hclust += 1;
  return Rcpp::List::create(Rcpp::Named("V") = V_final,
                            Rcpp::Named("V_length")= Rcpp::IntegerVector(V_length_final.begin(),V_length_final.end()),
                            Rcpp::Named("V_occurrences")=V_occurrences,
                            Rcpp::Named("V_frequencies")=V_frequencies_final,
                            Rcpp::Named("V_mean_diss")=V_mean_diss_final,
                            Rcpp::Named("R_motifs")= Rcpp::NumericVector(V_R_m.begin(),V_R_m.end()),
                            Rcpp::Named("index_final")=index_final);
}




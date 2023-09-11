// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include <numeric>
#include <ranges>
#include <algorithm>
using namespace Rcpp;
using namespace arma;
#include <algorithm>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

template<typename T>
decltype(auto) combn2(const T & y){ 
  int n = y.n_elem;
  uvec v(n,fill::zeros);  
  v(0) = 1;
  v(1) = 1;
  std::size_t l = 0;
  if constexpr(std::is_same<typename T::value_type, vec::value_type>::value) {
    mat result(2,n*(n-1)/2, fill::zeros); 
    std::size_t k;
    do {
      k = 0;
      auto filter_index = std::views::iota(0,n) 
        | std::views::filter([&v](int i){return v(i);});
      for(auto i: filter_index) 
        result(k++,l) = y(i); 
      l++;
    } while (std::prev_permutation(v.begin(), v.end())); 
    return result;
  } else if constexpr(std::is_same<typename T::value_type, uvec::value_type>::value){
    umat result(2,n*(n-1)/2, fill::zeros); 
    std::size_t k;
    do {
      k = 0;
      auto filter_index = std::views::iota(0,n) 
        | std::views::filter([&v](int i){return v(i);});
      for(auto i: filter_index) 
        result(k++,l) = y(i); 
      l++;
    } while (std::prev_permutation(v.begin(), v.end()));
    return result;
  } else if constexpr(std::is_same<typename T::value_type, ivec::value_type>::value){
    imat result(2,n*(n-1)/2, fill::zeros); 
    std::size_t k;
    do {
      k = 0;
      auto filter_index = std::views::iota(0,n) 
        | std::views::filter([&v](int i){return v(i);});
      for(auto i: filter_index) 
        result(k++,l) = y(i); 
      l++;
    } while (std::prev_permutation(v.begin(), v.end())); 
    return result;
  }
}

template<class T>
T repLem(const T & v,
         const ivec & times){
  T result(accu(times));
  std::size_t k = 0;
  std::size_t times_size = times.size();
  // check times_size == v.size()
  for (std::size_t i = 0; i < times_size; ++i) 
    for (std::size_t j = 0; j < times[i]; ++j)
      result(k++) = v(i);
  return result;
}

// [[Rcpp::export]]
List probKMA_silhouette_rcpp(const List & probKMA_results,
                             bool align){ // @TODO = set to FALSE by default
  // da ripensare con un qualche sorting su P_clean
  double alpha;
  bool use0, use1;
  unsigned int Y_size;
  
  std::string diss = as<std::string>(probKMA_results["diss"]);
  
  if (diss == "d0_L2" || diss == "d0_d1_L2") {
    const List & Y0 = probKMA_results["Y0"];
    Y_size = Y0.size();
  } else {
    const List & Y1 = probKMA_results["Y1"];
    Y_size = Y1.size();
  }
  
  List Y(Y_size);
  
  // to declare as external functions and used also in probKMA
  if (diss == "d0_L2") {
    alpha = 0;
    use0 = true;
    use1 = false;
    const List & Y0 = probKMA_results["Y0"];
    std::transform(Y0.cbegin(),
                   Y0.cend(),
                   Y.begin(),
                   [](const vec & y0) 
                     { return List::create(_["y0"] = y0,
                                           _["y1"] = R_NilValue); });
  } else if (diss == "d1_L2") {
    alpha = 1;
    use0 = false;
    use1 = true;
    const List & Y1 = probKMA_results["Y1"];
    std::transform(Y1.cbegin(),
                   Y1.cend(),
                   Y.begin(),
                   [](const vec & y1) 
                     { return List::create(_["y0"] = R_NilValue,
                                           _["y1"] = y1); });
  } else if (diss == "d0_d1_L2") {
    alpha = as<double>(probKMA_results["alpha"]);
    use0 = true;
    use1 = true;
    const List & Y0 = probKMA_results["Y0"];
    const List & Y1 = probKMA_results["Y1"];
    std::transform(Y0.cbegin(),
                   Y0.cend(),
                   Y1.cbegin(),
                   Y.begin(),
                   [](const vec & y0, const vec & y1) 
                     {  return List::create(_["y0"] = y0,
                                            _["y1"] = y1); });
  }
  
  const vec & w = probKMA_results["w"];
  const List & first_y = Y[0];
  const mat & first_y0 = use0 ? as<mat>(first_y[0]) : as<mat>(first_y[1]); 
  const unsigned int d = first_y0.n_cols;
  const unsigned int N = first_y0.n_rows;
  const unsigned int K = probKMA_results["K"];
  const List & motifs = use0 ? probKMA_results["V0"] : probKMA_results["V1"];
  
  std::vector<uvec> V_dom(K); // each element will be an uvec
  ivec V_length(K);
  unsigned int i = 0;
  unsigned int j;
  
  // V_dom for each centroid contains an uvec with 0 if all the elements of the row of the centroid are NA
  for(unsigned int i = 0; i < K; ++i){
    const mat & v = motifs[i];
    uvec v_dom_k(v.n_rows);
    j = 0;
    v.each_row([&v_dom_k,&j](const vec& v_row_j){ const uvec & nan_v_row_j = find_nan(v_row_j);
                                                  v_dom_k[j++] = nan_v_row_j.n_elem != v_row_j.n_elem;});
    V_length(i) = v.n_rows;
    V_dom[i] = v_dom_k;
  }
 
 // extract pieces of curves that fall in the different motifs
 const imat & P_clean = probKMA_results["P_clean"];
 List curves_in_motifs(P_clean.n_cols); //TODO: change List
 for (unsigned int k = 0; k < P_clean.n_cols; ++k) {
   const ivec & P_clean_k = P_clean.col(k);
   const uvec & P_clean_1 = find(P_clean_k == 1); // to check in case 1 doesn't exist in the col
   curves_in_motifs[k] = P_clean_1;
 }
 
 // if(!is.null(ncol(curves_in_motifs))) this if condition makes no sense since curves_in_motifs is a list
 //   curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs))) to be implemented and understood
 
 const ivec & curves_in_motifs_number=sum(P_clean,0).t();
 
 const imat & S_clean = probKMA_results["S_clean"];
 List S_clean_k(S_clean.n_cols); // TODO: change List
 for (unsigned int k=0; k < S_clean.n_cols; ++k){
  const ivec & col_S_clean = S_clean.col(k);
  const uvec & curves_in_motifs_k = curves_in_motifs[k];
  const ivec & indexed_S_clean = col_S_clean.elem(curves_in_motifs_k);
  S_clean_k[k] = indexed_S_clean;
 }
 
 // compute distances between pieces of curves
 List Y_in_motifs; // TODO : initialize size
 for (unsigned int i= 0; i < K; ++i){
   const uvec & curves_in_motif = curves_in_motifs[i];
   const ivec & s_k = S_clean_k[i];
   const uvec & v_dom = V_dom[i];
   for(unsigned int j = 0; j < curves_in_motif.n_elem; ++j){
     const List & y = Y[curves_in_motif[j]];
     const int s = s_k[j];
     List y_temp = List::create(_["y0"] = R_NilValue,_["y1"] = R_NilValue);
     const ivec & index = regspace<ivec>(1,v_dom.n_elem - std::max(0,1-s)) + std::max(1,s) - 1;
     int index_size = index.n_elem;
     if (use0){
       const mat & y0 = as<mat>(y["y0"]);
       const unsigned int d = y0.n_cols;
       const unsigned int y_len = y0.n_rows;
       mat y0_temp(std::max(0,1-s) + index_size,d);
       y0_temp.fill(datum::nan);
       auto filtered_j = std::views::iota(0,index_size) 
         | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
       for(int j : filtered_j)
         y0_temp.row(std::max(0,1-s) + j) =  y0.row(index[j] - 1);
       // $y0_temp[!v_dom,]=NA is the reason of v_dom(j)
       y_temp["y0"] = y0_temp;
     }
     if (use1) {
       const mat & y1 = as<mat>(y["y1"]);
       const unsigned int d = y1.n_cols;
       const unsigned int y_len = y1.n_rows;
       mat y1_temp(std::max(0,1-s) + index.n_elem,d);
       y1_temp.fill(datum::nan);
       auto filtered_j = std::views::iota(0,index_size)
         | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
       for(int j : filtered_j)
         y1_temp.row(std::max(0,1-s) + j) =  y1.row(index[j] - 1);
       // $y1_temp[!v_dom,]=NA is the reason of v_dom(j)
       y_temp["y1"] = y1_temp;
     }
     Y_in_motifs.push_back(y_temp);
   }
 }
 
 uvec Y_motifs = repLem<uvec>(regspace<uvec>(0,K-1),curves_in_motifs_number);
 const unsigned int Y_motifs_size = Y_motifs.size();
 
 unsigned int Y_in_motifs_size = Y_in_motifs.size();
 
 umat indeces_YY = combn2<uvec>(regspace<uvec>(0,Y_in_motifs_size-1)); //combn of indeces of YY
 
 ivec V_length_Y_motifs = repLem<ivec>(V_length,curves_in_motifs_number);
 
 const ivec& c = probKMA_results["c"];
 
 ivec c_Y_motifs = repLem<ivec>(c,curves_in_motifs_number);
 
 imat YY_lengths = combn2<ivec>(V_length_Y_motifs);
 
 const int YY_length_size = YY_lengths.n_cols;
 const ivec & YY_length_row0 = YY_lengths.row(0).t();
 const ivec & YY_length_row1 = YY_lengths.row(1).t();
 uvec swap = (YY_length_row0 < YY_length_row1); // const & to be add
 const uvec & equal_length = (YY_length_row0 == YY_length_row1);

 auto filtered_j_swap = std::views::iota(0,YY_length_size) 
   | std::views::filter([&swap](int j){return swap(j);});
 
 // to use std::swap o alternatives, to be checked
 for (int j : filtered_j_swap){
   std::swap(indeces_YY(0,j),indeces_YY(1,j));
 }
 
 Function find_min_diss(".find_min_diss");
 Function find_diss(".find_diss");
 
 // not checked the else condition 
 mat SD(2,YY_length_size);
 
 if(align){
   for(unsigned int i = 0; i < YY_length_size; ++i){
     const List & yy0 = Y_in_motifs[indeces_YY(0,i)];
     const List & yy1 = Y_in_motifs[indeces_YY(1,i)];
     bool equal_length_i = equal_length(i);
     NumericVector min_diss_aligned = as<NumericVector>(find_diss(yy0,    
                                                                  yy1,
                                                                  alpha,
                                                                  w,
                                                                  equal_length_i,
                                                                  d,
                                                                  use0,
                                                                  use1)); //use vec
     SD(0,i) = min_diss_aligned[0]; // set col(i) with the above vec
     SD(1,i) = min_diss_aligned[1];
   }
 }else{ //else condition to be checked
   imat c_Y_motifs_comb = combn2<ivec>(c_Y_motifs);  
   for(unsigned int i = 0; i < YY_length_size; ++i){
     const List & yy0 = Y_in_motifs[indeces_YY(0,i)];
     const List & yy1 = Y_in_motifs[indeces_YY(1,i)];
     const unsigned int cc_motifs_i = std::min(c_Y_motifs_comb(0,i),
                                               c_Y_motifs_comb(1,i));
     NumericVector min_diss = as<NumericVector>(find_min_diss(yy0,    
                                                              yy1,
                                                              alpha,
                                                              w,
                                                              cc_motifs_i,
                                                              d,
                                                              use0,
                                                              use1)); //use vec
     SD(0,i) = min_diss[0]; // set col(i) with the above vec
     SD(1,i) = min_diss[1];
   }
  
 }
 
 mat YY_D(Y_motifs_size,Y_motifs_size,fill::zeros);
 unsigned int k = 0;
 for (unsigned int j = 0; j < Y_motifs_size; ++j){
   for (unsigned int i = j+1; i <  Y_motifs_size; ++i){
     YY_D(i,j) = SD(1,k);
     k++;
   }
 }
 
 YY_D = YY_D + YY_D.t(); // check: there is a clever way to make symmetric in arma
   
 // compute intra-cluster distances
 umat intra(Y_motifs_size,Y_motifs_size,fill::zeros);
 for(unsigned int i = 0; i < Y_motifs_size; ++i){
   const uvec & temp = (Y_motifs == Y_motifs(i));
   intra.col(i) = temp;
 }
 intra.diag() *= 0;
 
 const vec & curves_in_motifs_number_Y_motifs = conv_to<vec>::from(curves_in_motifs_number.elem(Y_motifs));
 const vec & a = sum(intra%YY_D, 0).t()/(curves_in_motifs_number_Y_motifs - 1);

 // compute inter-cluster distances
 Y_motifs += 1;
 uvec Y_motifs_mod(Y_motifs_size);
 for(unsigned int i=0; i < Y_motifs_size; ++i)
     Y_motifs_mod(i) = Y_motifs(i)+1>K ? (Y_motifs(i)+1)%K - 1 : Y_motifs(i);
 const vec & curves_in_motifs_number_rep = conv_to<vec>::from(curves_in_motifs_number.elem(Y_motifs_mod));
 
 umat inter(Y_motifs_size,Y_motifs_size,fill::zeros);
 mat b_k(K-1,Y_motifs_size,fill::zeros);
 for(unsigned int k = 1; k <= K-1; ++k){
   for(unsigned int i = 0; i < Y_motifs_size; ++i){
    const int motif = Y_motifs(i);
    const unsigned int inter_motif = (motif+k)>K? (motif+k)%K : motif+k; 
    inter.col(i) = (Y_motifs== inter_motif);
   }
   b_k.row(k-1) = sum(inter%YY_D,0)/(curves_in_motifs_number_rep.t());
 }
 
 // to be tested because in my test b_k was already a matrix
 vec b(Y_motifs_size,1);
 if(b_k.n_rows > 1){
   b = min(b_k,1);
 }else{
   b = b_k.t();
 }
 
 //  compute silhouette
 vec silhouette= (b-a)/arma::max(a,b);
 for(auto & sil : silhouette){ 
   if(!is_finite(sil))
     sil = 0;
 }
 
 // compute average silhouette per cluster
 vec silhouette_average(K);
 silhouette_average.fill(datum::nan);
 for (int k = 0; k < K; k++) { 
   const uvec & indices = find(Y_motifs == k + 1);
   const vec & silhouette_k = silhouette.elem(indices);
   const uvec & sorted_indeces = sort_index(silhouette_k, "descend");
   uvec curves_in_motifs_k = curves_in_motifs[k];
   curves_in_motifs_k = curves_in_motifs_k.elem(sorted_indeces);
   curves_in_motifs_k += 1;
   curves_in_motifs[k] = wrap(curves_in_motifs_k);
   silhouette.elem(indices) = arma::sort(silhouette_k, "descend");
   silhouette_average(k) = mean(silhouette_k);
 }
 
 return List::create(silhouette,Y_motifs,curves_in_motifs,silhouette_average);
 
}
// cambiare ancora qualche lista, risolvere warning per i tipi 
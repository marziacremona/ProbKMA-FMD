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

IntegerVector myseq(int first, int last);

// [[Rcpp::export]]
NumericVector find_diss_aligned_rcpp(const List &y,
                                     const List &v,  
                                     const vec & w, 
                                     double alpha,
                                     bool aligned,
                                     unsigned int d,
                                     bool use0,
                                     bool use1)
{
  Function domain(".domain");
  Function select_domain(".select_domain");
  Function diss_d0_d1_L2(".diss_d0_d1_L2");
  
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
  int v_len = v_dom.size();
  int y_len = as<mat>(y[0]).n_rows;
  IntegerVector s_rep;
  if (aligned){
    s_rep = 1;
  } else {
    s_rep = myseq(1, y_len - v_len + 1);
  }
  std::size_t s_rep_size = s_rep.size();
  List y_rep(s_rep_size);
  
  // Convert y[0] and y[1] to mat objects
  const int index_size = v_len;
  mat temp_y0, temp_y1;
  if (use0) {
    temp_y0 = as<mat>(y[0]);
  }
  if (use1) {
    temp_y1 = as<mat>(y[1]);
  }
  
  auto index_range = std::views::iota(0,index_size);
  
  for (unsigned int i = 0; i < s_rep_size; ++i) {
    IntegerVector index = s_rep[i] - 1 + myseq(1,v_len);
    List y_rep_i = List::create(Named("y0") = R_NilValue, Named("y1") = R_NilValue);
    auto j_true = index_range
      | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));});
    if (use0) {
      mat new_y0(index_size, d);
      new_y0.fill(datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y0,&temp_y0,&index](int j){new_y0.row(j) = temp_y0.row(index[j] - 1);});
      y_rep_i["y0"] = new_y0;
    }
    if (use1) {
      mat new_y1(index_size, d);
      new_y1.fill(datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y1,&temp_y1,&index](int j){new_y1.row(j) = temp_y1.row(index[j] - 1);});
      y_rep_i["y1"] = new_y1;
      }
    y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
    y_rep[i] = y_rep_i;
    }
  
  NumericVector d_rep(s_rep_size);
  
  double min_d = std::numeric_limits<double>::max();
  int min_s = 0;
  
  for (unsigned int i = 0; i < s_rep_size; i++) {
    double dist = as<double>(diss_d0_d1_L2(y_rep[i], v_new, w, alpha));
    d_rep[i] = dist;
    if (dist < min_d){
      min_d = dist;
      min_s = s_rep[i];
    }
  }
  return NumericVector::create(min_s, min_d); 
}
  
// [[Rcpp::export]]
NumericVector find_diss(const List &y,const List &v,  
                        const vec & w, double alpha, unsigned int c_k,
                        unsigned int d,bool use0,bool use1)
{
  Function domain(".domain");
  Function select_domain(".select_domain");
  Function diss_d0_d1_L2(".diss_d0_d1_L2");
  
  // Convert domain and select_domain
  LogicalVector v_dom = as<LogicalVector>(domain(v,use0));
  Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
  int v_len = v_dom.size();
  int y_len = as<mat>(y[0]).n_rows;
  IntegerVector s_rep = myseq(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
  std::size_t s_rep_size = s_rep.size();
  List y_rep(s_rep_size);
  
  // Convert y[0] and y[1] to mat objects
  const int index_size = v_len;
  mat temp_y0, temp_y1;
  if (use0) {
    temp_y0 = as<mat>(y[0]);
  }
  if (use1) {
    temp_y1 = as<mat>(y[1]);
  }
  auto index_range = std::views::iota(0,index_size);
  
  for (unsigned int i = 0; i < s_rep_size; ++i) {
    IntegerVector index = s_rep[i] - 1 + seq_len(v_len);
    List y_rep_i = List::create(Named("y0") = R_NilValue, Named("y1") = R_NilValue);
    auto j_true = index_range
      | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));});
    if (use0) {
      mat new_y0(index_size, d);
      new_y0.fill(datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y0,&temp_y0,&index](int j){new_y0.row(j) = temp_y0.row(index[j] - 1);});
      y_rep_i["y0"] = new_y0;
    }
    if (use1) {
      mat new_y1(index_size, d);
      new_y1.fill(datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&new_y1,&temp_y1,&index](int j){new_y1.row(j) = temp_y1.row(index[j] - 1);});
      y_rep_i["y1"] = new_y1;
    }
    y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
    y_rep[i] = y_rep_i;
  }
  const unsigned int y_rep_size = y_rep.size();
  IntegerVector length_inter(y_rep_size);
  
  mat temp_y;
  int i = 0;
  for (const Rcpp::List& y_rep_i : y_rep) { // TODO: risolvere warning
    if (use0) 
      temp_y = as<mat>(y_rep_i["y0"]);
    else
      temp_y = as<mat>(y_rep_i["y1"]);
    uvec non_na_indices = find_nan(temp_y.col(0)); 
    length_inter[i] = temp_y.col(0).n_elem - non_na_indices.n_elem;
    ++i;
  }
  
  LogicalVector valid = length_inter >= c_k;
  if (sum(valid) == 0) {        
    valid[length_inter == max(length_inter)] = true; 
  }
  
  s_rep = s_rep[valid];
  y_rep = y_rep[valid];
  
  const unsigned int y_rep_size_valid = y_rep.size(); // TODO: aggiornare dimensione di y_rep
  NumericVector d_rep(y_rep_size_valid);
  
  double min_d = std::numeric_limits<double>::max();
  int min_s = 0;
  
  for (unsigned int i = 0; i < y_rep_size_valid; i++) {
    double dist = as<double>(diss_d0_d1_L2(y_rep[i], v_new, w, alpha));
    d_rep[i] = dist;
    if (dist < min_d){
      min_d = dist;
      min_s = s_rep[i];
    }
  }
  return NumericVector::create(min_s, min_d); 
}

IntegerVector myseq(int first, int last) {
  IntegerVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

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
    for (ivec::value_type j = 0; j < times[i]; ++j)
      result(k++) = v(i);
  return result;
}

// [[Rcpp::export]]
List probKMA_silhouette_rcpp(const List & probKMA_results,
                             bool align){ // @TODO = set to FALSE by default
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
  //const unsigned int N = first_y0.n_rows;
  const unsigned int K = probKMA_results["K"];
  const List & motifs = use0 ? probKMA_results["V0"] : probKMA_results["V1"];
  
  std::vector<uvec> V_dom(K); 
  ivec V_length(K);
  
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
 std::vector<uvec> curves_in_motifs(P_clean.n_cols); 
 for (unsigned int k = 0; k < P_clean.n_cols; ++k) {
   const ivec & P_clean_k = P_clean.col(k);
   const uvec & P_clean_1 = find(P_clean_k == 1); // to check in case 1 doesn't exist in the col
   curves_in_motifs[k] = P_clean_1;
 }
 
 // TODO: I don't understed what it does
 // if(!is.null(ncol(curves_in_motifs))) this if condition makes no sense since curves_in_motifs is a list
 //   curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs))) to be implemented and understood
 
 const ivec & curves_in_motifs_number=sum(P_clean,0).t();
 
 // for each centroid take the shift of the curve associated to that centroid 
 const imat & S_clean = probKMA_results["S_clean"];
 std::vector<ivec> S_clean_k(S_clean.n_cols);
 for (unsigned int k=0; k < S_clean.n_cols; ++k){
  const ivec & col_S_clean = S_clean.col(k);
  S_clean_k[k] = col_S_clean.elem(curves_in_motifs[k]);
 }
 
 // compute distances between pieces of curves
 std::vector<List> Y_in_motifs(Y_size);
 unsigned int l = 0;
 for (unsigned int i= 0; i < K; ++i){ // for each centroid
   const uvec & curves_in_motif = curves_in_motifs[i];
   const uvec & v_dom = V_dom[i];
   for(unsigned int j = 0; j < curves_in_motif.n_elem; ++j){ // for each curve assigned to that centroid
     const List & y = Y[curves_in_motif[j]]; //select the shifted part of the curve
     const int s = S_clean_k[i][j];
     List y_temp = List::create(_["y0"] = R_NilValue,_["y1"] = R_NilValue);
     const uvec & index = regspace<uvec>(1,v_dom.n_elem - std::max(0,1-s)) + std::max(1,s) - 1;
     int index_size = index.n_elem;
     if (use0){
       const mat & y0 = as<mat>(y["y0"]);
       const unsigned int d = y0.n_cols;
       const unsigned int y_len = y0.n_rows;
       mat y0_temp(std::max(0,1-s) + index_size,d);
       y0_temp.fill(datum::nan);
       auto filtered_j = std::views::iota(0,index_size) 
         | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
       for(int j : filtered_j) // maybe to use submatrix
         y0_temp.row(std::max(0,1-s) + j) =  y0.row(index[j] - 1);
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
       y_temp["y1"] = y1_temp;
     }
     Y_in_motifs[l++] = y_temp;
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
 const uvec & swap = (YY_length_row0 < YY_length_row1);
 const uvec & equal_length = (YY_length_row0 == YY_length_row1);

 auto filtered_j_swap = std::views::iota(0,YY_length_size) 
   | std::views::filter([&swap](int j){return swap(j);});
 
 for (int j : filtered_j_swap){
   std::swap(indeces_YY(0,j),indeces_YY(1,j));
 }
 
 vec SD(YY_length_size);
 
 if(align){
   
   for(int i = 0; i < YY_length_size; ++i){
     
     bool equal_length_i = equal_length(i);
     
     NumericVector min_diss_aligned = find_diss_aligned_rcpp(Y_in_motifs[indeces_YY(0,i)],    
                                                             Y_in_motifs[indeces_YY(1,i)],
                                                             w,
                                                             alpha,
                                                             equal_length_i,
                                                             d,
                                                             use0,
                                                             use1); 
                                                             
     SD(i) = min_diss_aligned[1];
  }
 }else{ 
   imat c_Y_motifs_comb = combn2<ivec>(c_Y_motifs);  
   for(int i = 0; i < YY_length_size; ++i){
     const unsigned int cc_motifs_i = std::min(c_Y_motifs_comb(0,i),
                                               c_Y_motifs_comb(1,i));
     NumericVector min_diss = find_diss(Y_in_motifs[indeces_YY(0,i)],    
                                        Y_in_motifs[indeces_YY(1,i)],
                                        w,
                                        alpha,
                                        cc_motifs_i,
                                        d,
                                        use0,
                                        use1);
     
     SD(i) = min_diss[1];
   }
  
 }
 
 mat YY_D(Y_motifs_size,Y_motifs_size,fill::zeros);
 unsigned int k = 0;
 for (unsigned int j = 0; j < Y_motifs_size; ++j){
   for (unsigned int i = j+1; i < Y_motifs_size; ++i){
     YY_D(i,j) = SD(k);
     k++;
   }
 }
 
 YY_D = YY_D + YY_D.t(); 
   
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
 
 // compute silhouette
 vec silhouette= (b-a)/arma::max(a,b);
 for(auto & sil : silhouette){ 
   if(!is_finite(sil))
     sil = 0;
 }
 
 // compute average silhouette per cluster
 vec silhouette_average(K);
 silhouette_average.fill(datum::nan);
 for (unsigned int k = 0; k < K; k++) { 
   const uvec & indices = find(Y_motifs == k + 1);
   const vec & silhouette_k = silhouette.elem(indices);
   const uvec & sorted_indeces = sort_index(silhouette_k, "descend");
   curves_in_motifs[k] = curves_in_motifs[k].elem(sorted_indeces);
   curves_in_motifs[k] += 1;
   silhouette.elem(indices) = arma::sort(silhouette_k, "descend");
   silhouette_average(k) = mean(silhouette_k);
 }
 
 return List::create(silhouette,Y_motifs,curves_in_motifs,silhouette_average);
 
}

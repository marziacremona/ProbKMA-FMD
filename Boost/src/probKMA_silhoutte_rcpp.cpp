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

List combn_rcpp(const List & y, // template version to be done
                int r) {
  List result; // the size can be set using factorial, not provided by c++
  const unsigned int n = y.size();
  uvec v(n, fill::zeros); 
  std::fill(v.begin(),v.begin() + r, 1); // the first two elements are true
  do {
    List temp(r);
    unsigned int k = 0; //ho bisogno di un puntatore che scorre su temp
    for (int i = 0; i < n; ++i) { //itero su y
      if (v[i]) {
        temp[k++] = y[i];
      }
    }
    result.push_back(temp);
  } while (std::prev_permutation(v.begin(), v.end()));
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
  
  List V_dom(K); // each element will be an uvec
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
 const mat & P_clean = probKMA_results["P_clean"];
 List curves_in_motifs(P_clean.n_cols);
 for (unsigned int k = 0; k < P_clean.n_cols; ++k) {
   const vec & P_clean_k = P_clean.col(k);
   const uvec & P_clean_1 = find(P_clean_k == 1); // to check in case 1 doesn't exist in the col
   curves_in_motifs[k] = P_clean_1;
 }
 
 // if(!is.null(ncol(curves_in_motifs))) this if condition makes no sense since curves_in_motifs is a list
 //   curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs))) to be implemented and understood
 
 const vec & curves_in_motifs_number=sum(P_clean,0).t();
 //return List::create(curves_in_motifs_number);
 
 const imat & S_clean = probKMA_results["S_clean"];
 List S_clean_k(S_clean.n_cols);
 // for sure there exists a clever way to select the elements indexed 
 for (unsigned int k=0; k < S_clean.n_cols; ++k){
  const ivec & col_S_clean = S_clean.col(k);
  const uvec & curves_in_motifs_k = curves_in_motifs[k];
  const ivec & indexed_S_clean = col_S_clean.elem(curves_in_motifs_k);
  S_clean_k[k] = indexed_S_clean;
 }
 
 // compute distances between pieces of curves
 List Y_in_motifs;
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
 

 ivec Y_motifs(accu(curves_in_motifs_number));
 unsigned int k = 0;
 for (unsigned int i = 0; i < K; ++i){
   for (unsigned int j = 0; j < curves_in_motifs_number[i]; ++j){
     Y_motifs(k++) = i;
   }
 }

 List YY = combn_rcpp(Y_in_motifs,2); 
 // YY=array(unlist(YY,recursive=FALSE),dim=c(2,length(YY))) not done, not usefull?!
 
 //sorting di V_length ?!
 List V_length_Y_motifs(accu(curves_in_motifs_number)); // potrebbe essere ivec se combn_rcpp fosse template
 List c_Y_motifs(accu(curves_in_motifs_number));
 const ivec& c = probKMA_results["c"];
 k = 0;
 for (unsigned int i = 0; i < K; ++i){
   for (unsigned int j = 0; j < curves_in_motifs_number[i]; ++j){
     const unsigned int V_length_i = V_length[i];
     const unsigned int c_i = c[i];
     c_Y_motifs[k] = c_i;
     V_length_Y_motifs[k] = V_length_i;
     k++;
   }
 }
 
 List YY_lengths = combn_rcpp(V_length_Y_motifs,2);
 
 const unsigned int YY_length_size = YY_lengths.size(); 
 uvec swap(YY_length_size);
 uvec equal_length(YY_length_size);
 for (unsigned int i = 0; i < YY_length_size; ++i){
   const List & yy_len = YY_lengths[i];
   const unsigned int yy_len0 = yy_len[0];
   const unsigned int yy_len1 = yy_len[1];
   swap(i) = yy_len0 < yy_len1;
   equal_length(i) = yy_len0 == yy_len1;
 }
 
 // to be checked swap part since swap is all zeros
 const int YY_size = YY.size();
 
 auto filtered_j_swap = std::views::iota(0,YY_size) 
   | std::views::filter([&swap](int j){return swap(j);});
 
 // to use std::swap o alternatives, to be checked
 for (int j : filtered_j_swap){
   List yy = YY[j]; // check: in this way I'm modifying also YY, I'm quite sure yes (this case test is not enough)
   const List && yy0 = yy[0];
   yy[0] = yy[1];
   yy[1] = yy0;
 }
 
 Function find_min_diss(".find_min_diss");
 Function find_diss(".find_diss");
 
 // not checked the else condition 
 mat SD(2,YY_length_size);
 
 if(align){
   for(unsigned int i = 0; i < YY_length_size; ++i){
     List yy = YY[i];
     List yy0 = yy[0];
     List yy1 = yy[1];
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
   List c_Y_motifs_comb = combn_rcpp(c_Y_motifs,2);
   for(unsigned int i = 0; i < YY_length_size; ++i){
     List yy = YY[i];
     List yy0 = yy[0];
     List yy1 = yy[1];
     const List & c_Y_motifs_comb_i = c_Y_motifs_comb[i];
     unsigned int c_Y_motifs_comb_i0 = c_Y_motifs_comb_i[0];
     unsigned int c_Y_motifs_comb_i1 = c_Y_motifs_comb_i[1];
     const unsigned int cc_motifs_i = std::min(c_Y_motifs_comb_i0,
                                               c_Y_motifs_comb_i1);
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
 
 const unsigned int Y_motifs_size = accu(curves_in_motifs_number); // to be declared above, multiple invocations to accu
 mat YY_D(Y_motifs_size,Y_motifs_size,fill::zeros);
 k = 0;
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
   const int Y_motifs_i = Y_motifs[i];
   const uvec & temp = (Y_motifs == Y_motifs_i);
   intra.col(i) = temp;
 }
 intra.diag() *= 0;
 vec curves_in_motifs_number_Y_motifs(Y_motifs_size); // TODO: change name
 for (unsigned int i = 0; i < Y_motifs_size; ++i){
   curves_in_motifs_number_Y_motifs(i) = curves_in_motifs_number(Y_motifs(i));
 }
 vec a = sum(intra%YY_D, 0).t();
 a /= (curves_in_motifs_number_Y_motifs - 1);
 
 // compute inter-cluster distances
 Y_motifs += 1;
 ivec Y_motifs_condition(Y_motifs_size);
 for(unsigned int i=0; i < Y_motifs_size; ++i){
   if (Y_motifs(i)+1>K)
     Y_motifs_condition(i) = (Y_motifs(i)+1)%K;
   else
     Y_motifs_condition(i) = Y_motifs(i)+1;
 }
 
 vec curves_in_motifs_number_Y_motifs_new(Y_motifs_size); //TODO : change name
 Y_motifs_condition -= 1;
 for (unsigned int i = 0; i < Y_motifs_size; ++i){
   curves_in_motifs_number_Y_motifs_new(i) = curves_in_motifs_number(Y_motifs_condition(i));
 }
 
 umat inter(Y_motifs_size,Y_motifs_size,fill::zeros);
 mat b_k(K-1,Y_motifs_size,fill::zeros);
 for(unsigned int k = 1; k <= K-1;++k){
   for(unsigned int i = 0; i < Y_motifs_size; ++i){
    const int motif = Y_motifs[i];
    const unsigned int temp_num = (motif+k)>K? (motif+k)%K : motif+k; //TODO : change name
    inter.col(i) = (Y_motifs== temp_num);
    b_k.row(k-1) = sum(inter%YY_D,0)/(curves_in_motifs_number_Y_motifs_new.t());
   }
 }
 
 // to be tested because in my test b_k was already a matrix
 vec b(Y_motifs_size,1);
 if(b_k.n_rows > 1){
   b = min(b_k,1);
 }else{
   b = b_k.t();
 }
 // or directly b = min(b_k,1); to be checked
 
 //  compute silhouette
 vec silhouette= (b-a)/arma::max(a,b);
 for(auto & sil : silhouette){ 
   if(!is_finite(sil))
     sil = 0;
 }
 
 // compute average silhouette per cluster
 vec silhouette_average(K);
 silhouette_average.fill(datum::nan);
 Y_motifs -= 1;
 for (int k = 0; k < K; k++) { 
   uvec indices = arma::find(Y_motifs == k);
   vec silhouette_k = silhouette.elem(indices);
   uvec sorted_indeces = sort_index(silhouette_k, "descend");
   uvec curves_in_motifs_k = curves_in_motifs[k];
   curves_in_motifs_k = curves_in_motifs_k.elem(sorted_indeces);
   curves_in_motifs_k += 1;
   curves_in_motifs[k] = wrap(curves_in_motifs_k);
   silhouette.elem(indices) = arma::sort(silhouette_k, "descend");
   silhouette_average(k) = arma::mean(silhouette_k);
 }
 
 return List::create(silhouette,silhouette_average,curves_in_motifs,Y_motifs);
}

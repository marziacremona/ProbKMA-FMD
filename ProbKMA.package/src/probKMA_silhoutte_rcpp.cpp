// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "silhoutte.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// [[Rcpp::export(.probKMA_silhouette_rcpp)]]
Rcpp::List probKMA_silhouette_rcpp(const Rcpp::List & probKMA_results,
                                   const Rcpp::Function & domain,
                                   const Rcpp::Function & select_domain,
                                   const Rcpp::Function & diss_d0_d1_L2,
                                   bool align){ 
  double alpha;
  bool use0, use1;
  unsigned int Y_size;
  
  const std::string& diss = Rcpp::as<std::string>(probKMA_results["diss"]);
  
  if (diss == "d0_L2" || diss == "d0_d1_L2") {
    const Rcpp::List & Y0 = probKMA_results["Y0"];
    Y_size = Y0.size();
  } else {
    const Rcpp::List & Y1 = probKMA_results["Y1"];
    Y_size = Y1.size();
  }
  
  Rcpp::List Y(Y_size);
  // @TODO: declare as external functions and used also in probKMA
  if (diss == "d0_L2") {
    alpha = 0;
    use0 = true;
    use1 = false;
    const Rcpp::List & Y0 = probKMA_results["Y0"];
    std::transform(Y0.cbegin(),
                   Y0.cend(),
                   Y.begin(),
                   [](const arma::mat & y0) 
                   { return Rcpp::List::create(Rcpp::_["y0"] = y0,
                                           Rcpp::_["y1"] = R_NilValue); });
  } else if (diss == "d1_L2") {
    alpha = 1;
    use0 = false;
    use1 = true;
    const Rcpp::List & Y1 = probKMA_results["Y1"];
    std::transform(Y1.cbegin(),
                   Y1.cend(),
                   Y.begin(),
                   [](const arma::mat & y1) 
                   { return Rcpp::List::create(Rcpp::_["y0"] = R_NilValue,
                                           Rcpp::_["y1"] = y1); });
  } else {
    alpha = Rcpp::as<double>(probKMA_results["alpha"]);
    use0 = true;
    use1 = true;
    const Rcpp::List & Y0 = probKMA_results["Y0"];
    const Rcpp::List & Y1 = probKMA_results["Y1"];
    std::transform(Y0.cbegin(),
                   Y0.cend(),
                   Y1.cbegin(),
                   Y.begin(),
                   [](const arma::mat & y0, const arma::mat & y1) 
                   {  return Rcpp::List::create(Rcpp::_["y0"] = y0,
                                            Rcpp::_["y1"] = y1); });
  }
  
  const arma::vec & w = probKMA_results["w"];
  const Rcpp::List & first_y = Y[0];
  const arma::mat & first_y0 = use0 ? Rcpp::as<arma::mat>(first_y[0]) : as<arma::mat>(first_y[1]); 
  const unsigned int d = first_y0.n_cols;
  //const unsigned int N = first_y0.n_rows;
  const unsigned int K = probKMA_results["K"];
  const Rcpp::List & motifs = use0 ? probKMA_results["V0"] : probKMA_results["V1"];

  
  std::vector<arma::uvec> V_dom(K); 
  arma::ivec V_length(K);
  
  // V_dom for each centroid contains an uvec with 0 if all the elements of the row of the centroid are NA
  // #TODO: testare con domini bucati
  for(unsigned int i = 0; i < K; ++i){
    const arma::mat & v = motifs[i];
    arma::uvec v_dom_k(v.n_rows);
    for (unsigned int j = 0; j < v.n_rows; ++j){
      const arma::uvec & nan_v_row_j = arma::find_nan(v.row(j));
      v_dom_k[j] = (nan_v_row_j.n_elem != d); // I am assuming that curves in R^d have centroids in R^d
    }
    V_length(i) = v.n_rows;
    V_dom[i] = v_dom_k;
  }
  
 // extract pieces of curves that fall in the different motifs
 const arma::imat & P_clean = probKMA_results["P_clean"];
 std::vector<arma::uvec> curves_in_motifs(P_clean.n_cols); 
 for (unsigned int k = 0; k < P_clean.n_cols; ++k) {
   const arma::ivec & P_clean_k = P_clean.col(k);
   const arma::uvec & P_clean_1 = find(P_clean_k == 1); // to check in case 1 doesn't exist in the col
   curves_in_motifs[k] = P_clean_1;
 }
 
 // @TODO: understand what it does
 // if(!is.null(ncol(curves_in_motifs))) this if condition makes no sense since curves_in_motifs is a list
 //   curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs))) 
 
 const arma::ivec & curves_in_motifs_number=arma::sum(P_clean,0).t();

 // for each centroid take the shift of the curve associated to that centroid 
 const arma::imat & S_clean = probKMA_results["S_clean"];
 std::vector<arma::ivec> S_clean_k(S_clean.n_cols);
 for (unsigned int k=0; k < S_clean.n_cols; ++k){
  const arma::ivec & col_S_clean = S_clean.col(k);
  S_clean_k[k] = col_S_clean.elem(curves_in_motifs[k]);
 }
 
 // compute distances between pieces of curves
 std::vector<Rcpp::List> Y_in_motifs(Y_size);
 unsigned int l = 0;
 for (unsigned int i= 0; i < K; ++i){ // for each centroid
   const arma::uvec & curves_in_motif = curves_in_motifs[i];
   const arma::uvec & v_dom = V_dom[i];
   for(unsigned int j = 0; j < curves_in_motif.n_elem; ++j){ // for each curve assigned to that centroid
     const Rcpp::List & y = Y[curves_in_motif[j]]; //select the shifted part of the curve
     const int s = S_clean_k[i][j];
     Rcpp::List y_temp = Rcpp::List::create(Rcpp::_["y0"] = R_NilValue,Rcpp::_["y1"] = R_NilValue);
     const arma::uvec & index = arma::regspace<arma::uvec>(1,v_dom.n_elem - std::max(0,1-s)) + std::max(1,s) - 1;
     int index_size = index.n_elem;
     if (use0){
       const arma::mat & y0 = Rcpp::as<arma::mat>(y["y0"]);
       const unsigned int d = y0.n_cols;
       const unsigned int y_len = y0.n_rows;
       arma::mat y0_temp(std::max(0,1-s) + index_size,d);
       y0_temp.fill(arma::datum::nan);
       auto filtered_j = std::views::iota(0,index_size) 
         | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
       for(int j : filtered_j) // maybe to use .rows provided by armadillo 
         y0_temp.row(std::max(0,1-s) + j) =  y0.row(index[j] - 1);
       y_temp["y0"] = y0_temp;
     }
     if (use1) {
       const arma::mat & y1 = Rcpp::as<arma::mat>(y["y1"]);
       const unsigned int d = y1.n_cols;
       const unsigned int y_len = y1.n_rows;
       arma::mat y1_temp(std::max(0,1-s) + index.n_elem,d);
       y1_temp.fill(arma::datum::nan);
       auto filtered_j = std::views::iota(0,index_size)
         | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
       for(int j : filtered_j)
         y1_temp.row(std::max(0,1-s) + j) =  y1.row(index[j] - 1);
       y_temp["y1"] = y1_temp;
     }
     Y_in_motifs[l++] = y_temp;
   }
 }

 arma::uvec Y_motifs = util::repLem<arma::uvec>(arma::regspace<arma::uvec>(0,K-1),curves_in_motifs_number);
 const unsigned int Y_motifs_size = Y_motifs.size();
 
 unsigned int Y_in_motifs_size = Y_in_motifs.size();
 
 arma::umat indeces_YY = util::combn2<arma::uword>(arma::regspace<arma::uvec>(0,Y_in_motifs_size-1)); //combn of indeces of Y_in_motifs
 
 arma::ivec V_length_Y_motifs = util::repLem<arma::ivec>(V_length,curves_in_motifs_number);
 
 const arma::ivec& c = probKMA_results["c"];
 
 arma::ivec c_Y_motifs = util::repLem<arma::ivec>(c,curves_in_motifs_number);
 
 arma::imat YY_lengths = util::combn2<arma::sword>(V_length_Y_motifs);
 
 const int YY_length_size = YY_lengths.n_cols;
 const arma::ivec & YY_length_row0 = YY_lengths.row(0).t();
 const arma::ivec & YY_length_row1 = YY_lengths.row(1).t();
 const arma::uvec & swap = (YY_length_row0 < YY_length_row1);
 const arma::uvec & equal_length = (YY_length_row0 == YY_length_row1);

 auto filtered_j_swap = std::views::iota(0,YY_length_size) 
   | std::views::filter([&swap](int j){return swap(j);});
 
 for (int j : filtered_j_swap){
   std::swap(indeces_YY(0,j),indeces_YY(1,j));
 }
 
 arma::vec SD(YY_length_size);
 
 if(align){
   
   for(int i = 0; i < YY_length_size; ++i){
     
     bool equal_length_i = equal_length(i);
     
     Rcpp::NumericVector min_diss_aligned = find_diss_aligned_rcpp(Y_in_motifs[indeces_YY(0,i)],    
                                                                   Y_in_motifs[indeces_YY(1,i)],
                                                                   w,
                                                                   alpha,
                                                                   equal_length_i,
                                                                   d,
                                                                   use0,
                                                                   use1,
                                                                   domain,
                                                                   select_domain,
                                                                   diss_d0_d1_L2); 
                                                             
     SD(i) = min_diss_aligned[1];
  }
 }else{ 
   arma::imat c_Y_motifs_comb = util::combn2<arma::sword>(c_Y_motifs);  
   for(int i = 0; i < YY_length_size; ++i){
     const unsigned int cc_motifs_i = std::min(c_Y_motifs_comb(0,i),
                                               c_Y_motifs_comb(1,i));
     Rcpp::NumericVector min_diss = find_diss(Y_in_motifs[indeces_YY(0,i)],    
                                              Y_in_motifs[indeces_YY(1,i)],
                                              w,
                                              alpha,
                                              cc_motifs_i,
                                              d,
                                              use0,
                                              use1,
                                              domain,
                                              select_domain,
                                              diss_d0_d1_L2);
     
     SD(i) = min_diss[1];
   }
  
 }
 
 arma::mat YY_D(Y_motifs_size,Y_motifs_size,arma::fill::zeros);
 unsigned int k = 0;
 for (unsigned int j = 0; j < Y_motifs_size; ++j){
   for (unsigned int i = j+1; i < Y_motifs_size; ++i){
     YY_D(i,j) = SD(k);
     k++;
   }
 }
 
 YY_D = YY_D + YY_D.t(); 
   
 // compute intra-cluster distances
 arma::umat intra(Y_motifs_size,Y_motifs_size,arma::fill::zeros);
 for(unsigned int i = 0; i < Y_motifs_size; ++i){
   const arma::uvec & temp = (Y_motifs == Y_motifs(i));
   intra.col(i) = temp;
 }
 intra.diag() *= 0;
 
 const arma::vec & curves_in_motifs_number_Y_motifs = arma::conv_to<arma::vec>::from(curves_in_motifs_number.elem(Y_motifs));
 const arma::vec & a = arma::sum(intra%YY_D, 0).t()/(curves_in_motifs_number_Y_motifs - 1);

 // compute inter-cluster distances
 Y_motifs += 1;
 arma::uvec Y_motifs_mod(Y_motifs_size);
 for(unsigned int i=0; i < Y_motifs_size; ++i)
     Y_motifs_mod(i) = Y_motifs(i)+1>K ? (Y_motifs(i)+1)%K - 1 : Y_motifs(i);
 const arma::vec & curves_in_motifs_number_rep = arma::conv_to<arma::vec>::from(curves_in_motifs_number.elem(Y_motifs_mod));
 
 arma::umat inter(Y_motifs_size,Y_motifs_size,arma::fill::zeros);
 arma::mat b_k(K-1,Y_motifs_size,arma::fill::zeros);
 for(unsigned int k = 1; k <= K-1; ++k){
   for(unsigned int i = 0; i < Y_motifs_size; ++i){
    const int motif = Y_motifs(i);
    const unsigned int inter_motif = (motif+k)>K? (motif+k)%K : motif+k; 
    inter.col(i) = (Y_motifs== inter_motif);
   }
   b_k.row(k-1) = sum(inter%YY_D,0)/(curves_in_motifs_number_rep.t());
 }
 
 // to be tested because in my test b_k was already a matrix
 arma::vec b(Y_motifs_size,1);
 if(b_k.n_rows > 1){
   b = min(b_k,1);
 }else{
   b = b_k.t();
 }
 
 // compute silhouette
 arma::vec silhouette= (b-a)/arma::max(a,b);
 for(auto & sil : silhouette){ 
   if(!arma::is_finite(sil))
     sil = 0;
 }
 
 // compute average silhouette per cluster
 arma::vec silhouette_average(K);
 silhouette_average.fill(arma::datum::nan);
 for (unsigned int k = 0; k < K; k++) { 
   const arma::uvec & indices = find(Y_motifs == k + 1);
   const arma::vec & silhouette_k = silhouette.elem(indices);
   const arma::uvec & sorted_indeces = arma::sort_index(silhouette_k, "descend");
   curves_in_motifs[k] = curves_in_motifs[k].elem(sorted_indeces);
   curves_in_motifs[k] += 1;
   silhouette.elem(indices) = arma::sort(silhouette_k, "descend");
   silhouette_average(k) = mean(silhouette_k);
 }
 
 return Rcpp::List::create(silhouette,Y_motifs,curves_in_motifs,silhouette_average,curves_in_motifs_number);
 
}

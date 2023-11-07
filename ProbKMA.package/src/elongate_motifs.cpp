#include "elongate_motifs.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]


void elongation_rcpp(Rcpp::List & V_new, 
                     Rcpp::List & V_dom,  
                     Rcpp::List & S_k, 
                     const arma::vec & p_k, 
                     const arma::ivec& len_elong_k, 
                     const arma::uvec& keep_k,  
                     double c, 
                     const Rcpp::Function& domain,
                     const Rcpp::Function& compute_motif,
                     const Rcpp::Function & select_domain,
                     const Rcpp::Function& diss_d0_d1_L2,
                     bool use0, bool use1,
                     const arma::vec& w, 
                     double alpha, double max_gap,  
                     const Rcpp::List& Y, int m, double deltaJk_elong,
                     const unsigned int index) 
{
  if(len_elong_k.empty()) return;
  
  const Rcpp::List & v_new_k = V_new[index];
  const arma::uvec & v_dom_k = V_dom[index];
  const arma::ivec & s_k = S_k[index];
  
  // new vec with zero at the top
  arma::ivec len_elong_k_zero(len_elong_k.size() + 1, arma::fill::zeros);
  std::copy(len_elong_k.begin(), len_elong_k.end(), len_elong_k_zero.begin() + 1);
  
  // create a matrix whose column_i contains the vector s_k - len_elong_k_zero[i]
  arma::uword len_elong_k_zero_size = len_elong_k_zero.size();
  arma::imat s_k_elong_left_right_temp(s_k.n_elem, len_elong_k_zero_size);
  
  for (arma::uword i=0; i < len_elong_k_zero_size;++i) {
    s_k_elong_left_right_temp.col(i) = s_k - len_elong_k_zero(i);
  }
  
  // create a sequence of integer from len_elong_k_zero.size() to 1
  arma::ivec reversedSequence = arma::regspace<arma::ivec>(len_elong_k_zero_size,-1,1);
  reversedSequence(0) -= 1;
  
  // repeat each col of s_k_elong_left_right a number of times specified by reversedSequence and return a list 
  std::vector<arma::ivec> s_k_elong_left_right = util::repeat_elements(s_k_elong_left_right_temp, reversedSequence);
  
  std::vector<arma::ivec> len_elong_k_right_list(len_elong_k_zero_size);
  const int max_len_elong_k = len_elong_k.back();
  unsigned int v_dom_elong_size = 0;
  //#pragma omp parallel for
  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i) {
    len_elong_k_right_list[i] = len_elong_k_zero.elem(find(len_elong_k_zero <=  max_len_elong_k - len_elong_k_zero(i)));
    v_dom_elong_size += len_elong_k_right_list[i].size();
  }
  
  //  v_dom_elong_left_right will be a vector of arma::uvec containing all the elongated version of v_dom_k
  std::vector<arma::uvec> v_dom_elong_left_right(v_dom_elong_size);  
  const int v_dom_k_len = v_dom_k.n_elem;
  unsigned int k = 0;
  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i){
    const arma::ivec & leng_elong_right_vector = len_elong_k_right_list[i];
    const unsigned int leng_elong_left = len_elong_k_zero(i);
    for (unsigned int j = 0; j <  leng_elong_right_vector.size(); ++j) {
        arma::uvec temp(leng_elong_left + v_dom_k_len + leng_elong_right_vector(j), arma::fill::ones);
        temp.rows(leng_elong_left, leng_elong_left + v_dom_k.n_elem - 1) = v_dom_k;
        v_dom_elong_left_right[k++] = temp;
    }
  }
  
  // create the list containing all the possible v_dom_k elongated using compute_motif
  const int v_elong_left_right_size = s_k_elong_left_right.size();
  Rcpp::List v_elong_left_right(v_elong_left_right_size); 
  for (int i = 0; i < v_elong_left_right_size; i++) {
    v_elong_left_right[i] = Rcpp::as<Rcpp::List>(compute_motif(Rcpp::as<Rcpp::LogicalVector>(Rcpp::wrap(v_dom_elong_left_right[i + 1])), s_k_elong_left_right[i], p_k, Y, m, use0, use1));
  }
  
  // create a LogicalVector start_with_NA whose elements are true iff the correspodent elemet of v_elong_left_right has lenght >2
  Rcpp::LogicalVector start_with_NA(v_elong_left_right_size);
  auto check_length = [](const Rcpp::List & v_elong_left_right_i){return(v_elong_left_right_i.size()>2);};
  std::transform(v_elong_left_right.begin(),v_elong_left_right.end(),start_with_NA.begin(),check_length);
  
  // filter centroid and shifts that are not in NA positions
  auto not_NA_index = std::views::iota(0,v_elong_left_right_size) 
    | std::views::filter([&start_with_NA](int index_j){return(!start_with_NA[index_j]);});
  auto filtered_elong = not_NA_index 
    | std::views::transform([&v_elong_left_right](int j){return v_elong_left_right[j];});
  auto filtered_s_k = not_NA_index 
    | std::views::transform([&s_k_elong_left_right](int j){return s_k_elong_left_right[j];});
  // I'm not filtering v_dom_elong_left_right 
  
  v_elong_left_right = Rcpp::List(filtered_elong.begin(),filtered_elong.end());
  s_k_elong_left_right = std::vector<arma::ivec>(filtered_s_k.begin(),filtered_s_k.end());
  
  // compute performance index before elongation
  double Jk_before = compute_Jk_rcpp(v_new_k, s_k, p_k, Y, alpha, w, m, use0, use1, domain, select_domain, diss_d0_d1_L2);
  
  // compute performance indexes for all possible elongations
  arma::vec c_k_after(v_elong_left_right_size);
  arma::vec Jk_after(v_elong_left_right_size);
  
  for (arma::uword i = 0; i < v_elong_left_right_size; i++) {
    const arma::uvec& domain_elong = Rcpp::as<arma::uvec>(domain(v_elong_left_right[i], use0));
    int c_i =  std::max(floor(domain_elong.n_elem*(1 - max_gap)),c); 
    c_k_after[i] = c_i; // avoid if, next line
    const Rcpp::List& v_i = v_elong_left_right[i];
    //const arma::ivec& s_i = s_k_elong_left_right[i];
    Jk_after[i] = compute_Jk_rcpp(v_i, s_k_elong_left_right[i], p_k, Y, alpha, w, m,use0 , use1,domain, select_domain, diss_d0_d1_L2, Rcpp::wrap(c_i),Rcpp::as<Rcpp::LogicalVector>(Rcpp::wrap(keep_k)));
  }
  
  // find the best elongation in terms of perf. index
  arma::vec diff_perc = ((Jk_after-Jk_before)/Jk_before);
  arma::uword best_elong = arma::index_min(diff_perc);
  
  // check that the min really exists
  bool elongate = false;
  if (best_elong < v_elong_left_right_size)
    elongate = diff_perc(best_elong) < deltaJk_elong;
  
  // evaluate if elongate or not
  if(elongate) {
    V_new[index] =  v_elong_left_right[best_elong];
    V_dom[index] =  Rcpp::as<Rcpp::LogicalVector>(Rcpp::wrap(v_dom_elong_left_right[best_elong + 1]));
    S_k[index] = s_k_elong_left_right[best_elong];
  } else {
    return;
  }
}


// function which returns V_new, V_dom, S_k after the elongation
//[[Rcpp::export(.elongate_motifs)]]
void elongate_motifs(Rcpp::List & V_new,
                     Rcpp::List & V_dom,
                     Rcpp::List & S_k,
                     const Rcpp::List & P_k,
                     const Rcpp::List & Y,
                     const arma::vec & w, 
                     int m, 
                     bool use0,
                     bool use1,
                     double alpha,
                     const arma::ivec & c,
                     const arma::ivec & c_max, // is a vector to be understood
                     double max_elong, 
                     double deltaJk_elong,
                     int trials_elong,
                     const arma::mat & D,
                     unsigned int K,
                     double max_gap,
                     const Rcpp::Function & compute_motif,
                     const Rcpp::Function & domain,
                     const Rcpp::Function & select_domain,
                     const Rcpp::Function & diss_d0_d1_L2)
{
  // with_gaps is a vector of indexes and contains all the indexes of the curves with gaps(any false) in the domain
  int V_dom_size = V_dom.size();
  Rcpp::IntegerVector len_dom(V_dom_size);
  
  auto filtered_j = std::views::iota(0,V_dom_size)
    | std::views::filter([&V_dom,&len_dom](int j){
      const arma::uvec& v_dom = V_dom[j];
      const arma::uvec& gaps_in_v_dom = arma::find(v_dom == 0);
      len_dom[j] = v_dom.size();
      return gaps_in_v_dom.n_elem > 0;});
  
  Rcpp::NumericVector with_gaps(filtered_j.begin(),filtered_j.end());
  int with_gaps_size = with_gaps.size(); 
  
  // @TODO: check this part for domains with Na
  if (with_gaps_size != 0){
    Rcpp::List V_filled(with_gaps_size);
    Rcpp::List V_dom_filled(with_gaps_size);
    arma::vec Jk_before(with_gaps_size);
    arma::vec Jk_after(with_gaps_size);
    
    // fill the domains of the motifs with gaps and recompute the motifs with the filled domains
    // and compute the perf.indexes before and after the filling
    for (unsigned int i = 0; i < with_gaps_size; ++i){
      const Rcpp::LogicalVector& v_dom = V_dom[with_gaps[i]];
      V_dom_filled[i] = Rcpp::rep(true,v_dom.size());
      V_filled[i] = Rcpp::as<Rcpp::List>(compute_motif(V_dom_filled[i],
                                           S_k[with_gaps[i]],
                                           P_k[with_gaps[i]],
                                           Y,
                                           m,
                                           use0,
                                           use1));
      Jk_before[i] = compute_Jk_rcpp(V_new[with_gaps[i]],
                                S_k[with_gaps[i]],
                                P_k[with_gaps[i]],
                                Y,
                                alpha,
                                w,
                                m,
                                use0,
                                use1,
                                domain,
                                select_domain,
                                diss_d0_d1_L2);
      Jk_after[i] = compute_Jk_rcpp(V_filled[i],
                               S_k[with_gaps[i]],
                               P_k[with_gaps[i]],
                               Y,
                               alpha,
                               w,
                               m,
                               use0,
                               use1,
                               domain,
                               select_domain,
                               diss_d0_d1_L2);
    }
    // if filling the domain improves the perf. index over a certain threshold replace the domain and the motifs with the filled one
    const arma::uvec& fill = (Jk_after-Jk_before)/Jk_before < deltaJk_elong;
    for (auto i: std::views::iota(0,with_gaps_size)
           | std::views::filter([&fill](auto j){return fill[j];})){
        V_dom[with_gaps[i]] = V_dom_filled[i];
        V_new[with_gaps[i]] = V_filled[i];
    }
  }

    Rcpp::List len_elong(V_dom_size);
    for (unsigned int i = 0; i < V_dom_size; ++i){
      const int len_max_elong_i = std::min<int>(std::floor(len_dom[i]*max_elong)
                                              ,c_max[i] - len_dom[i]);
      if (len_max_elong_i == 0){
        len_elong[i] = arma::ivec{};
      } 
      else{
        len_elong[i] =  (len_max_elong_i <= trials_elong) ? 
                        arma::regspace<arma::ivec>(1, len_max_elong_i): 
                        round(arma::linspace<arma::ivec>(1, len_max_elong_i, trials_elong));
      }
      
    }
    
    // vector of probabilities for the quantile function , to be checked this part
    arma::vec prob(1,arma::fill::value(0.25));
    // compute the quantile of the distance matrix
    arma::vec quantile = arma::quantile(vectorise(D), prob);
    // keep will be a matrix whose value i,j will be D(i,j) < quantile(0)
    arma::umat keep = D < quantile(0);
    // col-wise sum of the matrix keep
    const arma::uvec& col_sum_keep = (sum(keep, 0)).t();
    // vector of bool = true iff col_sum_keep[i]==0
    const arma::uvec& col_sum_keep_zero = (col_sum_keep==0);
    // empty_k stores the indexes of the col that have col_sum_keep = 0 
    const arma::uvec& empty_k = find(col_sum_keep_zero);
    
    for (auto k : empty_k){
      const unsigned int min_col_k_D = index_min(D.col(k));
      keep(min_col_k_D, k) = true;
    } 
    
    
    for (unsigned int i = 0; i < V_dom_size; ++i){
      
      const arma::vec& p_k_i = P_k[i];
      const arma::ivec& len_elong_k_i = len_elong[i];
      const arma::uvec& keep_k_i = keep.col(i);
      const int&  c_i = c[i];
    
      elongation_rcpp(V_new,
                      V_dom,
                      S_k,
                      p_k_i,
                      len_elong_k_i,
                      keep_k_i,
                      c_i,
                      domain,
                      compute_motif,
                      select_domain,
                      diss_d0_d1_L2,
                      use0,
                      use1,
                      w, 
                      alpha,
                      max_gap,
                      Y,
                      m,
                      deltaJk_elong,
                      i);

    }
}


#include "find_occurrences.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// [[Rcpp::export(.find_occurrences_cpp)]]
arma::mat find_occurrences_cpp(const Rcpp::List& v,
                               const Rcpp::List& Y,
                               const double R,
                               const double alpha,
                               const arma::vec& w,
                               const double c_k,
                               const bool use0,
                               const bool use1,
                               Rcpp::Function diss_d0_d1_L2,
                               Rcpp::Function domain,
                               Rcpp::Function select_domain)
{
  arma::uvec v_dom = Rcpp::as<arma::uvec>(domain(v,use0));
  std::size_t v_len = v_dom.size();
  arma::mat SD_motif(0,2);
  std::list<arma::uword> index_list;
  const Rcpp::List& temp_domain = select_domain(v,Rcpp::as<Rcpp::LogicalVector>(Rcpp::wrap(v_dom)),use0,use1);
  int t = 0;
  for(int i = 0,z = 1;i < Y.size();++i)
  {
    const Rcpp::List& Y_i = Y[i];
    std::size_t y_len = as<arma::mat>(Y_i[0]).n_rows;
    Rcpp::List y_rep(y_len-v_len+1);
    arma::ivec s_rep = arma::regspace<arma::ivec>(1,y_len-v_len+1);
    for(int j = 0;j<y_len-v_len+1;++j)
    {
      arma::mat y0;
      arma::mat y1;
      if(use0)
      {
        y0 = arma::conv_to<arma::mat>::from(Rcpp::as<arma::mat>(Y_i[0]).rows(j,v_len-1+j));
      }
      if(use1)
      {
        y1 = arma::conv_to<arma::mat>::from(as<arma::mat>(Y_i[1]).rows(j,v_len-1+j)); 
      }
      
      y_rep[j] = select_domain(Rcpp::List::create(Rcpp::Named("y0") = y0,
                                                  Rcpp::Named("y1") = y1),
                                                  Rcpp::as<Rcpp::LogicalVector>(Rcpp::wrap(v_dom)),use0,use1);
    }
    Rcpp::NumericVector valid; 
    for(int k = 0;k<y_rep.size();++k)
    {
      const arma::uvec& temp = Rcpp::as<arma::uvec>(domain(y_rep[k],use0));
      if(accu(temp) >= c_k) valid.push_back(k);
    }
    
    s_rep = arma::conv_to<arma::ivec>::from(s_rep.elem(Rcpp::as<arma::uvec>(Rcpp::wrap(valid))));
    y_rep = y_rep[valid];
    arma::vec d_rep(y_rep.size());
    for(int k = 0;k<y_rep.size();++k)
      d_rep[k] = Rcpp::as<double>(diss_d0_d1_L2(y_rep[k],temp_domain,w,alpha));
    
    arma::uvec d_rep_R(d_rep.size()+2,arma::fill::zeros);
    const arma::uvec logic = d_rep <= R;
    std::copy(logic.begin(),logic.end(),d_rep_R.begin()+1);
    arma::ivec diff_d_rep_R = arma::conv_to<arma::ivec>::from(arma::diff(d_rep_R));
    arma::uvec start = arma::find(diff_d_rep_R==1);
    arma::uvec end = arma::find(diff_d_rep_R == -1) - 1;
    
    t = 0; 
    arma::mat local_SD_motif(std::min(start.size(),end.size()),2);
    for(int k = 0;k<std::min(start.size(),end.size());++k) // NON SO SE SERVE STD::min NEL CASO GENERALE
    {
      const arma::vec& temp_d_rep = arma::conv_to<arma::vec>::from(d_rep(arma::span(start[k],end[k]))); 
      arma::uword index = arma::index_min(temp_d_rep);
      index = index + start[k];
      if(!start.empty() && !end.empty())
      {
        index_list.push_back(z);
        local_SD_motif(t,0) = s_rep[index];
        local_SD_motif(t++,1) = d_rep[index];
      }
    }
    local_SD_motif = arma::resize(local_SD_motif,t,2);
    SD_motif = arma::join_vert(SD_motif,local_SD_motif);
    ++z;
  }
  if(!SD_motif.empty())
  {
    arma::mat B(index_list.size(),1);
    int i = 0;
    for(auto it = index_list.begin();it != index_list.end();++it)
      B(i++,0) = *it;
    SD_motif = arma::join_horiz(B,SD_motif);
    return SD_motif; 
  }
  else
  {
    return arma::mat();
  }
}

// NEL CODICE R ANDRANNO RINOMINATE LE COLONNE



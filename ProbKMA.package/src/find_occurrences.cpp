#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <list>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


// [[Rcpp::export]]
arma::mat find_occurrences(const Rcpp::List& v,
                           const Rcpp::List& Y,
                           const double R,
                           const double alpha,
                           const double w,
                           const double c_k,
                           const bool use0,
                           const bool use1)
{
  Function diss_d0_d1_L2(".diss_d0_d1_L2");
  Function domain(".domain");
  Function select_domain(".select_domain");
  uvec v_dom = as<uvec>(domain(v,use0));
  std::size_t v_len = v_dom.size();
  arma::mat SD_motif(Y.size(),2);
  std::list<uword> index_list;
  const List& temp_domain = select_domain(v,Rcpp::as<Rcpp::LogicalVector>(wrap(v_dom)),use0,use1);
  int t = 0;
  for(int i = 0,z = 1;i < Y.size();++i)
  {
    const List& Y_i = Y[i];
    std::size_t y_len = as<mat>(Y_i[0]).n_rows;
    List y_rep(y_len-v_len+1);
    ivec s_rep = regspace<ivec>(1,y_len-v_len+1);
    for(int j = 0;j<y_len-v_len+1;++j)
    {
      arma::mat y0;
      arma::mat y1;
      if(use0)
        y0 = arma::conv_to<mat>::from(as<mat>(Y_i[0]).rows(j,v_len-1+j));
      if(use1)
        y1 = arma::conv_to<mat>::from(as<mat>(Y_i[1]).rows(j,v_len-1+j));
      
      y_rep[j] = select_domain(List::create(Named("y0") = y0,Named("y1") = y1),Rcpp::as<Rcpp::LogicalVector>(wrap(v_dom)),use0,use1);
    }
    NumericVector valid; // TODO: da cambiare
    for(int k = 0;k<y_rep.size();++k)
    {
      const uvec& temp = as<uvec>(domain(y_rep[k],use0));
      if(accu(temp) >= c_k) valid.push_back(k);
    }
    s_rep = arma::conv_to<ivec>::from(s_rep.elem(as<uvec>(wrap(valid))));
    y_rep = y_rep[valid];
    vec d_rep(y_rep.size());
    for(int k = 0;k<y_rep.size();++k)
      d_rep[k] = as<double>(diss_d0_d1_L2(y_rep[k],temp_domain,w,alpha));
    
    uvec d_rep_R(d_rep.size()+2,fill::zeros);
    const uvec logic = d_rep <= R;
    std::copy(logic.begin(),logic.end(),d_rep_R.begin()+1);
    ivec diff_d_rep_R = conv_to<ivec>::from(diff(d_rep_R));
    uvec start = find(diff_d_rep_R==1);
    uvec end = find(diff_d_rep_R == -1) - 1;
    for(int k = 0;k<std::min(start.size(),end.size());++k) // NON SO SE SERVE STD::min NEL CASO GENERALE
    {
      const vec& temp_d_rep = arma::conv_to<vec>::from(d_rep(span(start[k],end[k]))); 
      uword index = index_min(temp_d_rep);
      index = index + start[k];
      if(!start.empty() && !end.empty())
      {
        index_list.push_back(z);
        SD_motif(t,0) = s_rep[index];
        SD_motif(t++,1) = d_rep[index]; 
      }
    }
    ++z;
  }
  if(!SD_motif.empty())
  {
    SD_motif = resize(SD_motif,t,2);
    mat B(index_list.size(),1);
    int i = 0;
    for(auto it = index_list.begin();it != index_list.end();++it)
      B(i++,0) = *it;
    SD_motif = join_horiz(B,SD_motif);
    return SD_motif;
  }
  else
  {
    return arma::mat();
  }
}
// NEL CODICE R ANDRANNO RINOMINATE LE COLONNE



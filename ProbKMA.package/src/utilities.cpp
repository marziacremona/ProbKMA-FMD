#include "utilities.h"

Rcpp::List repeat_elements(const arma::imat& A,const arma::ivec & times) {
  arma::uword times_size = times.n_elem;
  Rcpp::List result((times_size*(times_size+1))/2 - 1);
  std::size_t i = 0;
  for(arma::uword j = 0;j < times_size;++j)
  {
    const arma::imat& B = repmat(A.col(j),1,times[j]);
    B.each_col([&result,&i](const arma::ivec& v){result[i++] = v;});
  }
  return result;
}
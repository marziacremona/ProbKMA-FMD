#include "RcppArmadillo.h"
using namespace Rcpp;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// Abstract class for dissimilarities
class Dissimilarity
{
public:
  
  Dissimilarity(): w(arma::vec()) {}; // in realtà anche arma::vec potrebbe essere di sola competenza delle distanze della prof e quindi di L2 e H1
  Dissimilarity(const arma::vec& w_): w(w_) {};
  
  // I am assuming that the data structure are vector<pair<mat,mat>>
  virtual double compute(const std::pair<arma::mat,arma::mat> &v,
                         const std::pair<arma::mat,arma::mat> &y)=0; // to be declared as const

protected:
  arma::vec w;
};

// L2 distance (change name)
class L2: public Dissimilarity
{
public:
  L2():Dissimilarity() {};  
  L2(const arma::vec& w_):Dissimilarity(w_) {};  // constructor taking arma::vec w_ as input
  
  // sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y)
  double distance(const arma::mat& v,
                  const arma::mat& y)
                  {
                    arma::mat diff = arma::square(y - v); //(y-v)^2
                    
                    diff.replace(arma::datum::nan,0);
                    
                    const arma::rowvec & col_sum = arma::sum(diff,0); //colSums
                    
                    arma::urowvec length_dom(y.n_cols,arma::fill::zeros); //length of the domain
                    
                    for(arma::uword i = 0; i < y.n_rows; ++i)
                      if(is_finite(y.row(i)))
                        length_dom += 1;
                  
                    return sum(col_sum/length_dom%w)/y.n_cols;
                  };
  
  virtual double compute(const std::pair<arma::mat,arma::mat> &v,
                         const std::pair<arma::mat,arma::mat> &y) override 
                         {
                         return this->distance(v.first,y.first);
                         }
};

class H1 final: public L2 //change name
{
public:
  
  H1():L2(), alpha(0.0) {};  
  H1(const arma::vec& w_, double alpha_):L2(w_), alpha(alpha_) {};  // constructor taking arma::vec w_ and double alpha_ as input
  
  // Override the compute method
  virtual double compute(const std::pair<arma::mat,arma::mat> &v,
                         const std::pair<arma::mat,arma::mat> &y) override
                         {
                          return alpha == 1? distance(v.second,y.second) : 
                                 (1-alpha)*distance(v.first,y.first) + 
                                 alpha*distance(v.second,y.second);
                         }
  
protected:
  
  double alpha;  
};


// [[Rcpp::export]]
double test_dissimilarities_class(const List & v,
                                  const List & y,
                                  const arma::vec &w,
                                  double alpha){
  arma::mat y0 = as<arma::mat>(y[0]);
  arma::mat y1 = as<arma::mat>(y[1]);
  arma::mat v0 = as<arma::mat>(v[0]);
  arma::mat v1 = as<arma::mat>(v[1]);
  std::pair<arma::mat,arma::mat> y_pair(std::make_pair<arma::mat,arma::mat>(std::move(y0),std::move(y1)));
  std::pair<arma::mat,arma::mat> v_pair(std::make_pair<arma::mat,arma::mat>(std::move(v0),std::move(v1)));
  H1 dist(w,alpha);
  return dist.compute(v_pair,y_pair);
}
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "RcppArmadillo.h"
#include <vector>
#include <ranges>
#include <algorithm>

namespace util
{
  /// Repeat columns inside a matrix
  //** 
  //* @param[A] matrix whose columns have to be repeated.
  //* @param[times] vector that contain how many times each column vector has 
  //*               to be repeated
  //* 
  //* @return a vector of armadillo vectors that contain the repeated columns 
  std::vector<arma::ivec> repeat_elements(const arma::imat& A,const arma::ivec & times);

  /// Create a uniform sequence starting from first to last included.
  //**
  //* @param[first] start of the sequence 
  //* @param[last] end of the sequence 
  //* 
  //* @return a sequence of number spanned by one from first to last 
  Rcpp::IntegerVector myseq(int first, int last);
  
  
  template<typename T,typename Comparator>
  arma::vec avg_rank(const T& x)
  {
    R_xlen_t sz = x.size();
    Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
    std::sort(w.begin(), w.end(), Comparator(x));
    
    arma::vec r;
    r.set_size(sz);
    for (R_xlen_t n, i = 0; i < sz; i += n) {
      n = 1;
      while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
      for (R_xlen_t k = 0; k < n; k++) {
        r[w[i + k]] = i + (n + 1) / 2.;
      }
    }
    return r;
  }
  
  template<typename V,typename T>
  T order2(const V& x, bool desc = false) {
    std::size_t n = x.size();
    T idx;
    if constexpr(std::is_same<V, arma::uvec>::value)
    {
      idx.set_size(n); 
    }
    else
    {
      idx = T(n);
    }
    std::iota(idx.begin(), idx.end(), static_cast<size_t>(1));
    if (desc) {
      auto comparator = [&x](size_t a, size_t b){ return x[a - 1] > x[b - 1]; };
      std::stable_sort(idx.begin(), idx.end(), comparator);
    } else {
      auto comparator = [&x](size_t a, size_t b){ return x[a - 1] < x[b - 1]; };
      std::stable_sort(idx.begin(), idx.end(), comparator);
      // simulate na.last
      size_t nas = 0;
      for (size_t i = 0; i < n; ++i, ++nas)
        if (!Rcpp::Vector<REALSXP>::is_na(x[idx[i] - 1])) break;
        std::rotate(idx.begin(), idx.begin() + nas, idx.end());
    }
    return idx;
  }
  
  
  template<typename T>
  arma::Mat<T> combn2(const arma::Col<T> & y){  // Col<double> = vec, Col<uword> = uvec, Col<sword> = ivec
    int n = y.n_elem;
    arma::uvec v(n,arma::fill::zeros);  
    v(0) = 1;
    v(1) = 1;
    arma::uword l = 0;
    arma::Mat<T> result(2,n*(n-1)/2, arma::fill::zeros); 
    arma::uword k;
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
  
  template<class T>
  T repLem(const T & v,
           const arma::ivec & times){
    T result(arma::accu(times));
    std::size_t k = 0;
    std::size_t times_size = times.size();
    // @TODO: check times_size == v.size()
    for (std::size_t i = 0; i < times_size; ++i) 
      for (arma::ivec::value_type j = 0; j < times[i]; ++j)
        result(k++) = v(i);
    return result;
  }
  
}



#endif //UTILITIES_HPP